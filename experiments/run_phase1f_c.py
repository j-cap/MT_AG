from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import yaml

from mt_ag.aacopf import run_literal_small_aacopf
from mt_ag.geometry import wrap_angle
from mt_ag.paper_pf import generate_paper_trajectory
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory

ROOT = Path(__file__).resolve().parents[1]


def terminal_hold_time(mask, dt, hold_steps):
    mask = np.asarray(mask, dtype=bool)
    false_indices = np.flatnonzero(~mask)
    start = 0 if len(false_indices) == 0 else int(false_indices[-1] + 1)
    if len(mask) - start < hold_steps:
        return None
    return float(start * dt)


def distribution(values):
    x = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "std": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def controlled_particles(trajectory, auxiliary, cloud_cfg, correct_fraction, rng):
    n_particles = cloud_cfg["n_particles"]
    n_correct = round(correct_fraction * n_particles)
    n_wrong = n_particles - n_correct
    position_std = cloud_cfg["position_std_m"]
    yaw_std = np.deg2rad(cloud_cfg["yaw_std_deg"])

    truth_position = trajectory.state[0, :2]
    truth_yaw = trajectory.state[0, 2]
    vector = truth_position - auxiliary[0]
    range0 = np.linalg.norm(vector)
    truth_bearing = np.arctan2(vector[1], vector[0])
    wrong_bearing = truth_bearing + np.deg2rad(cloud_cfg["wrong_bearing_offset_deg"])
    wrong_center = auxiliary[0] + range0 * np.array(
        [np.cos(wrong_bearing), np.sin(wrong_bearing)]
    )
    wrong_yaw = wrap_angle(truth_yaw + np.deg2rad(cloud_cfg["wrong_yaw_offset_deg"]))

    correct = np.column_stack(
        [
            rng.normal(truth_position[0], position_std, n_correct),
            rng.normal(truth_position[1], position_std, n_correct),
            wrap_angle(rng.normal(truth_yaw, yaw_std, n_correct)),
        ]
    )
    wrong = np.column_stack(
        [
            rng.normal(wrong_center[0], position_std, n_wrong),
            rng.normal(wrong_center[1], position_std, n_wrong),
            wrap_angle(rng.normal(wrong_yaw, yaw_std, n_wrong)),
        ]
    )
    particles = np.vstack([correct, wrong])
    rng.shuffle(particles, axis=0)
    return particles


def build_datasets(cfg, trajectory, auxiliary, split_name, quick):
    split_cfg = cfg[split_name]
    cloud_cfg = cfg["particle_cloud"].copy()
    if quick:
        cloud_cfg["n_particles"] = min(120, cloud_cfg["n_particles"])
    n_seeds = min(2, split_cfg["n_seeds"]) if quick else split_cfg["n_seeds"]
    datasets = []
    for scenario_index, (scenario_name, scenario_cfg) in enumerate(cfg["scenarios"].items()):
        for offset in range(n_seeds):
            seed = split_cfg["seed_start"] + offset
            range_rng = np.random.default_rng(
                cfg["randomness"]["range_seed_offset"] + 10000 * scenario_index + seed
            )
            ranges = generate_uwb_ranges(
                trajectory.state[:, :2],
                auxiliary,
                range_rng,
                cfg["uwb"]["sigma_range_m"],
            )
            init_rng = np.random.default_rng(
                cfg["randomness"]["initializer_seed_offset"] + 10000 * scenario_index + seed
            )
            particles = controlled_particles(
                trajectory,
                auxiliary,
                cloud_cfg,
                scenario_cfg["correct_fraction"],
                init_rng,
            )
            datasets.append(
                {
                    "scenario": scenario_name,
                    "seed": seed,
                    "ranges": ranges,
                    "particles": particles,
                }
            )
    return datasets, cloud_cfg


def evaluate_dataset(dataset, setting, cfg, trajectory, auxiliary, cloud_cfg):
    yaw_threshold = np.deg2rad(cloud_cfg["correct_mode_yaw_threshold_deg"])
    result = run_literal_small_aacopf(
        trajectory.increments,
        dataset["ranges"],
        auxiliary,
        dataset["particles"],
        sigma_uwb_m=cfg["uwb"]["sigma_range_m"],
        propagation_convention="pre_turn",
        alpha=setting["alpha"],
        beta=setting["beta"],
        c_lambda=setting["c_lambda"],
        epsilon_distance=cfg["aco"]["epsilon_distance"],
        epsilon_weight=cfg["aco"]["epsilon_weight"],
        initial_particles_conditioned_on_z0=True,
        truth_state=trajectory.state,
        correct_mode_position_threshold_m=cloud_cfg["correct_mode_position_threshold_m"],
        correct_mode_yaw_threshold_rad=yaw_threshold,
    )

    hold_steps = max(1, round(cloud_cfg["hold_time_s"] / trajectory.dt))
    mode_fraction = result.correct_mode_fraction_post_aco
    correct_time = terminal_hold_time(
        mode_fraction >= cloud_cfg["concentration_fraction"],
        trajectory.dt,
        hold_steps,
    )
    wrong_time = terminal_hold_time(
        mode_fraction <= cloud_cfg["wrong_lock_fraction"],
        trajectory.dt,
        hold_steps,
    )
    update_slice = slice(1, None)
    unique_parent = result.unique_parent_fraction[update_slice]
    max_multiplicity = result.max_destination_multiplicity[update_slice]
    n_particles = len(dataset["particles"])
    catastrophic = bool(
        np.min(unique_parent) < cloud_cfg["catastrophic_unique_parent_fraction"]
    )
    dominant_clone = bool(
        np.max(max_multiplicity) >= cloud_cfg["dominant_clone_fraction"] * n_particles
    )

    position_error = np.linalg.norm(result.estimate[:, :2] - trajectory.state[:, :2], axis=1)
    yaw_error = np.abs(wrap_angle(result.estimate[:, 2] - trajectory.state[:, 2]))
    return {
        "scenario": dataset["scenario"],
        "seed": dataset["seed"],
        "correct_lock": correct_time is not None,
        "correct_lock_time_s": correct_time,
        "wrong_lock": wrong_time is not None,
        "wrong_lock_time_s": wrong_time,
        "catastrophic_ancestry_collapse": catastrophic,
        "dominant_clone_event": dominant_clone,
        "final_correct_mode_fraction": float(mode_fraction[-1]),
        "final_correct_mode_mass_pre_aco": float(result.correct_mode_mass_pre_aco[-1]),
        "minimum_unique_parent_fraction": float(np.min(unique_parent)),
        "mean_unique_parent_fraction": float(np.mean(unique_parent)),
        "maximum_destination_multiplicity_fraction": float(
            np.max(max_multiplicity) / n_particles
        ),
        "mean_moved_fraction": float(np.mean(result.moved_fraction[update_slice])),
        "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
        "yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error**2)))),
        "transition_runtime_s": float(np.sum(result.transition_runtime_s[update_slice])),
    }


def aggregate_runs(runs):
    scenario_names = sorted({run["scenario"] for run in runs})

    def aggregate_subset(subset):
        successful_times = [
            run["correct_lock_time_s"]
            for run in subset
            if run["correct_lock_time_s"] is not None
        ]
        return {
            "n": len(subset),
            "correct_lock_fraction": float(np.mean([run["correct_lock"] for run in subset])),
            "wrong_lock_fraction": float(np.mean([run["wrong_lock"] for run in subset])),
            "catastrophic_collapse_fraction": float(
                np.mean([run["catastrophic_ancestry_collapse"] for run in subset])
            ),
            "dominant_clone_fraction": float(
                np.mean([run["dominant_clone_event"] for run in subset])
            ),
            "correct_lock_time_s": distribution(successful_times) if successful_times else None,
            "final_correct_mode_fraction": distribution(
                [run["final_correct_mode_fraction"] for run in subset]
            ),
            "minimum_unique_parent_fraction": distribution(
                [run["minimum_unique_parent_fraction"] for run in subset]
            ),
            "maximum_destination_multiplicity_fraction": distribution(
                [run["maximum_destination_multiplicity_fraction"] for run in subset]
            ),
            "mean_moved_fraction": distribution([run["mean_moved_fraction"] for run in subset]),
            "position_rmse_m": distribution([run["position_rmse_m"] for run in subset]),
            "yaw_rmse_deg": distribution([run["yaw_rmse_deg"] for run in subset]),
            "transition_runtime_s": distribution(
                [run["transition_runtime_s"] for run in subset]
            ),
        }

    return {
        "overall": aggregate_subset(runs),
        "by_scenario": {
            name: aggregate_subset([run for run in runs if run["scenario"] == name])
            for name in scenario_names
        },
    }


def setting_rank_key(entry, success_floor):
    metrics = entry["aggregate"]["overall"]
    eligible = metrics["correct_lock_fraction"] >= success_floor
    if eligible:
        return (
            0,
            metrics["wrong_lock_fraction"],
            metrics["catastrophic_collapse_fraction"],
            metrics["dominant_clone_fraction"],
            -metrics["correct_lock_fraction"],
            -metrics["final_correct_mode_fraction"]["mean"],
            -metrics["minimum_unique_parent_fraction"]["mean"],
        )
    return (
        1,
        -metrics["correct_lock_fraction"],
        metrics["wrong_lock_fraction"],
        metrics["catastrophic_collapse_fraction"],
        metrics["dominant_clone_fraction"],
        -metrics["final_correct_mode_fraction"]["mean"],
        -metrics["minimum_unique_parent_fraction"]["mean"],
    )


def compact_setting(entry):
    overall = entry["aggregate"]["overall"]
    return {
        **entry["setting"],
        "correct_lock_fraction": overall["correct_lock_fraction"],
        "wrong_lock_fraction": overall["wrong_lock_fraction"],
        "catastrophic_collapse_fraction": overall["catastrophic_collapse_fraction"],
        "dominant_clone_fraction": overall["dominant_clone_fraction"],
        "final_correct_mode_fraction_mean": overall["final_correct_mode_fraction"]["mean"],
        "minimum_unique_parent_fraction_mean": overall[
            "minimum_unique_parent_fraction"
        ]["mean"],
        "mean_moved_fraction": overall["mean_moved_fraction"]["mean"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1f_c.yaml").read_text())
    sim_cfg = cfg["simulation"].copy()
    if args.quick:
        sim_cfg["duration"] = min(3.0, sim_cfg["duration"])
    trajectory = generate_paper_trajectory(**sim_cfg)
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")

    development, dev_cloud_cfg = build_datasets(
        cfg, trajectory, auxiliary, "development", args.quick
    )
    alpha_values = cfg["sweep"]["alpha"]
    beta_values = cfg["sweep"]["beta"]
    c_values = cfg["sweep"]["c_lambda"]
    if args.quick:
        alpha_values = alpha_values[:2]
        beta_values = beta_values[:2]
        c_values = c_values[:3]

    entries = []
    for alpha, beta, c_lambda in itertools.product(alpha_values, beta_values, c_values):
        setting = {
            "alpha": float(alpha),
            "beta": float(beta),
            "c_lambda": float(c_lambda),
        }
        runs = [
            evaluate_dataset(dataset, setting, cfg, trajectory, auxiliary, dev_cloud_cfg)
            for dataset in development
        ]
        entries.append(
            {
                "setting": setting,
                "aggregate": aggregate_runs(runs),
                "runs": runs,
            }
        )

    success_floor = cfg["selection"]["minimum_development_correct_lock_fraction"]
    ranked = sorted(entries, key=lambda entry: setting_rank_key(entry, success_floor))
    selected = ranked[0]

    validation, val_cloud_cfg = build_datasets(
        cfg, trajectory, auxiliary, "validation", args.quick
    )
    validation_runs = [
        evaluate_dataset(
            dataset,
            selected["setting"],
            cfg,
            trajectory,
            auxiliary,
            val_cloud_cfg,
        )
        for dataset in validation
    ]
    validation_aggregate = aggregate_runs(validation_runs)
    val_overall = validation_aggregate["overall"]
    accepted = bool(
        val_overall["correct_lock_fraction"]
        >= cfg["selection"]["validation_correct_lock_target"]
        and val_overall["wrong_lock_fraction"]
        <= cfg["selection"]["validation_wrong_lock_maximum"]
        and val_overall["catastrophic_collapse_fraction"]
        <= cfg["selection"]["validation_catastrophic_collapse_maximum"]
    )

    raw = {
        "design": {
            "purpose": "P1F-C safety-aware alpha/beta/c_lambda tuning",
            "development_cases": list(cfg["scenarios"]),
            "development_run_count_per_setting": len(development),
            "validation_run_count": len(validation),
            "selection_policy": cfg["selection"]["policy"],
            "quick": args.quick,
        },
        "development_ranked": ranked,
        "selected_setting": selected["setting"],
        "validation": {
            "aggregate": validation_aggregate,
            "runs": validation_runs,
        },
        "freeze_accepted": accepted,
    }

    summary = {
        "selected_setting": selected["setting"],
        "development": compact_setting(selected),
        "validation": {
            "overall": {
                "correct_lock_fraction": val_overall["correct_lock_fraction"],
                "wrong_lock_fraction": val_overall["wrong_lock_fraction"],
                "catastrophic_collapse_fraction": val_overall[
                    "catastrophic_collapse_fraction"
                ],
                "dominant_clone_fraction": val_overall["dominant_clone_fraction"],
                "final_correct_mode_fraction_mean": val_overall[
                    "final_correct_mode_fraction"
                ]["mean"],
                "minimum_unique_parent_fraction_mean": val_overall[
                    "minimum_unique_parent_fraction"
                ]["mean"],
                "mean_moved_fraction": val_overall["mean_moved_fraction"]["mean"],
                "position_rmse_mean_m": val_overall["position_rmse_m"]["mean"],
                "yaw_rmse_mean_deg": val_overall["yaw_rmse_deg"]["mean"],
            },
            "by_scenario": {
                name: {
                    "correct_lock_fraction": values["correct_lock_fraction"],
                    "wrong_lock_fraction": values["wrong_lock_fraction"],
                    "catastrophic_collapse_fraction": values[
                        "catastrophic_collapse_fraction"
                    ],
                    "final_correct_mode_fraction_mean": values[
                        "final_correct_mode_fraction"
                    ]["mean"],
                }
                for name, values in validation_aggregate["by_scenario"].items()
            },
        },
        "freeze_accepted": accepted,
        "top_10_development": [compact_setting(entry) for entry in ranked[:10]],
    }

    out_dir = ROOT / "results/phase1/p1f_c"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
