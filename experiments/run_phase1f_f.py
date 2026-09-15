from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml

from mt_ag.aacopf import run_literal_small_aacopf
from mt_ag.geometry import wrap_angle
from mt_ag.paper_pf import generate_paper_trajectory, run_paper_bootstrap_pf
from mt_ag.simulation import auxiliary_trajectory

ROOT = Path(__file__).resolve().parents[1]


def distribution(values):
    x = np.asarray(values, dtype=float)
    if len(x) == 0:
        return None
    return {
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "std": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def terminal_hold_time(mask, dt, hold_steps):
    mask = np.asarray(mask, dtype=bool)
    false_indices = np.flatnonzero(~mask)
    start = 0 if len(false_indices) == 0 else int(false_indices[-1] + 1)
    if len(mask) - start < hold_steps:
        return None
    return float(start * dt)


def recovery_time_after_index(signal, threshold, start_index, dt, hold_steps):
    if start_index is None:
        return None
    signal = np.asarray(signal, dtype=float)
    for index in range(start_index, len(signal) - hold_steps + 1):
        if np.all(signal[index : index + hold_steps] >= threshold):
            return float((index - start_index) * dt)
    return None


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


def build_ranges(cfg, condition_name, trajectory, auxiliary, seed):
    sigma = cfg["uwb"]["nominal_sigma_range_m"]
    ideal = np.linalg.norm(trajectory.state[:, :2] - auxiliary, axis=1)
    base_rng = np.random.default_rng(cfg["randomness"]["nominal_range_seed_offset"] + seed)
    base_noise = base_rng.normal(0.0, sigma, len(ideal))
    additional = np.zeros_like(ideal)
    condition = cfg["ranging_conditions"][condition_name]
    family = condition["family"]
    metadata = {"family": family, "indices": []}

    if family == "clean":
        pass
    elif family == "impulses":
        rng = np.random.default_rng(cfg["randomness"]["impulse_seed_offset"] + seed)
        eligible = np.arange(1, len(ideal))
        indices = np.sort(rng.choice(eligible, size=condition["n_events"], replace=False))
        if condition.get("symmetric_sign", False):
            signs = rng.choice(np.array([-1.0, 1.0]), size=len(indices))
        else:
            signs = np.ones(len(indices))
        additional[indices] = condition["magnitude_m"] * signs
        metadata.update(
            {
                "indices": indices.tolist(),
                "event_values_m": additional[indices].tolist(),
            }
        )
    elif family == "mixture":
        rng = np.random.default_rng(cfg["randomness"]["mixture_seed_offset"] + seed)
        mask = rng.uniform(size=len(ideal)) < condition["probability"]
        mask[0] = False
        standard_gross = rng.normal(size=len(ideal))
        additional[mask] = condition["gross_sigma_m"] * standard_gross[mask]
        metadata.update(
            {
                "indices": np.flatnonzero(mask).tolist(),
                "gross_sigma_m": condition["gross_sigma_m"],
            }
        )
    elif family == "nlos_burst":
        rng = np.random.default_rng(cfg["randomness"]["nlos_seed_offset"] + seed)
        max_duration_s = max(
            item.get("duration_s", 0.0)
            for item in cfg["ranging_conditions"].values()
            if item["family"] == "nlos_burst"
        )
        max_duration_steps = max(1, round(max_duration_s / trajectory.dt))
        earliest = max(1, round(1.0 / trajectory.dt))
        latest = len(ideal) - max_duration_steps - 1
        start = int(rng.integers(earliest, max(earliest + 1, latest + 1)))
        duration_steps = max(1, round(condition["duration_s"] / trajectory.dt))
        stop = min(len(ideal), start + duration_steps)
        indices = np.arange(start, stop)
        standard_extra = rng.normal(size=len(ideal))
        additional[indices] = condition["positive_bias_m"] + condition["extra_sigma_m"] * standard_extra[indices]
        metadata.update(
            {
                "indices": indices.tolist(),
                "start_time_s": float(start * trajectory.dt),
                "end_time_s": float((stop - 1) * trajectory.dt),
                "positive_bias_m": condition["positive_bias_m"],
            }
        )
    else:
        raise ValueError(f"Unknown ranging family: {family}")

    indices = np.asarray(metadata["indices"], dtype=int)
    metadata["n_corrupted"] = int(len(indices))
    metadata["first_corrupted_index"] = int(indices[0]) if len(indices) else None
    metadata["last_corrupted_index"] = int(indices[-1]) if len(indices) else None
    metadata["max_abs_additional_error_m"] = float(np.max(np.abs(additional)))
    return ideal + base_noise + additional, metadata


def pose_metrics(estimate, truth, cloud_cfg, dt):
    position_error = np.linalg.norm(estimate[:, :2] - truth[:, :2], axis=1)
    yaw_error = np.abs(wrap_angle(estimate[:, 2] - truth[:, 2]))
    hold_steps = max(1, round(cloud_cfg["hold_time_s"] / dt))
    pose_ok = (position_error < cloud_cfg["correct_mode_position_threshold_m"]) & (
        yaw_error < np.deg2rad(cloud_cfg["correct_mode_yaw_threshold_deg"])
    )
    return {
        "pose_success": terminal_hold_time(pose_ok, dt, hold_steps) is not None,
        "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
        "yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error**2)))),
        "final_position_error_m": float(position_error[-1]),
        "final_yaw_error_deg": float(np.rad2deg(yaw_error[-1])),
    }


def mode_metrics(signal, metadata, cfg, dt):
    cloud_cfg = cfg["particle_cloud"]
    hold_steps = max(1, round(cloud_cfg["hold_time_s"] / dt))
    recovery_hold_steps = max(1, round(cfg["recovery"]["hold_time_s"] / dt))
    signal = np.asarray(signal, dtype=float)
    correct_time = terminal_hold_time(
        signal >= cloud_cfg["concentration_fraction"], dt, hold_steps
    )
    wrong_time = terminal_hold_time(signal <= cloud_cfg["wrong_lock_fraction"], dt, hold_steps)
    first = metadata["first_corrupted_index"]
    last = metadata["last_corrupted_index"]
    after_first = signal[first:] if first is not None else signal
    return {
        "correct_lock": correct_time is not None,
        "correct_lock_time_s": correct_time,
        "wrong_lock": wrong_time is not None,
        "wrong_lock_time_s": wrong_time,
        "final_correct_mode_signal": float(signal[-1]),
        "minimum_correct_mode_signal_after_first_corruption": float(np.min(after_first)),
        "recovery_time_after_last_corruption_s": recovery_time_after_index(
            signal,
            cfg["recovery"]["threshold_fraction"],
            last,
            dt,
            recovery_hold_steps,
        ),
    }


def evaluate_one(cfg, condition_name, initial_scenario, seed):
    trajectory = generate_paper_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, cfg["geometry"]["mode"])
    ranges, corruption = build_ranges(cfg, condition_name, trajectory, auxiliary, seed)
    cloud_cfg = cfg["particle_cloud"]
    scenario_cfg = cfg["initial_mode_scenarios"][initial_scenario]
    init_rng = np.random.default_rng(
        cfg["randomness"]["initializer_seed_offset"]
        + 10000 * list(cfg["initial_mode_scenarios"]).index(initial_scenario)
        + seed
    )
    particles = controlled_particles(
        trajectory,
        auxiliary,
        cloud_cfg,
        scenario_cfg["correct_fraction"],
        init_rng,
    )
    yaw_threshold = np.deg2rad(cloud_cfg["correct_mode_yaw_threshold_deg"])

    pf_rng = np.random.default_rng(cfg["randomness"]["pf_seed_offset"] + seed)
    pf = run_paper_bootstrap_pf(
        trajectory.increments,
        ranges,
        auxiliary,
        particles.copy(),
        pf_rng,
        sigma_uwb_m=cfg["uwb"]["assumed_likelihood_sigma_m"],
        resample_fraction=cfg["conventional_pf"]["resample_fraction"],
        propagation_convention=cfg["paper_model"]["propagation_convention"],
        sigma_delta_l_m=cfg["paper_model"]["sigma_delta_l_m"],
        sigma_delta_phi_rad=cfg["paper_model"]["sigma_delta_phi_rad"],
        initial_particles_conditioned_on_z0=True,
        truth_state=trajectory.state,
        correct_mode_position_threshold_m=cloud_cfg["correct_mode_position_threshold_m"],
        correct_mode_yaw_threshold_rad=yaw_threshold,
    )
    pf_metrics = pose_metrics(pf.estimate, trajectory.state, cloud_cfg, trajectory.dt)
    pf_metrics.update(mode_metrics(pf.correct_mode_mass, corruption, cfg, trajectory.dt))
    pf_metrics.update(
        {
            "resampling_events": int(
                np.sum(pf.unique_fraction_post_transition < 1.0 - 1e-15)
            ),
            "minimum_unique_fraction": float(np.min(pf.unique_fraction_post_transition)),
        }
    )

    aco_cfg = cfg["frozen_aacopf"]
    aco = run_literal_small_aacopf(
        trajectory.increments,
        ranges,
        auxiliary,
        particles.copy(),
        sigma_uwb_m=cfg["uwb"]["assumed_likelihood_sigma_m"],
        propagation_convention=cfg["paper_model"]["propagation_convention"],
        alpha=aco_cfg["alpha"],
        beta=aco_cfg["beta"],
        c_lambda=aco_cfg["c_lambda"],
        epsilon_distance=aco_cfg["epsilon_distance"],
        epsilon_weight=aco_cfg["epsilon_weight"],
        initial_particles_conditioned_on_z0=True,
        truth_state=trajectory.state,
        correct_mode_position_threshold_m=cloud_cfg["correct_mode_position_threshold_m"],
        correct_mode_yaw_threshold_rad=yaw_threshold,
    )
    aco_metrics = pose_metrics(aco.estimate, trajectory.state, cloud_cfg, trajectory.dt)
    aco_metrics.update(
        mode_metrics(aco.correct_mode_fraction_post_aco, corruption, cfg, trajectory.dt)
    )
    update_slice = slice(1, None)
    min_parent = float(np.min(aco.unique_parent_fraction[update_slice]))
    max_multiplicity = int(np.max(aco.max_destination_multiplicity[update_slice]))
    aco_metrics.update(
        {
            "final_correct_mode_mass_pre_aco": float(aco.correct_mode_mass_pre_aco[-1]),
            "minimum_unique_parent_fraction": min_parent,
            "catastrophic_collapse": bool(
                min_parent < cloud_cfg["catastrophic_unique_parent_fraction"]
            ),
            "dominant_clone": bool(
                max_multiplicity >= cloud_cfg["dominant_clone_fraction"] * len(particles)
            ),
            "maximum_destination_multiplicity_fraction": float(
                max_multiplicity / len(particles)
            ),
            "mean_moved_fraction": float(np.mean(aco.moved_fraction[update_slice])),
            "transition_runtime_s": float(np.sum(aco.transition_runtime_s[update_slice])),
        }
    )

    return {
        "condition": condition_name,
        "initial_scenario": initial_scenario,
        "seed": seed,
        "corruption": corruption,
        "pf": pf_metrics,
        "aacopf": aco_metrics,
    }


def aggregate_method(runs, method):
    values = [run[method] for run in runs]
    recovery = [
        value["recovery_time_after_last_corruption_s"]
        for value in values
        if value["recovery_time_after_last_corruption_s"] is not None
    ]
    result = {
        "n": len(values),
        "correct_lock_fraction": float(np.mean([v["correct_lock"] for v in values])),
        "wrong_lock_fraction": float(np.mean([v["wrong_lock"] for v in values])),
        "pose_success_fraction": float(np.mean([v["pose_success"] for v in values])),
        "position_rmse_m": distribution([v["position_rmse_m"] for v in values]),
        "yaw_rmse_deg": distribution([v["yaw_rmse_deg"] for v in values]),
        "final_correct_mode_signal": distribution(
            [v["final_correct_mode_signal"] for v in values]
        ),
        "minimum_correct_mode_signal_after_first_corruption": distribution(
            [v["minimum_correct_mode_signal_after_first_corruption"] for v in values]
        ),
        "recovery_time_after_last_corruption_s": distribution(recovery),
        "recovery_fraction": float(
            np.mean([v["recovery_time_after_last_corruption_s"] is not None for v in values])
        ),
    }
    if method == "pf":
        result.update(
            {
                "resampling_events": distribution([v["resampling_events"] for v in values]),
                "minimum_unique_fraction": distribution(
                    [v["minimum_unique_fraction"] for v in values]
                ),
            }
        )
    else:
        result.update(
            {
                "catastrophic_collapse_fraction": float(
                    np.mean([v["catastrophic_collapse"] for v in values])
                ),
                "dominant_clone_fraction": float(np.mean([v["dominant_clone"] for v in values])),
                "minimum_unique_parent_fraction": distribution(
                    [v["minimum_unique_parent_fraction"] for v in values]
                ),
                "maximum_destination_multiplicity_fraction": distribution(
                    [v["maximum_destination_multiplicity_fraction"] for v in values]
                ),
                "mean_moved_fraction": distribution([v["mean_moved_fraction"] for v in values]),
            }
        )
    return result


def aggregate_runs(runs):
    pf_success = np.array([run["pf"]["correct_lock"] for run in runs], dtype=bool)
    aco_success = np.array([run["aacopf"]["correct_lock"] for run in runs], dtype=bool)
    return {
        "n": len(runs),
        "corruption_count": distribution([run["corruption"]["n_corrupted"] for run in runs]),
        "max_abs_additional_error_m": distribution(
            [run["corruption"]["max_abs_additional_error_m"] for run in runs]
        ),
        "pf": aggregate_method(runs, "pf"),
        "aacopf": aggregate_method(runs, "aacopf"),
        "paired_lock_outcomes": {
            "both_correct_lock": int(np.sum(pf_success & aco_success)),
            "aacopf_only_correct_lock": int(np.sum(~pf_success & aco_success)),
            "pf_only_correct_lock": int(np.sum(pf_success & ~aco_success)),
            "neither_correct_lock": int(np.sum(~pf_success & ~aco_success)),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--condition", required=True)
    parser.add_argument("--initial-scenario", required=True)
    parser.add_argument("--seed-start", type=int, default=None)
    parser.add_argument("--n-seeds", type=int, default=None)
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1f_f.yaml").read_text())
    if args.condition not in cfg["ranging_conditions"]:
        raise ValueError(f"Unknown condition: {args.condition}")
    if args.initial_scenario not in cfg["initial_mode_scenarios"]:
        raise ValueError(f"Unknown initial scenario: {args.initial_scenario}")

    validation = cfg["validation"]
    seed_start = validation["seed_start"] if args.seed_start is None else args.seed_start
    n_seeds = validation["n_seeds"] if args.n_seeds is None else args.n_seeds
    runs = [
        evaluate_one(cfg, args.condition, args.initial_scenario, seed)
        for seed in range(seed_start, seed_start + n_seeds)
    ]
    aggregate = aggregate_runs(runs)
    payload = {
        "design": {
            "condition": args.condition,
            "initial_scenario": args.initial_scenario,
            "seed_start": seed_start,
            "n_seeds": n_seeds,
            "n_particles": cfg["particle_cloud"]["n_particles"],
            "frozen_aacopf": {
                "alpha": cfg["frozen_aacopf"]["alpha"],
                "beta": cfg["frozen_aacopf"]["beta"],
                "c_lambda": cfg["frozen_aacopf"]["c_lambda"],
            },
            "claim_boundary": cfg["claim_boundary"]["description"],
        },
        "aggregate": aggregate,
        "runs": runs,
    }
    out_dir = ROOT / "results/phase1/p1f_f/shards"
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / f"{args.condition}__{args.initial_scenario}.json"
    output.write_text(json.dumps(payload, indent=2))
    print(json.dumps({"output": str(output), "aggregate": aggregate}, indent=2))


if __name__ == "__main__":
    main()
