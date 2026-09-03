from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml

from mt_ag.aacopf import literal_all_pairs_aco_transition, run_literal_small_aacopf
from mt_ag.geometry import wrap_angle
from mt_ag.paper_pf import generate_paper_trajectory
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory

ROOT = Path(__file__).resolve().parents[1]


def distribution(values):
    x = np.asarray(values, dtype=float)
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


def controlled_bimodal_particles(trajectory, auxiliary, cfg, rng):
    n_particles = cfg["n_particles"]
    n_correct = int(round(cfg["correct_fraction"] * n_particles))
    n_wrong = n_particles - n_correct
    position_std = cfg["position_std_m"]
    yaw_std = np.deg2rad(cfg["yaw_std_deg"])

    truth_position = trajectory.state[0, :2]
    truth_yaw = trajectory.state[0, 2]
    vector = truth_position - auxiliary[0]
    range0 = np.linalg.norm(vector)
    truth_bearing = np.arctan2(vector[1], vector[0])
    wrong_bearing = truth_bearing + np.deg2rad(cfg["wrong_bearing_offset_deg"])
    wrong_center = auxiliary[0] + range0 * np.array(
        [np.cos(wrong_bearing), np.sin(wrong_bearing)]
    )
    wrong_yaw = wrap_angle(truth_yaw + np.deg2rad(cfg["wrong_yaw_offset_deg"]))

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
    return particles, wrong_center, wrong_yaw


def run_controlled_campaign(cfg, quick):
    sim_cfg = cfg["simulation"].copy()
    if quick:
        sim_cfg["duration"] = min(sim_cfg["duration"], 3.0)
    trajectory = generate_paper_trajectory(**sim_cfg)
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    uwb_cfg = cfg["uwb"]
    aco_cfg = cfg["aco"]
    case_cfg = cfg["controlled_bimodal"].copy()
    if quick:
        case_cfg["n_seeds"] = min(case_cfg["n_seeds"], 2)
        case_cfg["n_particles"] = min(case_cfg["n_particles"], 120)
    randomness = cfg["randomness"]
    yaw_threshold = np.deg2rad(case_cfg["correct_mode_yaw_threshold_deg"])
    hold_steps = max(1, round(1.0 / trajectory.dt))

    runs = []
    for seed in range(case_cfg["n_seeds"]):
        range_rng = np.random.default_rng(randomness["range_seed_offset"] + seed)
        ranges = generate_uwb_ranges(
            trajectory.state[:, :2],
            auxiliary,
            range_rng,
            uwb_cfg["sigma_range_m"],
        )
        init_rng = np.random.default_rng(randomness["initializer_seed_offset"] + seed)
        particles, wrong_center, wrong_yaw = controlled_bimodal_particles(
            trajectory,
            auxiliary,
            case_cfg,
            init_rng,
        )
        result = run_literal_small_aacopf(
            trajectory.increments,
            ranges,
            auxiliary,
            particles,
            sigma_uwb_m=uwb_cfg["sigma_range_m"],
            propagation_convention="pre_turn",
            alpha=aco_cfg["alpha"],
            beta=aco_cfg["beta"],
            c_lambda=aco_cfg["c_lambda_mechanism_reference"],
            epsilon_distance=aco_cfg["epsilon_distance"],
            epsilon_weight=aco_cfg["epsilon_weight"],
            initial_particles_conditioned_on_z0=True,
            truth_state=trajectory.state,
            correct_mode_position_threshold_m=case_cfg[
                "correct_mode_position_threshold_m"
            ],
            correct_mode_yaw_threshold_rad=yaw_threshold,
        )

        position_error = np.linalg.norm(
            result.estimate[:, :2] - trajectory.state[:, :2], axis=1
        )
        yaw_error = np.abs(wrap_angle(result.estimate[:, 2] - trajectory.state[:, 2]))
        mode_fraction = result.correct_mode_fraction_post_aco
        mode_concentration_time = terminal_hold_time(
            mode_fraction >= 0.9,
            trajectory.dt,
            hold_steps,
        )
        update_slice = slice(1, None)
        map_lineage_error = abs(
            wrap_angle(result.map_initial_yaw[-1] - trajectory.state[0, 2])
        )
        runs.append(
            {
                "seed": seed,
                "n_particles": case_cfg["n_particles"],
                "wrong_mode_center_m": wrong_center.tolist(),
                "wrong_mode_yaw_deg": float(np.rad2deg(wrong_yaw)),
                "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
                "yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error**2)))),
                "final_position_error_m": float(position_error[-1]),
                "final_yaw_error_deg": float(np.rad2deg(yaw_error[-1])),
                "final_map_initial_yaw_error_deg": float(np.rad2deg(map_lineage_error)),
                "initial_correct_mode_fraction_post_aco": float(mode_fraction[0]),
                "final_correct_mode_fraction_post_aco": float(mode_fraction[-1]),
                "final_correct_mode_mass_pre_aco": float(
                    result.correct_mode_mass_pre_aco[-1]
                ),
                "mode_concentration_time_s": mode_concentration_time,
                "mean_moved_fraction": float(np.mean(result.moved_fraction[update_slice])),
                "mean_unique_parent_fraction": float(
                    np.mean(result.unique_parent_fraction[update_slice])
                ),
                "min_unique_parent_fraction": float(
                    np.min(result.unique_parent_fraction[update_slice])
                ),
                "mean_unique_destination_count": float(
                    np.mean(result.unique_destination_count[update_slice])
                ),
                "mean_destination_multiplicity": float(
                    np.mean(result.mean_destination_multiplicity[update_slice])
                ),
                "max_destination_multiplicity": int(
                    np.max(result.max_destination_multiplicity[update_slice])
                ),
                "mean_candidate_count": float(
                    np.mean(result.mean_candidate_count[update_slice])
                ),
                "mean_max_transition_probability": float(
                    np.mean(result.mean_max_probability[update_slice])
                ),
                "max_transition_probability": float(
                    np.max(result.max_probability[update_slice])
                ),
                "mean_neff_fraction_pre_aco": float(
                    np.mean(result.neff_pre_aco[update_slice]) / case_cfg["n_particles"]
                ),
                "mean_transition_runtime_ms": float(
                    1000.0 * np.mean(result.transition_runtime_s[update_slice])
                ),
                "total_transition_runtime_s": float(
                    np.sum(result.transition_runtime_s[update_slice])
                ),
            }
        )

    successful_times = [
        run["mode_concentration_time_s"]
        for run in runs
        if run["mode_concentration_time_s"] is not None
    ]
    aggregate = {
        "n": len(runs),
        "n_particles": case_cfg["n_particles"],
        "mode_concentration_fraction": float(
            np.mean([run["mode_concentration_time_s"] is not None for run in runs])
        ),
        "mode_concentration_time_s": (
            distribution(successful_times) if successful_times else None
        ),
        "position_rmse_m": distribution([run["position_rmse_m"] for run in runs]),
        "yaw_rmse_deg": distribution([run["yaw_rmse_deg"] for run in runs]),
        "final_correct_mode_fraction_post_aco": distribution(
            [run["final_correct_mode_fraction_post_aco"] for run in runs]
        ),
        "final_correct_mode_mass_pre_aco": distribution(
            [run["final_correct_mode_mass_pre_aco"] for run in runs]
        ),
        "mean_moved_fraction": distribution(
            [run["mean_moved_fraction"] for run in runs]
        ),
        "mean_unique_parent_fraction": distribution(
            [run["mean_unique_parent_fraction"] for run in runs]
        ),
        "max_destination_multiplicity": distribution(
            [run["max_destination_multiplicity"] for run in runs]
        ),
        "mean_max_transition_probability": distribution(
            [run["mean_max_transition_probability"] for run in runs]
        ),
        "mean_neff_fraction_pre_aco": distribution(
            [run["mean_neff_fraction_pre_aco"] for run in runs]
        ),
        "mean_transition_runtime_ms": distribution(
            [run["mean_transition_runtime_ms"] for run in runs]
        ),
    }
    return {"aggregate": aggregate, "runs": runs}


def run_runtime_scaling(cfg, quick):
    aco_cfg = cfg["aco"]
    runtime_cfg = cfg["runtime_scaling"].copy()
    counts = runtime_cfg["particle_counts"]
    repeats = runtime_cfg["repeats"]
    warmups = runtime_cfg["warmup_repeats"]
    if quick:
        counts = [50, 100, 200]
        repeats = 2
        warmups = 1
    rng = np.random.default_rng(cfg["randomness"]["runtime_seed"])
    rows = []

    for n_particles in counts:
        particles = np.column_stack(
            [
                rng.normal(0.0, 5.0, n_particles),
                rng.normal(0.0, 5.0, n_particles),
                rng.uniform(-np.pi, np.pi, n_particles),
            ]
        )
        weights = rng.random(n_particles) + 1.0e-6
        weights /= np.sum(weights)
        lineage = particles[:, 2].copy()

        for _ in range(warmups):
            literal_all_pairs_aco_transition(
                particles,
                weights,
                lineage,
                alpha=aco_cfg["alpha"],
                beta=aco_cfg["beta"],
                c_lambda=aco_cfg["c_lambda_mechanism_reference"],
                epsilon_distance=aco_cfg["epsilon_distance"],
                epsilon_weight=aco_cfg["epsilon_weight"],
            )

        runtimes = []
        last_diag = None
        for _ in range(repeats):
            _, _, _, last_diag = literal_all_pairs_aco_transition(
                particles,
                weights,
                lineage,
                alpha=aco_cfg["alpha"],
                beta=aco_cfg["beta"],
                c_lambda=aco_cfg["c_lambda_mechanism_reference"],
                epsilon_distance=aco_cfg["epsilon_distance"],
                epsilon_weight=aco_cfg["epsilon_weight"],
            )
            runtimes.append(last_diag.runtime_s)

        rows.append(
            {
                "n_particles": n_particles,
                "dense_pair_count": last_diag.dense_pair_count,
                "candidate_score_count": last_diag.candidate_score_count,
                "runtime_s": distribution(runtimes),
            }
        )

    log_n = np.log([row["n_particles"] for row in rows])
    log_runtime = np.log([row["runtime_s"]["median"] for row in rows])
    slope, intercept = np.polyfit(log_n, log_runtime, 1)
    return {
        "rows": rows,
        "empirical_loglog_runtime_slope": float(slope),
        "empirical_loglog_runtime_intercept": float(intercept),
        "theoretical_dense_pair_scaling": "N^2",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1f_b.yaml").read_text())
    controlled = run_controlled_campaign(cfg, args.quick)
    runtime_scaling = run_runtime_scaling(cfg, args.quick)

    raw = {
        "design": {
            "purpose": "literal-small ACO mechanism validation, not tuning",
            "state": "[x,y,phi]",
            "inputs": "[delta_L,delta_phi]",
            "propagation": "pre_turn Equation (3)",
            "aco_update": "synchronous all-pairs strictly-higher-weight destinations",
            "post_transition_weights": "uniform",
            "alpha": cfg["aco"]["alpha"],
            "beta": cfg["aco"]["beta"],
            "c_lambda": cfg["aco"]["c_lambda_mechanism_reference"],
            "c_lambda_status": "provisional mechanism reference; not tuned or frozen",
            "quick": args.quick,
        },
        "controlled_bimodal": controlled,
        "runtime_scaling": runtime_scaling,
    }

    aggregate = controlled["aggregate"]
    summary = {
        "controlled_bimodal": {
            "n": aggregate["n"],
            "n_particles": aggregate["n_particles"],
            "mode_concentration_fraction": aggregate["mode_concentration_fraction"],
            "mode_concentration_time_median_s": (
                None
                if aggregate["mode_concentration_time_s"] is None
                else aggregate["mode_concentration_time_s"]["median"]
            ),
            "position_rmse_mean_m": aggregate["position_rmse_m"]["mean"],
            "yaw_rmse_mean_deg": aggregate["yaw_rmse_deg"]["mean"],
            "final_correct_mode_fraction_mean": aggregate[
                "final_correct_mode_fraction_post_aco"
            ]["mean"],
            "final_correct_mode_mass_pre_aco_mean": aggregate[
                "final_correct_mode_mass_pre_aco"
            ]["mean"],
            "mean_moved_fraction": aggregate["mean_moved_fraction"]["mean"],
            "mean_unique_parent_fraction": aggregate[
                "mean_unique_parent_fraction"
            ]["mean"],
            "max_destination_multiplicity_mean": aggregate[
                "max_destination_multiplicity"
            ]["mean"],
            "mean_max_transition_probability": aggregate[
                "mean_max_transition_probability"
            ]["mean"],
            "mean_transition_runtime_ms": aggregate[
                "mean_transition_runtime_ms"
            ]["mean"],
        },
        "runtime_scaling": runtime_scaling,
    }

    out_dir = ROOT / "results/phase1/p1f_b"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
