from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.aacopf import run_literal_small_aacopf
from mt_ag.geometry import wrap_angle
from mt_ag.paper_pf import (
    generate_paper_trajectory,
    initialize_random_annulus_particles,
    run_paper_bootstrap_pf,
)
from mt_ag.sensors import generate_uwb_ranges
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
        "q95": float(np.quantile(x, 0.95)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def terminal_convergence_time(mask, dt, hold_steps):
    mask = np.asarray(mask, dtype=bool)
    false_indices = np.flatnonzero(~mask)
    start = 0 if len(false_indices) == 0 else int(false_indices[-1] + 1)
    if len(mask) - start < hold_steps:
        return None
    return float(start * dt)


def initial_support_count(particles, truth0, position_threshold_m, yaw_threshold_rad):
    position_error = np.linalg.norm(particles[:, :2] - truth0[:2], axis=1)
    yaw_error = np.abs(wrap_angle(particles[:, 2] - truth0[2]))
    return int(np.sum((position_error < position_threshold_m) & (yaw_error < yaw_threshold_rad)))


def terminal_metrics(estimate, truth, conv_cfg, dt):
    position_error = np.linalg.norm(estimate[:, :2] - truth[:, :2], axis=1)
    yaw_error = np.abs(wrap_angle(estimate[:, 2] - truth[:, 2]))
    yaw_threshold = np.deg2rad(conv_cfg["yaw_threshold_deg"])
    hold_steps = round(conv_cfg["hold_time_s"] / dt)
    late_steps = round(conv_cfg["late_window_s"] / dt)
    late = slice(max(0, len(position_error) - late_steps), None)
    pose_ok = (position_error < conv_cfg["position_threshold_m"]) & (yaw_error < yaw_threshold)
    return {
        "pose_convergence_s": terminal_convergence_time(pose_ok, dt, hold_steps),
        "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
        "late_position_rmse_m": float(np.sqrt(np.mean(position_error[late] ** 2))),
        "final_position_error_m": float(position_error[-1]),
        "yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error**2)))),
        "late_yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error[late] ** 2)))),
        "final_yaw_error_deg": float(np.rad2deg(yaw_error[-1])),
    }


def aggregate_method(runs, method):
    values = [run[method] for run in runs]
    convergence_times = [
        value["pose_convergence_s"]
        for value in values
        if value["pose_convergence_s"] is not None
    ]
    aggregate = {
        "n": len(values),
        "pose_success_fraction": float(
            np.mean([value["pose_convergence_s"] is not None for value in values])
        ),
        "pose_convergence_s": distribution(convergence_times),
        "late_position_rmse_m": distribution(
            [value["late_position_rmse_m"] for value in values]
        ),
        "late_yaw_rmse_deg": distribution([value["late_yaw_rmse_deg"] for value in values]),
        "final_position_error_m": distribution(
            [value["final_position_error_m"] for value in values]
        ),
        "runtime_s": distribution([value["runtime_s"] for value in values]),
    }
    if method == "pf":
        aggregate.update(
            {
                "final_correct_mode_mass": distribution(
                    [value["final_correct_mode_mass"] for value in values]
                ),
                "resampling_events": distribution(
                    [value["resampling_events"] for value in values]
                ),
                "minimum_unique_fraction": distribution(
                    [value["minimum_unique_fraction"] for value in values]
                ),
            }
        )
    else:
        aggregate.update(
            {
                "final_correct_mode_mass_pre_aco": distribution(
                    [value["final_correct_mode_mass_pre_aco"] for value in values]
                ),
                "final_correct_mode_fraction_post_aco": distribution(
                    [value["final_correct_mode_fraction_post_aco"] for value in values]
                ),
                "minimum_unique_parent_fraction": distribution(
                    [value["minimum_unique_parent_fraction"] for value in values]
                ),
                "catastrophic_collapse_fraction": float(
                    np.mean([value["catastrophic_collapse"] for value in values])
                ),
                "dominant_clone_fraction": float(
                    np.mean([value["dominant_clone"] for value in values])
                ),
                "mean_moved_fraction": distribution(
                    [value["mean_moved_fraction"] for value in values]
                ),
                "maximum_destination_multiplicity": distribution(
                    [value["maximum_destination_multiplicity"] for value in values]
                ),
                "transition_runtime_s": distribution(
                    [value["transition_runtime_s"] for value in values]
                ),
            }
        )
    return aggregate


def support_group(runs, predicate):
    selected = [run for run in runs if predicate(run["initial_support_count"])]
    if not selected:
        return None
    return {
        "n": len(selected),
        "support_count": distribution([run["initial_support_count"] for run in selected]),
        "pf_pose_success_fraction": float(
            np.mean([run["pf"]["pose_convergence_s"] is not None for run in selected])
        ),
        "aacopf_pose_success_fraction": float(
            np.mean([run["aacopf"]["pose_convergence_s"] is not None for run in selected])
        ),
    }


def paired_summary(runs):
    pf_success = np.array([run["pf"]["pose_convergence_s"] is not None for run in runs])
    aco_success = np.array([run["aacopf"]["pose_convergence_s"] is not None for run in runs])
    pos_delta = np.array(
        [
            run["aacopf"]["late_position_rmse_m"] - run["pf"]["late_position_rmse_m"]
            for run in runs
        ]
    )
    yaw_delta = np.array(
        [run["aacopf"]["late_yaw_rmse_deg"] - run["pf"]["late_yaw_rmse_deg"] for run in runs]
    )
    return {
        "both_success": int(np.sum(pf_success & aco_success)),
        "aacopf_only_success": int(np.sum(~pf_success & aco_success)),
        "pf_only_success": int(np.sum(pf_success & ~aco_success)),
        "neither_success": int(np.sum(~pf_success & ~aco_success)),
        "aacopf_minus_pf_late_position_rmse_m": distribution(pos_delta),
        "aacopf_minus_pf_late_yaw_rmse_deg": distribution(yaw_delta),
        "aacopf_lower_late_position_rmse_fraction": float(np.mean(pos_delta < 0.0)),
        "aacopf_lower_late_yaw_rmse_fraction": float(np.mean(yaw_delta < 0.0)),
    }


def compact_budget(value):
    pf = value["pf"]
    aco = value["aacopf"]
    paired = value["paired"]
    return {
        "n": value["n"],
        "n_particles": value["n_particles"],
        "initial_support_positive_fraction": value["initial_support_positive_fraction"],
        "initial_support_count_mean": value["initial_support_count"]["mean"],
        "pf_pose_success_fraction": pf["pose_success_fraction"],
        "aacopf_pose_success_fraction": aco["pose_success_fraction"],
        "pf_late_position_rmse_mean_m": pf["late_position_rmse_m"]["mean"],
        "aacopf_late_position_rmse_mean_m": aco["late_position_rmse_m"]["mean"],
        "pf_late_yaw_rmse_mean_deg": pf["late_yaw_rmse_deg"]["mean"],
        "aacopf_late_yaw_rmse_mean_deg": aco["late_yaw_rmse_deg"]["mean"],
        "aacopf_only_success": paired["aacopf_only_success"],
        "pf_only_success": paired["pf_only_success"],
        "both_success": paired["both_success"],
        "neither_success": paired["neither_success"],
        "aacopf_catastrophic_collapse_fraction": aco["catastrophic_collapse_fraction"],
        "aacopf_dominant_clone_fraction": aco["dominant_clone_fraction"],
        "pf_runtime_mean_s": pf["runtime_s"]["mean"],
        "aacopf_runtime_mean_s": aco["runtime_s"]["mean"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1f_d.yaml").read_text())
    sim_cfg = dict(cfg["simulation"])
    comparison_cfg = cfg["comparison"]
    if args.quick:
        sim_cfg["duration"] = min(sim_cfg["duration"], 5.0)
    trajectory = generate_paper_trajectory(**sim_cfg)
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    sigma_uwb = cfg["uwb"]["sigma_range_m"]
    delta_d = cfg["initialization"]["delta_d_sigma_factor"] * sigma_uwb
    conv_cfg = cfg["convergence"]
    support_cfg = cfg["support_diagnostics"]
    ancestry_cfg = cfg["ancestry_diagnostics"]
    aco_cfg = cfg["frozen_aacopf"]
    pf_cfg = cfg["conventional_pf"]
    model_cfg = cfg["paper_model"]
    seed_cfg = cfg["randomness"]

    n_seeds = min(2, comparison_cfg["n_seeds"]) if args.quick else comparison_cfg["n_seeds"]
    budgets = [100, 200] if args.quick else comparison_cfg["particle_budgets"]
    support_yaw_threshold = np.deg2rad(support_cfg["correct_region_yaw_threshold_deg"])
    aco_yaw_threshold = np.deg2rad(conv_cfg["yaw_threshold_deg"])

    results_by_budget = {}
    for budget_index, n_particles in enumerate(budgets):
        runs = []
        for seed in range(n_seeds):
            range_rng = np.random.default_rng(seed_cfg["range_seed_offset"] + seed)
            ranges = generate_uwb_ranges(
                trajectory.state[:, :2], auxiliary, range_rng, sigma_uwb
            )
            init_seed = seed_cfg["initializer_seed_offset"] + 1000 * budget_index + seed
            init_rng = np.random.default_rng(init_seed)
            initial_particles = initialize_random_annulus_particles(
                ranges[0], auxiliary[0], n_particles, init_rng, delta_d
            )
            support_count = initial_support_count(
                initial_particles,
                trajectory.state[0],
                support_cfg["correct_region_position_threshold_m"],
                support_yaw_threshold,
            )

            pf_rng = np.random.default_rng(
                seed_cfg["pf_seed_offset"] + 1000 * budget_index + seed
            )
            start = time.perf_counter()
            pf_result = run_paper_bootstrap_pf(
                trajectory.increments,
                ranges,
                auxiliary,
                initial_particles.copy(),
                pf_rng,
                sigma_uwb_m=sigma_uwb,
                resample_fraction=pf_cfg["resample_fraction"],
                propagation_convention=model_cfg["propagation_convention"],
                sigma_delta_l_m=model_cfg["sigma_delta_l_m"],
                sigma_delta_phi_rad=model_cfg["sigma_delta_phi_rad"],
                initial_particles_conditioned_on_z0=True,
                truth_state=trajectory.state,
                correct_mode_position_threshold_m=conv_cfg["position_threshold_m"],
                correct_mode_yaw_threshold_rad=aco_yaw_threshold,
            )
            pf_runtime = time.perf_counter() - start
            pf_metrics = terminal_metrics(pf_result.estimate, trajectory.state, conv_cfg, trajectory.dt)
            pf_resampled = pf_result.unique_fraction_post_transition < 1.0 - 1e-15
            pf_metrics.update(
                {
                    "runtime_s": float(pf_runtime),
                    "final_correct_mode_mass": float(pf_result.correct_mode_mass[-1]),
                    "resampling_events": int(np.sum(pf_resampled)),
                    "minimum_unique_fraction": float(
                        np.min(pf_result.unique_fraction_post_transition)
                    ),
                }
            )

            start = time.perf_counter()
            aco_result = run_literal_small_aacopf(
                trajectory.increments,
                ranges,
                auxiliary,
                initial_particles.copy(),
                sigma_uwb_m=sigma_uwb,
                propagation_convention=model_cfg["propagation_convention"],
                alpha=aco_cfg["alpha"],
                beta=aco_cfg["beta"],
                c_lambda=aco_cfg["c_lambda"],
                epsilon_distance=aco_cfg["epsilon_distance"],
                epsilon_weight=aco_cfg["epsilon_weight"],
                initial_particles_conditioned_on_z0=True,
                truth_state=trajectory.state,
                correct_mode_position_threshold_m=conv_cfg["position_threshold_m"],
                correct_mode_yaw_threshold_rad=aco_yaw_threshold,
            )
            aco_runtime = time.perf_counter() - start
            aco_metrics = terminal_metrics(
                aco_result.estimate, trajectory.state, conv_cfg, trajectory.dt
            )
            update_slice = slice(1, None)
            min_parent = float(np.min(aco_result.unique_parent_fraction[update_slice]))
            max_multiplicity = int(
                np.max(aco_result.max_destination_multiplicity[update_slice])
            )
            aco_metrics.update(
                {
                    "runtime_s": float(aco_runtime),
                    "transition_runtime_s": float(
                        np.sum(aco_result.transition_runtime_s[update_slice])
                    ),
                    "final_correct_mode_mass_pre_aco": float(
                        aco_result.correct_mode_mass_pre_aco[-1]
                    ),
                    "final_correct_mode_fraction_post_aco": float(
                        aco_result.correct_mode_fraction_post_aco[-1]
                    ),
                    "minimum_unique_parent_fraction": min_parent,
                    "catastrophic_collapse": bool(
                        min_parent < ancestry_cfg["catastrophic_unique_parent_fraction"]
                    ),
                    "dominant_clone": bool(
                        max_multiplicity >= ancestry_cfg["dominant_clone_fraction"] * n_particles
                    ),
                    "maximum_destination_multiplicity": max_multiplicity,
                    "mean_moved_fraction": float(
                        np.mean(aco_result.moved_fraction[update_slice])
                    ),
                }
            )

            runs.append(
                {
                    "seed": seed,
                    "n_particles": n_particles,
                    "initial_support_count": support_count,
                    "initial_support_fraction": support_count / n_particles,
                    "pf": pf_metrics,
                    "aacopf": aco_metrics,
                }
            )

        budget_result = {
            "n": len(runs),
            "n_particles": n_particles,
            "initial_support_count": distribution(
                [run["initial_support_count"] for run in runs]
            ),
            "initial_support_positive_fraction": float(
                np.mean([run["initial_support_count"] > 0 for run in runs])
            ),
            "pf": aggregate_method(runs, "pf"),
            "aacopf": aggregate_method(runs, "aacopf"),
            "paired": paired_summary(runs),
            "support_groups": {
                "zero": support_group(runs, lambda count: count == 0),
                "one_to_two": support_group(runs, lambda count: 1 <= count <= 2),
                "three_or_more": support_group(runs, lambda count: count >= 3),
                "positive": support_group(runs, lambda count: count > 0),
            },
            "runs": runs,
        }
        results_by_budget[str(n_particles)] = budget_result

    primary_budget = (
        budgets[-1] if args.quick else int(comparison_cfg["primary_particle_budget"])
    )
    raw = {
        "design": {
            "state": "[x,y,phi]",
            "inputs": "[delta_L,delta_phi]",
            "geometry": "informative moving globally-known auxiliary",
            "initialization": "matched random first-range annulus + random yaw",
            "dr_increments": "exact",
            "uwb_sigma_m": sigma_uwb,
            "delta_d_m": delta_d,
            "frozen_aacopf": {
                "alpha": aco_cfg["alpha"],
                "beta": aco_cfg["beta"],
                "c_lambda": aco_cfg["c_lambda"],
            },
            "particle_budgets": budgets,
            "n_seeds": n_seeds,
            "primary_particle_budget": primary_budget,
            "quick": args.quick,
        },
        "budgets": results_by_budget,
    }
    summary = {
        "design": raw["design"],
        "budgets": {
            key: compact_budget(value) for key, value in results_by_budget.items()
        },
        "primary": compact_budget(results_by_budget[str(primary_budget)]),
    }

    out_dir = ROOT / "results/phase1/p1f_d"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
