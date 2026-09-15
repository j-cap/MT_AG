from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.aacopf import run_literal_small_aacopf
from mt_ag.geometry import wrap_angle
from mt_ag.observability import (
    observability_history_paper_state,
    relative_geometry_metrics,
)
from mt_ag.paper_pf import (
    generate_paper_trajectory,
    initialize_random_annulus_particles,
    run_paper_bootstrap_pf,
)
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


def method_metrics(estimate, map_state, truth, auxiliary, conv_cfg, dt):
    position_error = np.linalg.norm(estimate[:, :2] - truth[:, :2], axis=1)
    yaw_error = np.abs(wrap_angle(estimate[:, 2] - truth[:, 2]))
    yaw_threshold = np.deg2rad(conv_cfg["yaw_threshold_deg"])
    hold_steps = max(1, round(conv_cfg["hold_time_s"] / dt))
    late_steps = max(1, round(conv_cfg["late_window_s"] / dt))
    late = slice(max(0, len(position_error) - late_steps), None)
    pose_ok = (position_error < conv_cfg["position_threshold_m"]) & (yaw_error < yaw_threshold)

    true_ranges = np.linalg.norm(truth[:, :2] - auxiliary, axis=1)
    map_ranges = np.linalg.norm(map_state[:, :2] - auxiliary, axis=1)
    map_range_error = map_ranges - true_ranges
    return {
        "pose_convergence_s": terminal_convergence_time(pose_ok, dt, hold_steps),
        "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
        "late_position_rmse_m": float(np.sqrt(np.mean(position_error[late] ** 2))),
        "final_position_error_m": float(position_error[-1]),
        "yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error**2)))),
        "late_yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error[late] ** 2)))),
        "final_yaw_error_deg": float(np.rad2deg(yaw_error[-1])),
        "map_range_fit_rmse_m": float(np.sqrt(np.mean(map_range_error**2))),
        "late_map_range_fit_rmse_m": float(np.sqrt(np.mean(map_range_error[late] ** 2))),
    }


def aggregate_method(runs, method):
    values = [run[method] for run in runs]
    convergence = [v["pose_convergence_s"] for v in values if v["pose_convergence_s"] is not None]
    out = {
        "n": len(values),
        "pose_success_fraction": float(np.mean([v["pose_convergence_s"] is not None for v in values])),
        "pose_convergence_s": distribution(convergence),
        "late_position_rmse_m": distribution([v["late_position_rmse_m"] for v in values]),
        "late_yaw_rmse_deg": distribution([v["late_yaw_rmse_deg"] for v in values]),
        "final_position_error_m": distribution([v["final_position_error_m"] for v in values]),
        "late_map_range_fit_rmse_m": distribution([v["late_map_range_fit_rmse_m"] for v in values]),
        "final_position_spread_m": distribution([v["final_position_spread_m"] for v in values]),
        "final_yaw_resultant": distribution([v["final_yaw_resultant"] for v in values]),
        "runtime_s": distribution([v["runtime_s"] for v in values]),
    }
    if method == "pf":
        out.update(
            {
                "final_correct_mode_mass": distribution([v["final_correct_mode_mass"] for v in values]),
                "resampling_events": distribution([v["resampling_events"] for v in values]),
                "minimum_unique_fraction": distribution([v["minimum_unique_fraction"] for v in values]),
            }
        )
    else:
        out.update(
            {
                "final_correct_mode_mass_pre_aco": distribution(
                    [v["final_correct_mode_mass_pre_aco"] for v in values]
                ),
                "final_correct_mode_fraction_post_aco": distribution(
                    [v["final_correct_mode_fraction_post_aco"] for v in values]
                ),
                "minimum_unique_parent_fraction": distribution(
                    [v["minimum_unique_parent_fraction"] for v in values]
                ),
                "catastrophic_collapse_fraction": float(
                    np.mean([v["catastrophic_collapse"] for v in values])
                ),
                "dominant_clone_fraction": float(np.mean([v["dominant_clone"] for v in values])),
                "mean_moved_fraction": distribution([v["mean_moved_fraction"] for v in values]),
                "maximum_destination_multiplicity": distribution(
                    [v["maximum_destination_multiplicity"] for v in values]
                ),
                "transition_runtime_s": distribution([v["transition_runtime_s"] for v in values]),
            }
        )
    return out


def compact_geometry(result):
    return {
        "n": result["n"],
        "n_particles": result["n_particles"],
        "observability_rank": result["observability"]["final_rank"],
        "observability_sigma_ratio": result["observability"]["final_sigma_ratio"],
        "initial_support_count_mean": result["initial_support_count"]["mean"],
        "pf_pose_success_fraction": result["pf"]["pose_success_fraction"],
        "aacopf_pose_success_fraction": result["aacopf"]["pose_success_fraction"],
        "pf_late_position_rmse_mean_m": result["pf"]["late_position_rmse_m"]["mean"],
        "aacopf_late_position_rmse_mean_m": result["aacopf"]["late_position_rmse_m"]["mean"],
        "pf_late_yaw_rmse_mean_deg": result["pf"]["late_yaw_rmse_deg"]["mean"],
        "aacopf_late_yaw_rmse_mean_deg": result["aacopf"]["late_yaw_rmse_deg"]["mean"],
        "pf_late_map_range_fit_rmse_mean_m": result["pf"]["late_map_range_fit_rmse_m"]["mean"],
        "aacopf_late_map_range_fit_rmse_mean_m": result["aacopf"]["late_map_range_fit_rmse_m"]["mean"],
        "pf_final_position_spread_mean_m": result["pf"]["final_position_spread_m"]["mean"],
        "aacopf_final_position_spread_mean_m": result["aacopf"]["final_position_spread_m"]["mean"],
        "pf_final_yaw_resultant_mean": result["pf"]["final_yaw_resultant"]["mean"],
        "aacopf_final_yaw_resultant_mean": result["aacopf"]["final_yaw_resultant"]["mean"],
        "aacopf_catastrophic_collapse_fraction": result["aacopf"]["catastrophic_collapse_fraction"],
        "aacopf_dominant_clone_fraction": result["aacopf"]["dominant_clone_fraction"],
    }


def observability_summary(trajectory, auxiliary, cfg):
    history = observability_history_paper_state(
        trajectory.state,
        trajectory.increments,
        auxiliary,
        convention=cfg["paper_model"]["propagation_convention"],
        rank_rtol=cfg["observability"]["rank_rtol"],
    )
    return {
        "final_singular_values": [float(v) for v in history["singular_values"][-1]],
        "final_rank": int(history["rank"][-1]),
        "final_sigma_ratio": float(history["sigma_ratio"][-1]),
        "geometry": relative_geometry_metrics(trajectory.state[:, :2], auxiliary),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1f_e.yaml").read_text())
    sim_cfg = dict(cfg["simulation"])
    if args.quick:
        sim_cfg["duration"] = min(5.0, sim_cfg["duration"])
    trajectory = generate_paper_trajectory(**sim_cfg)

    geometries = {}
    for name in cfg["geometry_cases"]:
        geometries[name] = auxiliary_trajectory(
            trajectory.t,
            name,
            target_positions=trajectory.state[:, :2] if name == "constant_bearing" else None,
        )
    moving_reference = auxiliary_trajectory(trajectory.t, "moving")

    observability = {
        name: observability_summary(trajectory, auxiliary, cfg)
        for name, auxiliary in geometries.items()
    }
    observability["moving_reference"] = observability_summary(trajectory, moving_reference, cfg)

    n_seeds = min(2, cfg["comparison"]["n_seeds"]) if args.quick else cfg["comparison"]["n_seeds"]
    n_particles = 200 if args.quick else cfg["comparison"]["n_particles"]
    sigma_uwb = cfg["uwb"]["sigma_range_m"]
    delta_d = cfg["initialization"]["delta_d_sigma_factor"] * sigma_uwb
    support_yaw = np.deg2rad(cfg["support_diagnostics"]["correct_region_yaw_threshold_deg"])
    conv_yaw = np.deg2rad(cfg["convergence"]["yaw_threshold_deg"])
    seed_cfg = cfg["randomness"]
    aco_cfg = cfg["frozen_aacopf"]

    runs_by_geometry = {name: [] for name in geometries}
    for seed in range(n_seeds):
        noise_rng = np.random.default_rng(seed_cfg["range_seed_offset"] + seed)
        shared_noise = noise_rng.normal(0.0, sigma_uwb, len(trajectory.t))
        ranges_by_geometry = {
            name: np.linalg.norm(trajectory.state[:, :2] - auxiliary, axis=1) + shared_noise
            for name, auxiliary in geometries.items()
        }
        z0_values = [values[0] for values in ranges_by_geometry.values()]
        if not np.allclose(z0_values, z0_values[0], atol=1e-12, rtol=0.0):
            raise RuntimeError("negative-control geometries must share the same initial range")

        init_rng = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        initial_particles = initialize_random_annulus_particles(
            z0_values[0],
            next(iter(geometries.values()))[0],
            n_particles,
            init_rng,
            delta_d,
        )
        support_count = initial_support_count(
            initial_particles,
            trajectory.state[0],
            cfg["support_diagnostics"]["correct_region_position_threshold_m"],
            support_yaw,
        )

        for geometry_index, (name, auxiliary) in enumerate(geometries.items()):
            ranges = ranges_by_geometry[name]
            pf_rng = np.random.default_rng(
                seed_cfg["pf_seed_offset"] + 1000 * geometry_index + seed
            )
            start = time.perf_counter()
            pf_result = run_paper_bootstrap_pf(
                trajectory.increments,
                ranges,
                auxiliary,
                initial_particles.copy(),
                pf_rng,
                sigma_uwb_m=sigma_uwb,
                resample_fraction=cfg["conventional_pf"]["resample_fraction"],
                propagation_convention=cfg["paper_model"]["propagation_convention"],
                sigma_delta_l_m=cfg["paper_model"]["sigma_delta_l_m"],
                sigma_delta_phi_rad=cfg["paper_model"]["sigma_delta_phi_rad"],
                initial_particles_conditioned_on_z0=True,
                truth_state=trajectory.state,
                correct_mode_position_threshold_m=cfg["convergence"]["position_threshold_m"],
                correct_mode_yaw_threshold_rad=conv_yaw,
            )
            pf_runtime = time.perf_counter() - start
            pf_metrics = method_metrics(
                pf_result.estimate,
                pf_result.map_state,
                trajectory.state,
                auxiliary,
                cfg["convergence"],
                trajectory.dt,
            )
            pf_metrics.update(
                {
                    "runtime_s": float(pf_runtime),
                    "final_correct_mode_mass": float(pf_result.correct_mode_mass[-1]),
                    "resampling_events": int(
                        np.sum(pf_result.unique_fraction_post_transition < 1.0 - 1e-15)
                    ),
                    "minimum_unique_fraction": float(
                        np.min(pf_result.unique_fraction_post_transition)
                    ),
                    "final_position_spread_m": float(pf_result.position_spread[-1]),
                    "final_yaw_resultant": float(pf_result.yaw_resultant[-1]),
                }
            )

            start = time.perf_counter()
            aco_result = run_literal_small_aacopf(
                trajectory.increments,
                ranges,
                auxiliary,
                initial_particles.copy(),
                sigma_uwb_m=sigma_uwb,
                propagation_convention=cfg["paper_model"]["propagation_convention"],
                alpha=aco_cfg["alpha"],
                beta=aco_cfg["beta"],
                c_lambda=aco_cfg["c_lambda"],
                epsilon_distance=aco_cfg["epsilon_distance"],
                epsilon_weight=aco_cfg["epsilon_weight"],
                initial_particles_conditioned_on_z0=True,
                truth_state=trajectory.state,
                correct_mode_position_threshold_m=cfg["convergence"]["position_threshold_m"],
                correct_mode_yaw_threshold_rad=conv_yaw,
            )
            aco_runtime = time.perf_counter() - start
            aco_metrics = method_metrics(
                aco_result.estimate,
                aco_result.map_state,
                trajectory.state,
                auxiliary,
                cfg["convergence"],
                trajectory.dt,
            )
            update_slice = slice(1, None)
            min_parent = float(np.min(aco_result.unique_parent_fraction[update_slice]))
            max_multiplicity = int(np.max(aco_result.max_destination_multiplicity[update_slice]))
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
                        min_parent < cfg["ancestry_diagnostics"]["catastrophic_unique_parent_fraction"]
                    ),
                    "dominant_clone": bool(
                        max_multiplicity
                        >= cfg["ancestry_diagnostics"]["dominant_clone_fraction"] * n_particles
                    ),
                    "maximum_destination_multiplicity": max_multiplicity,
                    "mean_moved_fraction": float(
                        np.mean(aco_result.moved_fraction[update_slice])
                    ),
                    "final_position_spread_m": float(aco_result.position_spread[-1]),
                    "final_yaw_resultant": float(aco_result.yaw_resultant[-1]),
                }
            )

            runs_by_geometry[name].append(
                {
                    "seed": seed,
                    "initial_support_count": support_count,
                    "pf": pf_metrics,
                    "aacopf": aco_metrics,
                }
            )

    raw_geometry = {}
    summary_geometry = {}
    for name, runs in runs_by_geometry.items():
        result = {
            "n": len(runs),
            "n_particles": n_particles,
            "observability": observability[name],
            "initial_support_count": distribution([r["initial_support_count"] for r in runs]),
            "initial_support_positive_fraction": float(
                np.mean([r["initial_support_count"] > 0 for r in runs])
            ),
            "pf": aggregate_method(runs, "pf"),
            "aacopf": aggregate_method(runs, "aacopf"),
            "runs": runs,
        }
        raw_geometry[name] = result
        summary_geometry[name] = compact_geometry(result)

    summary = {
        "frozen_aacopf": {
            "alpha": aco_cfg["alpha"],
            "beta": aco_cfg["beta"],
            "c_lambda": aco_cfg["c_lambda"],
        },
        "negative_controls": summary_geometry,
        "moving_observability_reference": observability["moving_reference"],
        "acceptance": {
            "criterion": "Neither PF nor frozen AACOPF should show systematic full-pose recovery in locally rank-deficient geometries.",
            "aacopf_any_negative_control_pose_success": bool(
                any(v["aacopf_pose_success_fraction"] > 0.0 for v in summary_geometry.values())
            ),
        },
    }
    raw = {
        "design": {
            "purpose": "P1F-E rank-deficient geometry negative controls",
            "geometries": list(geometries),
            "n_seeds": n_seeds,
            "n_particles": n_particles,
            "same_initial_cloud_across_geometries": True,
            "same_additive_range_noise_across_geometries": True,
            "quick": args.quick,
        },
        "observability": observability,
        "negative_controls": raw_geometry,
        "summary": summary,
    }

    out_dir = ROOT / "results/phase1/p1f_e"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
