from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.geometry import wrap_angle
from mt_ag.imu import simulate_imu_measurements
from mt_ag.initialization import initialize_range_conditioned_particles
from mt_ag.observability import (
    observability_history_planar_imu,
    relative_geometry_metrics,
    singular_values_planar_imu,
)
from mt_ag.particle_filter import run_imu_bootstrap_pf_from_particles
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory, generate_curved_trajectory

ROOT = Path(__file__).resolve().parents[1]


def distribution(values):
    values = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        "q95": float(np.quantile(values, 0.95)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def terminal_convergence_time(mask, dt, hold_steps):
    mask = np.asarray(mask, dtype=bool)
    false_indices = np.flatnonzero(~mask)
    start = 0 if len(false_indices) == 0 else int(false_indices[-1] + 1)
    if len(mask) - start < hold_steps:
        return None
    return float(start * dt)


def aggregate_runs(runs):
    position_times = [
        run["position_convergence_s"]
        for run in runs
        if run["position_convergence_s"] is not None
    ]
    pose_times = [
        run["pose_convergence_s"]
        for run in runs
        if run["pose_convergence_s"] is not None
    ]
    return {
        "n": len(runs),
        "position_success_fraction": float(
            np.mean([run["position_convergence_s"] is not None for run in runs])
        ),
        "pose_success_fraction": float(
            np.mean([run["pose_convergence_s"] is not None for run in runs])
        ),
        "position_convergence_s": (
            distribution(position_times) if position_times else None
        ),
        "pose_convergence_s": distribution(pose_times) if pose_times else None,
        "position_rmse_m": distribution([run["position_rmse_m"] for run in runs]),
        "late_position_rmse_m": distribution(
            [run["late_position_rmse_m"] for run in runs]
        ),
        "final_position_error_m": distribution(
            [run["final_position_error_m"] for run in runs]
        ),
        "late_heading_rmse_deg": distribution(
            [run["late_heading_rmse_deg"] for run in runs]
        ),
        "final_heading_error_deg": distribution(
            [run["final_heading_error_deg"] for run in runs]
        ),
        "final_position_spread_m": distribution(
            [run["final_position_spread_m"] for run in runs]
        ),
        "final_yaw_resultant": distribution(
            [run["final_yaw_resultant"] for run in runs]
        ),
        "runtime_s": distribution([run["runtime_s"] for run in runs]),
    }


def compact(aggregate):
    return {
        "n": aggregate["n"],
        "position_success_fraction": aggregate["position_success_fraction"],
        "pose_success_fraction": aggregate["pose_success_fraction"],
        "position_convergence_median_s": (
            None
            if aggregate["position_convergence_s"] is None
            else aggregate["position_convergence_s"]["median"]
        ),
        "pose_convergence_median_s": (
            None
            if aggregate["pose_convergence_s"] is None
            else aggregate["pose_convergence_s"]["median"]
        ),
        "position_rmse_mean_m": aggregate["position_rmse_m"]["mean"],
        "late_position_rmse_mean_m": aggregate["late_position_rmse_m"]["mean"],
        "final_position_error_mean_m": aggregate["final_position_error_m"]["mean"],
        "late_heading_rmse_mean_deg": aggregate["late_heading_rmse_deg"]["mean"],
        "final_heading_error_mean_deg": aggregate["final_heading_error_deg"]["mean"],
        "final_position_spread_mean_m": aggregate["final_position_spread_m"]["mean"],
        "final_yaw_resultant_mean": aggregate["final_yaw_resultant"]["mean"],
        "runtime_mean_s": aggregate["runtime_s"]["mean"],
    }


def first_observable_time(history, dt, ratio_threshold):
    mask = (history["rank"] == 5) & (history["sigma_ratio"] >= ratio_threshold)
    indices = np.flatnonzero(mask)
    return None if len(indices) == 0 else float(indices[0] * dt)


def build_auxiliary(trajectory, mode):
    return auxiliary_trajectory(
        trajectory.t,
        mode,
        target_positions=trajectory.state[:, :2],
    )


def write_snapshot(path, snapshot):
    particles = snapshot["particles"]
    weights = snapshot["weights"]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["p_x_m", "p_y_m", "v_x_mps", "v_y_mps", "yaw_rad", "weight"])
        for particle, weight in zip(particles, weights, strict=True):
            writer.writerow([*particle, weight])


def write_summary_markdown(path, summary):
    lines = [
        "# Phase P1D — Geometry and observability",
        "",
        "## Core geometry comparison",
        "",
        "| Geometry | Bearing span | Final obs. ratio | Realistic pose success | Ideal pose success | Realistic late position RMSE |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for mode, result in summary["cases"].items():
        geometry = result["geometry"]
        realistic = result["realistic"]
        ideal = result["idealized"]
        lines.append(
            f"| {mode} | {geometry['bearing_span_deg']:.2f} deg | "
            f"{result['observability']['final_sigma_ratio']:.3e} | "
            f"{realistic['pose_success_fraction']:.0%} | "
            f"{ideal['pose_success_fraction']:.0%} | "
            f"{realistic['late_position_rmse_mean_m']:.3f} m |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            summary["interpretation"],
            "",
            "The observability ratio is a local finite-horizon diagnostic, not a complete nonlinear observability proof or an application threshold.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1d.yaml").read_text())
    trajectory = generate_curved_trajectory(**cfg["simulation"])
    dt = trajectory.dt
    out = ROOT / "results/phase1/p1d"
    out.mkdir(parents=True, exist_ok=True)

    conv_cfg = cfg["convergence"]
    hold_steps = round(conv_cfg["hold_time_s"] / dt)
    late_steps = round(conv_cfg["late_window_s"] / dt)
    yaw_threshold = np.deg2rad(conv_cfg["yaw_threshold_deg"])
    obs_cfg = cfg["observability"]
    init_cfg = cfg["initialization"]
    pf_cfg = cfg["particle_filter"]
    seed_cfg = cfg["randomness"]
    representative_seed = cfg["realistic"]["representative_seed"]
    snapshot_steps = [
        round(value / dt) for value in cfg["realistic"]["snapshot_times_s"]
    ]

    cases = {}
    diagnostic_rows = []

    for mode in cfg["geometry_cases"]:
        auxiliary = build_auxiliary(trajectory, mode)
        geometry = relative_geometry_metrics(trajectory.state[:, :2], auxiliary)
        history = observability_history_planar_imu(
            trajectory.state,
            trajectory.ideal_imu,
            auxiliary,
            dt,
            rank_rtol=obs_cfg["rank_rtol"],
        )
        final_singular = singular_values_planar_imu(
            trajectory.state, trajectory.ideal_imu, auxiliary, dt
        )
        observability = {
            "final_singular_values": [float(value) for value in final_singular],
            "final_rank": int(history["rank"][-1]),
            "final_sigma_ratio": float(history["sigma_ratio"][-1]),
            "first_ratio_threshold_time_s": first_observable_time(
                history, dt, obs_cfg["sigma_ratio_threshold"]
            ),
            "ratio_threshold": float(obs_cfg["sigma_ratio_threshold"]),
        }

        realistic_runs = []
        time_errors = []
        time_spreads = []
        time_yaw_resultants = []
        time_neff_fractions = []

        for seed in range(cfg["realistic"]["n_seeds"]):
            rng_sensor = np.random.default_rng(seed_cfg["sensor_seed_offset"] + seed)
            imu = simulate_imu_measurements(
                trajectory.ideal_imu,
                dt,
                rng_sensor,
                **cfg["imu"],
            )
            ranges = generate_uwb_ranges(
                trajectory.state[:, :2],
                auxiliary,
                rng_sensor,
                cfg["uwb"]["sigma_range_m"],
            )
            rng_init = np.random.default_rng(
                seed_cfg["initializer_seed_offset"] + seed
            )
            particles = initialize_range_conditioned_particles(
                ranges[0],
                auxiliary[0],
                init_cfg["n_bearing"],
                init_cfg["n_yaw"],
                rng_init,
                radial_std_m=cfg["uwb"]["initial_radial_std_m"],
                yaw_mode=init_cfg["yaw_mode"],
                velocity_mode=init_cfg["velocity_mode"],
                speed_mean_mps=init_cfg["speed_mean_mps"],
            )
            rng_pf = np.random.default_rng(seed_cfg["pf_seed_offset"] + seed)
            start = time.perf_counter()
            pf = run_imu_bootstrap_pf_from_particles(
                imu.measured,
                ranges,
                auxiliary,
                particles,
                dt,
                rng_pf,
                sigma_process_accel_mps2=pf_cfg["sigma_process_accel_mps2"],
                sigma_process_gyro_rps=pf_cfg["sigma_process_gyro_rps"],
                sigma_uwb_m=cfg["uwb"]["sigma_range_m"],
                resample_fraction=pf_cfg["resample_fraction"],
                initial_particles_conditioned_on_z0=True,
                snapshot_steps=snapshot_steps if seed == representative_seed else None,
            )
            runtime = time.perf_counter() - start
            position_error = np.linalg.norm(
                pf.estimate[:, :2] - trajectory.state[:, :2], axis=1
            )
            yaw_error = np.abs(wrap_angle(pf.estimate[:, 4] - trajectory.state[:, 4]))
            position_ok = position_error < conv_cfg["position_threshold_m"]
            pose_ok = position_ok & (yaw_error < yaw_threshold)
            late = slice(max(0, len(position_error) - late_steps), None)
            realistic_runs.append(
                {
                    "seed": seed,
                    "position_convergence_s": terminal_convergence_time(
                        position_ok, dt, hold_steps
                    ),
                    "pose_convergence_s": terminal_convergence_time(
                        pose_ok, dt, hold_steps
                    ),
                    "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
                    "late_position_rmse_m": float(
                        np.sqrt(np.mean(position_error[late] ** 2))
                    ),
                    "final_position_error_m": float(position_error[-1]),
                    "late_heading_rmse_deg": float(
                        np.rad2deg(np.sqrt(np.mean(yaw_error[late] ** 2)))
                    ),
                    "final_heading_error_deg": float(np.rad2deg(yaw_error[-1])),
                    "final_position_spread_m": float(pf.position_spread[-1]),
                    "final_yaw_resultant": float(pf.yaw_resultant[-1]),
                    "runtime_s": float(runtime),
                }
            )
            time_errors.append(position_error)
            time_spreads.append(pf.position_spread)
            time_yaw_resultants.append(pf.yaw_resultant)
            time_neff_fractions.append(pf.neff / len(particles))

            if seed == representative_seed:
                for step, snapshot in pf.snapshots.items():
                    write_snapshot(
                        out / f"snapshot_{mode}_seed{seed}_t{step * dt:.1f}s.csv",
                        snapshot,
                    )

        errors = np.asarray(time_errors)
        spreads = np.asarray(time_spreads)
        yaw_resultants = np.asarray(time_yaw_resultants)
        neff_fractions = np.asarray(time_neff_fractions)
        for k, t in enumerate(trajectory.t):
            diagnostic_rows.append(
                {
                    "geometry": mode,
                    "time_s": float(t),
                    "obs_sigma_min": float(history["sigma_min"][k]),
                    "obs_sigma_ratio": float(history["sigma_ratio"][k]),
                    "obs_rank": int(history["rank"][k]),
                    "position_error_mean_m": float(np.mean(errors[:, k])),
                    "position_error_median_m": float(np.median(errors[:, k])),
                    "position_spread_mean_m": float(np.mean(spreads[:, k])),
                    "yaw_resultant_mean": float(np.mean(yaw_resultants[:, k])),
                    "neff_fraction_mean": float(np.mean(neff_fractions[:, k])),
                }
            )

        ideal_runs = []
        ideal_cfg = cfg["idealized"]
        exact_ranges = generate_uwb_ranges(trajectory.state[:, :2], auxiliary)
        for seed in range(ideal_cfg["n_seeds"]):
            rng_init = np.random.default_rng(
                seed_cfg["initializer_seed_offset"] + seed
            )
            particles = initialize_range_conditioned_particles(
                exact_ranges[0],
                auxiliary[0],
                ideal_cfg["n_bearing"],
                ideal_cfg["n_yaw"],
                rng_init,
                radial_std_m=0.0,
                yaw_mode="uniform",
                velocity_mode="aligned_fixed_speed",
                speed_mean_mps=init_cfg["speed_mean_mps"],
            )
            rng_pf = np.random.default_rng(seed_cfg["pf_seed_offset"] + seed)
            start = time.perf_counter()
            pf = run_imu_bootstrap_pf_from_particles(
                trajectory.ideal_imu,
                exact_ranges,
                auxiliary,
                particles,
                dt,
                rng_pf,
                sigma_process_accel_mps2=0.0,
                sigma_process_gyro_rps=0.0,
                sigma_uwb_m=ideal_cfg["uwb_likelihood_sigma_m"],
                resample_fraction=0.0,
                initial_particles_conditioned_on_z0=True,
            )
            runtime = time.perf_counter() - start
            position_error = np.linalg.norm(
                pf.estimate[:, :2] - trajectory.state[:, :2], axis=1
            )
            yaw_error = np.abs(wrap_angle(pf.estimate[:, 4] - trajectory.state[:, 4]))
            position_ok = position_error < conv_cfg["position_threshold_m"]
            pose_ok = position_ok & (yaw_error < yaw_threshold)
            late = slice(max(0, len(position_error) - late_steps), None)
            ideal_runs.append(
                {
                    "seed": seed,
                    "position_convergence_s": terminal_convergence_time(
                        position_ok, dt, hold_steps
                    ),
                    "pose_convergence_s": terminal_convergence_time(
                        pose_ok, dt, hold_steps
                    ),
                    "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
                    "late_position_rmse_m": float(
                        np.sqrt(np.mean(position_error[late] ** 2))
                    ),
                    "final_position_error_m": float(position_error[-1]),
                    "late_heading_rmse_deg": float(
                        np.rad2deg(np.sqrt(np.mean(yaw_error[late] ** 2)))
                    ),
                    "final_heading_error_deg": float(np.rad2deg(yaw_error[-1])),
                    "final_position_spread_m": float(pf.position_spread[-1]),
                    "final_yaw_resultant": float(pf.yaw_resultant[-1]),
                    "runtime_s": float(runtime),
                }
            )

        cases[mode] = {
            "geometry": geometry,
            "observability": observability,
            "realistic": compact(aggregate_runs(realistic_runs)),
            "idealized": compact(aggregate_runs(ideal_runs)),
            "realistic_runs": realistic_runs,
            "idealized_runs": ideal_runs,
        }

    with (out / "diagnostics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(diagnostic_rows[0]))
        writer.writeheader()
        writer.writerows(diagnostic_rows)

    informative = cases["moving"]
    stationary = cases["stationary"]
    constant_bearing = cases["constant_bearing"]
    interpretation = (
        "The controlled geometry study separates range availability from geometric "
        "information. The stationary and constant-bearing cases retain poor local "
        "observability compared with the informative moving auxiliary. The idealized "
        "particle filter is used to show whether this geometry difference is reflected "
        "in global-pose recovery independently of realistic sensor uncertainty; the "
        "realistic runs then show how much of that advantage survives with the P1C "
        "noise model. Detailed numerical values are given in the tables above."
    )
    summary = {
        "design": {
            "geometry_cases": cfg["geometry_cases"],
            "realistic_n_seeds": cfg["realistic"]["n_seeds"],
            "idealized_n_seeds": cfg["idealized"]["n_seeds"],
            "realistic_particles": init_cfg["n_bearing"] * init_cfg["n_yaw"],
            "idealized_particles": (
                cfg["idealized"]["n_bearing"] * cfg["idealized"]["n_yaw"]
            ),
            "matched_noise_and_particle_randomness_across_geometries": True,
        },
        "cases": cases,
        "comparison": {
            "moving_vs_stationary_final_obs_ratio": (
                informative["observability"]["final_sigma_ratio"]
                / max(stationary["observability"]["final_sigma_ratio"], 1e-300)
            ),
            "moving_vs_constant_bearing_final_obs_ratio": (
                informative["observability"]["final_sigma_ratio"]
                / max(
                    constant_bearing["observability"]["final_sigma_ratio"], 1e-300
                )
            ),
        },
        "interpretation": interpretation,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    write_summary_markdown(out / "summary.md", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
