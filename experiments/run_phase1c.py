from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.geometry import wrap_angle
from mt_ag.imu import simulate_imu_measurements
from mt_ag.initialization import initialize_range_conditioned_particles
from mt_ag.particle_filter import run_imu_bootstrap_pf_from_particles
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory, generate_curved_trajectory

ROOT = Path(__file__).resolve().parents[1]


def distribution(values):
    x = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "std": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "q95": float(np.quantile(x, 0.95)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def convergence_time(mask, dt, hold_steps):
    run = 0
    for k, ok in enumerate(mask):
        run = run + 1 if ok else 0
        if run >= hold_steps:
            return float((k - hold_steps + 1) * dt)
    return None


def aggregate(runs):
    pos_times = [r["position_convergence_s"] for r in runs if r["position_convergence_s"] is not None]
    pose_times = [r["pose_convergence_s"] for r in runs if r["pose_convergence_s"] is not None]
    return {
        "n": len(runs),
        "n_particles": runs[0]["n_particles"],
        "position_success_fraction": float(np.mean([r["position_convergence_s"] is not None for r in runs])),
        "pose_success_fraction": float(np.mean([r["pose_convergence_s"] is not None for r in runs])),
        "position_convergence_s": distribution(pos_times) if pos_times else None,
        "pose_convergence_s": distribution(pose_times) if pose_times else None,
        "position_rmse_m": distribution([r["position_rmse_m"] for r in runs]),
        "late_position_rmse_m": distribution([r["late_position_rmse_m"] for r in runs]),
        "final_position_error_m": distribution([r["final_position_error_m"] for r in runs]),
        "late_heading_rmse_deg": distribution([r["late_heading_rmse_deg"] for r in runs]),
        "final_heading_error_deg": distribution([r["final_heading_error_deg"] for r in runs]),
        "runtime_s": distribution([r["runtime_s"] for r in runs]),
        "final_position_spread_m": distribution([r["final_position_spread_m"] for r in runs]),
        "final_yaw_resultant": distribution([r["final_yaw_resultant"] for r in runs]),
    }


def compact(result):
    return {
        "n": result["n"],
        "n_particles": result["n_particles"],
        "position_success_fraction": result["position_success_fraction"],
        "pose_success_fraction": result["pose_success_fraction"],
        "position_convergence_median_s": (
            None if result["position_convergence_s"] is None else result["position_convergence_s"]["median"]
        ),
        "pose_convergence_median_s": (
            None if result["pose_convergence_s"] is None else result["pose_convergence_s"]["median"]
        ),
        "position_rmse_mean_m": result["position_rmse_m"]["mean"],
        "late_position_rmse_mean_m": result["late_position_rmse_m"]["mean"],
        "final_position_error_mean_m": result["final_position_error_m"]["mean"],
        "late_heading_rmse_mean_deg": result["late_heading_rmse_deg"]["mean"],
        "final_heading_error_mean_deg": result["final_heading_error_deg"]["mean"],
        "runtime_mean_s": result["runtime_s"]["mean"],
        "final_position_spread_mean_m": result["final_position_spread_m"]["mean"],
        "final_yaw_resultant_mean": result["final_yaw_resultant"]["mean"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1c.yaml").read_text())
    trajectory = generate_curved_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    seed_cfg = cfg["randomness"]
    uwb_cfg = cfg["uwb"]
    pf_cfg = cfg["particle_filter"]
    conv_cfg = cfg["convergence"]
    hold_steps = round(conv_cfg["hold_time_s"] / trajectory.dt)
    late_steps = round(conv_cfg["late_window_s"] / trajectory.dt)
    yaw_threshold = np.deg2rad(conv_cfg["yaw_threshold_deg"])
    sensor_cache = {}

    def seed_count(section):
        n = cfg[section]["n_seeds"]
        return min(3, n) if args.quick else n

    def measurements(seed):
        if seed not in sensor_cache:
            rng = np.random.default_rng(seed_cfg["sensor_seed_offset"] + seed)
            imu = simulate_imu_measurements(trajectory.ideal_imu, trajectory.dt, rng, **cfg["imu"])
            ranges = generate_uwb_ranges(
                trajectory.state[:, :2], auxiliary, rng, uwb_cfg["sigma_range_m"]
            )
            sensor_cache[seed] = imu, ranges
        return sensor_cache[seed]

    def run_case(seed, initializer, n_bearing, n_yaw):
        imu, ranges = measurements(seed)
        init = dict(initializer)
        init.pop("n_seeds", None)
        init.pop("n_bearing", None)
        init.pop("n_yaw", None)
        yaw_mode = init.pop("yaw_mode")
        velocity_mode = init.pop("velocity_mode")
        if yaw_mode == "known":
            init["known_yaw_rad"] = trajectory.state[0, 4]

        rng_init = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        particles = initialize_range_conditioned_particles(
            ranges[0],
            auxiliary[0],
            n_bearing,
            n_yaw,
            rng_init,
            radial_std_m=uwb_cfg["initial_radial_std_m"],
            yaw_mode=yaw_mode,
            velocity_mode=velocity_mode,
            **init,
        )
        rng_pf = np.random.default_rng(seed_cfg["pf_seed_offset"] + seed)
        start = time.perf_counter()
        pf = run_imu_bootstrap_pf_from_particles(
            imu.measured,
            ranges,
            auxiliary,
            particles,
            trajectory.dt,
            rng_pf,
            sigma_process_accel_mps2=pf_cfg["sigma_process_accel_mps2"],
            sigma_process_gyro_rps=pf_cfg["sigma_process_gyro_rps"],
            sigma_uwb_m=uwb_cfg["sigma_range_m"],
            resample_fraction=pf_cfg["resample_fraction"],
            initial_particles_conditioned_on_z0=True,
        )
        runtime = time.perf_counter() - start

        pos_error = np.linalg.norm(pf.estimate[:, :2] - trajectory.state[:, :2], axis=1)
        yaw_error = np.abs(wrap_angle(pf.estimate[:, 4] - trajectory.state[:, 4]))
        pos_ok = pos_error < conv_cfg["position_threshold_m"]
        pose_ok = pos_ok & (yaw_error < yaw_threshold)
        late = slice(max(0, len(pos_error) - late_steps), None)
        return {
            "seed": seed,
            "n_particles": len(particles),
            "position_convergence_s": convergence_time(pos_ok, trajectory.dt, hold_steps),
            "pose_convergence_s": convergence_time(pose_ok, trajectory.dt, hold_steps),
            "position_rmse_m": float(np.sqrt(np.mean(pos_error**2))),
            "late_position_rmse_m": float(np.sqrt(np.mean(pos_error[late] ** 2))),
            "final_position_error_m": float(pos_error[-1]),
            "late_heading_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error[late] ** 2)))),
            "final_heading_error_deg": float(np.rad2deg(yaw_error[-1])),
            "runtime_s": float(runtime),
            "final_position_spread_m": float(pf.position_spread[-1]),
            "final_yaw_resultant": float(pf.yaw_resultant[-1]),
        }

    raw = {
        "design": {
            "trajectory": "P1A/P1B deterministic 60 s curved trajectory",
            "auxiliary": "moving and globally known",
            "uwb": "Gaussian sigma 0.12 m, synchronous 10 Hz",
            "first_range_handling": "initial ring conditioned on z0; filter likelihood starts at z1",
            "position_threshold_m": conv_cfg["position_threshold_m"],
            "yaw_threshold_deg": conv_cfg["yaw_threshold_deg"],
            "hold_time_s": conv_cfg["hold_time_s"],
            "quick": args.quick,
        },
        "known_yaw_control": {},
        "core_unknown_pose": {},
        "velocity_prior": {},
        "grid_resolution": {},
    }

    control = cfg["known_yaw_control"]
    control_runs = [
        run_case(seed, control, 100 if args.quick else control["n_bearing"], control["n_yaw"])
        for seed in range(seed_count("known_yaw_control"))
    ]
    raw["known_yaw_control"] = {"aggregate": aggregate(control_runs), "runs": control_runs}

    core = cfg["core_unknown_pose"]
    core_grid = 40 if args.quick else core["n_bearing"]
    core_runs = [
        run_case(seed, core, core_grid, core_grid)
        for seed in range(seed_count("core_unknown_pose"))
    ]
    raw["core_unknown_pose"] = {"aggregate": aggregate(core_runs), "runs": core_runs}

    velocity = cfg["velocity_prior"]
    velocity_grid = 40 if args.quick else velocity["n_bearing"]
    for name, case in velocity["cases"].items():
        initializer = {"yaw_mode": "uniform", **case}
        runs = [
            run_case(seed, initializer, velocity_grid, velocity_grid)
            for seed in range(seed_count("velocity_prior"))
        ]
        raw["velocity_prior"][name] = {
            "configuration": case,
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    grids = [25, 40] if args.quick else cfg["grid_resolution"]["values"]
    grid_initializer = {
        "yaw_mode": "uniform",
        "velocity_mode": "aligned_fixed_speed",
        "speed_mean_mps": 0.75,
    }
    for grid in grids:
        runs = [
            run_case(seed, grid_initializer, grid, grid)
            for seed in range(seed_count("grid_resolution"))
        ]
        raw["grid_resolution"][str(grid)] = {
            "n_bearing": grid,
            "n_yaw": grid,
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    summary = {
        "known_yaw_control": compact(raw["known_yaw_control"]["aggregate"]),
        "core_unknown_pose": compact(raw["core_unknown_pose"]["aggregate"]),
        "velocity_prior": {
            name: compact(value["aggregate"]) for name, value in raw["velocity_prior"].items()
        },
        "grid_resolution": {
            name: compact(value["aggregate"]) for name, value in raw["grid_resolution"].items()
        },
    }
    out = ROOT / "results/phase1/p1c"
    out.mkdir(parents=True, exist_ok=True)
    (out / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
