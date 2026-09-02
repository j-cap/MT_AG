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
        "q05": float(np.quantile(x, 0.05)),
        "q95": float(np.quantile(x, 0.95)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def sustained_convergence_time(mask, dt, hold_steps):
    run = 0
    for k, value in enumerate(mask):
        run = run + 1 if value else 0
        if run >= hold_steps:
            return float((k - hold_steps + 1) * dt)
    return None


def aggregate(runs):
    successful_position_times = [
        r["position_convergence_s"]
        for r in runs
        if r["position_convergence_s"] is not None
    ]
    successful_pose_times = [
        r["pose_convergence_s"]
        for r in runs
        if r["pose_convergence_s"] is not None
    ]
    return {
        "n": len(runs),
        "n_particles": int(runs[0]["n_particles"]),
        "position_success_fraction": float(
            np.mean([r["position_convergence_s"] is not None for r in runs])
        ),
        "pose_success_fraction": float(
            np.mean([r["pose_convergence_s"] is not None for r in runs])
        ),
        "position_convergence_s": (
            distribution(successful_position_times) if successful_position_times else None
        ),
        "pose_convergence_s": (
            distribution(successful_pose_times) if successful_pose_times else None
        ),
        "position_rmse_m": distribution([r["position_rmse_m"] for r in runs]),
        "late_position_rmse_m": distribution(
            [r["late_position_rmse_m"] for r in runs]
        ),
        "final_position_error_m": distribution(
            [r["final_position_error_m"] for r in runs]
        ),
        "heading_rmse_deg": distribution([r["heading_rmse_deg"] for r in runs]),
        "late_heading_rmse_deg": distribution(
            [r["late_heading_rmse_deg"] for r in runs]
        ),
        "final_heading_error_deg": distribution(
            [r["final_heading_error_deg"] for r in runs]
        ),
        "runtime_s": distribution([r["runtime_s"] for r in runs]),
        "final_position_spread_m": distribution(
            [r["final_position_spread_m"] for r in runs]
        ),
        "final_yaw_resultant": distribution(
            [r["final_yaw_resultant"] for r in runs]
        ),
    }


def compact(a):
    return {
        "n": a["n"],
        "n_particles": a["n_particles"],
        "position_success_fraction": a["position_success_fraction"],
        "pose_success_fraction": a["pose_success_fraction"],
        "position_convergence_median_s": (
            None
            if a["position_convergence_s"] is None
            else a["position_convergence_s"]["median"]
        ),
        "pose_convergence_median_s": (
            None
            if a["pose_convergence_s"] is None
            else a["pose_convergence_s"]["median"]
        ),
        "position_rmse_mean_m": a["position_rmse_m"]["mean"],
        "late_position_rmse_mean_m": a["late_position_rmse_m"]["mean"],
        "final_position_error_mean_m": a["final_position_error_m"]["mean"],
        "late_heading_rmse_mean_deg": a["late_heading_rmse_deg"]["mean"],
        "final_heading_error_mean_deg": a["final_heading_error_deg"]["mean"],
        "runtime_mean_s": a["runtime_s"]["mean"],
        "final_position_spread_mean_m": a["final_position_spread_m"]["mean"],
        "final_yaw_resultant_mean": a["final_yaw_resultant"]["mean"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run at most three seeds per condition and smaller grids.",
    )
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1c.yaml").read_text())
    trajectory = generate_curved_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    seed_cfg = cfg["randomness"]
    imu_cfg = cfg["imu"]
    uwb_cfg = cfg["uwb"]
    pf_cfg = cfg["particle_filter"]
    conv_cfg = cfg["convergence"]
    hold_steps = int(round(conv_cfg["hold_time_s"] / trajectory.dt))
    late_steps = int(round(conv_cfg["late_window_s"] / trajectory.dt))
    yaw_threshold = np.deg2rad(conv_cfg["yaw_threshold_deg"])

    def n_seeds(section):
        value = int(cfg[section]["n_seeds"])
        return min(3, value) if args.quick else value

    sensor_cache = {}

    def sensor_data(seed):
        if seed not in sensor_cache:
            rng_sensor = np.random.default_rng(seed_cfg["sensor_seed_offset"] + seed)
            imu = simulate_imu_measurements(
                trajectory.ideal_imu,
                trajectory.dt,
                rng_sensor,
                **imu_cfg,
            )
            ranges = generate_uwb_ranges(
                trajectory.state[:, :2],
                auxiliary,
                rng_sensor,
                float(uwb_cfg["sigma_range_m"]),
            )
            sensor_cache[seed] = (imu, ranges)
        return sensor_cache[seed]

    def run_case(seed, init_cfg, *, n_bearing, n_yaw):
        imu, ranges = sensor_data(seed)
        rng_init = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        init_kwargs = dict(init_cfg)
        yaw_mode = init_kwargs.pop("yaw_mode")
        velocity_mode = init_kwargs.pop("velocity_mode")
        init_kwargs.pop("n_seeds", None)
        init_kwargs.pop("n_bearing", None)
        init_kwargs.pop("n_yaw", None)
        if yaw_mode == "known":
            init_kwargs["known_yaw_rad"] = float(trajectory.state[0, 4])
        particles = initialize_range_conditioned_particles(
            float(ranges[0]),
            auxiliary[0],
            int(n_bearing),
            int(n_yaw),
            rng_init,
            radial_std_m=float(uwb_cfg["initial_radial_std_m"]),
            yaw_mode=yaw_mode,
            velocity_mode=velocity_mode,
            **init_kwargs,
        )
        rng_pf = np.random.default_rng(seed_cfg["pf_seed_offset"] + seed)
        start = time.perf_counter()
        result = run_imu_bootstrap_pf_from_particles(
            imu.measured,
            ranges,
            auxiliary,
            particles,
            trajectory.dt,
            rng_pf,
            sigma_process_accel_mps2=float(pf_cfg["sigma_process_accel_mps2"]),
            sigma_process_gyro_rps=float(pf_cfg["sigma_process_gyro_rps"]),
            sigma_uwb_m=float(uwb_cfg["sigma_range_m"]),
            resample_fraction=float(pf_cfg["resample_fraction"]),
            initial_particles_conditioned_on_z0=True,
        )
        runtime = time.perf_counter() - start

        pos_error = np.linalg.norm(
            result.estimate[:, :2] - trajectory.state[:, :2], axis=1
        )
        yaw_error = np.abs(
            wrap_angle(result.estimate[:, 4] - trajectory.state[:, 4])
        )
        position_mask = pos_error < float(conv_cfg["position_threshold_m"])
        pose_mask = position_mask & (yaw_error < yaw_threshold)
        position_convergence = sustained_convergence_time(
            position_mask, trajectory.dt, hold_steps
        )
        pose_convergence = sustained_convergence_time(
            pose_mask, trajectory.dt, hold_steps
        )
        late = slice(max(0, len(pos_error) - late_steps), len(pos_error))
        return {
            "seed": int(seed),
            "n_particles": int(len(particles)),
            "position_convergence_s": position_convergence,
            "pose_convergence_s": pose_convergence,
            "position_rmse_m": float(np.sqrt(np.mean(pos_error**2))),
            "late_position_rmse_m": float(
                np.sqrt(np.mean(pos_error[late] ** 2))
            ),
            "final_position_error_m": float(pos_error[-1]),
            "heading_rmse_deg": float(
                np.rad2deg(np.sqrt(np.mean(yaw_error**2)))
            ),
            "late_heading_rmse_deg": float(
                np.rad2deg(np.sqrt(np.mean(yaw_error[late] ** 2)))
            ),
            "final_heading_error_deg": float(np.rad2deg(yaw_error[-1])),
            "runtime_s": float(runtime),
            "final_position_spread_m": float(result.position_spread[-1]),
            "final_yaw_resultant": float(result.yaw_resultant[-1]),
            "initial_range_measurement_m": float(ranges[0]),
        }

    raw = {
        "design": {
            "trajectory": "P1A/P1B deterministic 60 s curved trajectory",
            "auxiliary": "moving, globally known",
            "uwb": "Gaussian sigma 0.12 m, synchronous 10 Hz",
            "initialization": (
                "position on first-range ring; z0 not reused in PF likelihood"
            ),
            "position_convergence_threshold_m": conv_cfg["position_threshold_m"],
            "yaw_convergence_threshold_deg": conv_cfg["yaw_threshold_deg"],
            "hold_time_s": conv_cfg["hold_time_s"],
            "quick_mode": bool(args.quick),
        },
        "known_yaw_control": {},
        "core_unknown_pose": {},
        "velocity_prior": {},
        "grid_resolution": {},
    }

    control_cfg = dict(cfg["known_yaw_control"])
    control_runs = [
        run_case(
            seed,
            control_cfg,
            n_bearing=(100 if args.quick else control_cfg["n_bearing"]),
            n_yaw=control_cfg["n_yaw"],
        )
        for seed in range(n_seeds("known_yaw_control"))
    ]
    raw["known_yaw_control"] = {
        "aggregate": aggregate(control_runs),
        "runs": control_runs,
    }

    core_cfg = dict(cfg["core_unknown_pose"])
    core_runs = [
        run_case(
            seed,
            core_cfg,
            n_bearing=(40 if args.quick else core_cfg["n_bearing"]),
            n_yaw=(40 if args.quick else core_cfg["n_yaw"]),
        )
        for seed in range(n_seeds("core_unknown_pose"))
    ]
    raw["core_unknown_pose"] = {
        "aggregate": aggregate(core_runs),
        "runs": core_runs,
    }

    velocity_cfg = cfg["velocity_prior"]
    for name, case in velocity_cfg["cases"].items():
        init_cfg = {"yaw_mode": "uniform", **case}
        runs = [
            run_case(
                seed,
                init_cfg,
                n_bearing=(40 if args.quick else velocity_cfg["n_bearing"]),
                n_yaw=(40 if args.quick else velocity_cfg["n_yaw"]),
            )
            for seed in range(n_seeds("velocity_prior"))
        ]
        raw["velocity_prior"][name] = {
            "configuration": case,
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    grid_cfg = cfg["grid_resolution"]
    values = [25, 40] if args.quick else grid_cfg["values"]
    for grid_size in values:
        init_cfg = {
            "yaw_mode": "uniform",
            "velocity_mode": "aligned_fixed_speed",
            "speed_mean_mps": 0.75,
        }
        runs = [
            run_case(
                seed,
                init_cfg,
                n_bearing=int(grid_size),
                n_yaw=int(grid_size),
            )
            for seed in range(n_seeds("grid_resolution"))
        ]
        raw["grid_resolution"][str(grid_size)] = {
            "n_bearing": int(grid_size),
            "n_yaw": int(grid_size),
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    summary = {
        "known_yaw_control": compact(raw["known_yaw_control"]["aggregate"]),
        "core_unknown_pose": compact(raw["core_unknown_pose"]["aggregate"]),
        "velocity_prior": {
            name: compact(value["aggregate"])
            for name, value in raw["velocity_prior"].items()
        },
        "grid_resolution": {
            name: compact(value["aggregate"])
            for name, value in raw["grid_resolution"].items()
        },
    }

    out = ROOT / "results/phase1/p1c"
    out.mkdir(parents=True, exist_ok=True)
    (out / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
