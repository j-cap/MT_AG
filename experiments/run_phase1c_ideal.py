from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.geometry import wrap_angle
from mt_ag.initialization import initialize_range_conditioned_particles
from mt_ag.particle_filter import run_imu_bootstrap_pf_from_particles
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory, generate_curved_trajectory

ROOT = Path(__file__).resolve().parents[1]


def terminal_convergence_time(mask, dt, hold_steps):
    mask = np.asarray(mask, dtype=bool)
    false_indices = np.flatnonzero(~mask)
    start = 0 if len(false_indices) == 0 else int(false_indices[-1] + 1)
    if len(mask) - start < hold_steps:
        return None
    return float(start * dt)


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1c.yaml").read_text())
    ideal_cfg = cfg["ideal_sanity"]
    conv_cfg = cfg["convergence"]
    seed_cfg = cfg["randomness"]

    trajectory = generate_curved_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    ranges = generate_uwb_ranges(trajectory.state[:, :2], auxiliary)
    hold_steps = round(conv_cfg["hold_time_s"] / trajectory.dt)
    late_steps = round(conv_cfg["late_window_s"] / trajectory.dt)
    yaw_threshold = np.deg2rad(conv_cfg["yaw_threshold_deg"])

    runs = []
    for seed in range(ideal_cfg["n_seeds"]):
        rng_init = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        particles = initialize_range_conditioned_particles(
            ranges[0],
            auxiliary[0],
            ideal_cfg["n_bearing"],
            ideal_cfg["n_yaw"],
            rng_init,
            radial_std_m=0.0,
            yaw_mode="uniform",
            velocity_mode="aligned_fixed_speed",
            speed_mean_mps=0.75,
        )
        rng_pf = np.random.default_rng(seed_cfg["pf_seed_offset"] + seed)
        start = time.perf_counter()
        pf = run_imu_bootstrap_pf_from_particles(
            trajectory.ideal_imu,
            ranges,
            auxiliary,
            particles,
            trajectory.dt,
            rng_pf,
            sigma_process_accel_mps2=0.0,
            sigma_process_gyro_rps=0.0,
            sigma_uwb_m=ideal_cfg["uwb_likelihood_sigma_m"],
            resample_fraction=0.0,
            initial_particles_conditioned_on_z0=True,
        )
        runtime = time.perf_counter() - start

        pos_error = np.linalg.norm(pf.estimate[:, :2] - trajectory.state[:, :2], axis=1)
        yaw_error = np.abs(wrap_angle(pf.estimate[:, 4] - trajectory.state[:, 4]))
        pos_ok = pos_error < conv_cfg["position_threshold_m"]
        pose_ok = pos_ok & (yaw_error < yaw_threshold)
        late = slice(max(0, len(pos_error) - late_steps), None)
        runs.append(
            {
                "seed": seed,
                "n_particles": len(particles),
                "position_convergence_s": terminal_convergence_time(
                    pos_ok, trajectory.dt, hold_steps
                ),
                "pose_convergence_s": terminal_convergence_time(
                    pose_ok, trajectory.dt, hold_steps
                ),
                "position_rmse_m": float(np.sqrt(np.mean(pos_error**2))),
                "late_position_rmse_m": float(
                    np.sqrt(np.mean(pos_error[late] ** 2))
                ),
                "final_position_error_m": float(pos_error[-1]),
                "late_heading_rmse_deg": float(
                    np.rad2deg(np.sqrt(np.mean(yaw_error[late] ** 2)))
                ),
                "final_heading_error_deg": float(np.rad2deg(yaw_error[-1])),
                "final_position_spread_m": float(pf.position_spread[-1]),
                "final_yaw_resultant": float(pf.yaw_resultant[-1]),
                "runtime_s": float(runtime),
            }
        )

    position_success = [r["position_convergence_s"] is not None for r in runs]
    pose_success = [r["pose_convergence_s"] is not None for r in runs]
    summary = {
        "design": {
            "n_seeds": ideal_cfg["n_seeds"],
            "n_bearing": ideal_cfg["n_bearing"],
            "n_yaw": ideal_cfg["n_yaw"],
            "n_particles": ideal_cfg["n_bearing"] * ideal_cfg["n_yaw"],
            "exact_imu": True,
            "exact_uwb": True,
            "initial_radial_std_m": 0.0,
            "process_noise": 0.0,
            "resampling": False,
            "likelihood_sigma_m": ideal_cfg["uwb_likelihood_sigma_m"],
        },
        "position_success_fraction": float(np.mean(position_success)),
        "pose_success_fraction": float(np.mean(pose_success)),
        "position_rmse_mean_m": float(np.mean([r["position_rmse_m"] for r in runs])),
        "late_position_rmse_mean_m": float(
            np.mean([r["late_position_rmse_m"] for r in runs])
        ),
        "final_position_error_mean_m": float(
            np.mean([r["final_position_error_m"] for r in runs])
        ),
        "late_heading_rmse_mean_deg": float(
            np.mean([r["late_heading_rmse_deg"] for r in runs])
        ),
        "final_heading_error_mean_deg": float(
            np.mean([r["final_heading_error_deg"] for r in runs])
        ),
        "runs": runs,
    }

    out = ROOT / "results/phase1/p1c"
    out.mkdir(parents=True, exist_ok=True)
    (out / "ideal_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
