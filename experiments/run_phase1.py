import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mt_ag.imu import mechanize_planar, simulate_imu_measurements
from mt_ag.metrics import summary_metrics
from mt_ag.observability import singular_values_planar_imu
from mt_ag.particle_filter import run_imu_bootstrap_pf
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory, generate_curved_trajectory


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1.yaml").read_text())
    seed = cfg["seed"]
    rng_sensor = np.random.default_rng(seed)
    rng_pf = np.random.default_rng(seed + 1)
    trajectory = generate_curved_trajectory(**cfg["simulation"])
    ideal_dr = mechanize_planar(trajectory.state[0], trajectory.ideal_imu, trajectory.dt)
    imu_data = simulate_imu_measurements(trajectory.ideal_imu, trajectory.dt, rng_sensor, **cfg["imu"])
    noisy_dr = mechanize_planar(trajectory.state[0], imu_data.measured, trajectory.dt)
    auxiliary_moving = auxiliary_trajectory(trajectory.t, "moving")
    uwb_sigma = cfg["uwb"]["sigma_range_m"]
    ranges = generate_uwb_ranges(trajectory.state[:, :2], auxiliary_moving, rng_sensor, uwb_sigma)
    pf_cfg = cfg["particle_filter"]
    pf = run_imu_bootstrap_pf(
        imu_data.measured, ranges, auxiliary_moving, trajectory.state[0], pf_cfg["initial_std"],
        trajectory.dt, rng_pf, n_particles=pf_cfg["n_particles"],
        sigma_process_accel_mps2=pf_cfg["sigma_process_accel_mps2"],
        sigma_process_gyro_rps=pf_cfg["sigma_process_gyro_rps"], sigma_uwb_m=uwb_sigma,
        resample_fraction=pf_cfg["resample_fraction"],
    )
    auxiliary_stationary = auxiliary_trajectory(trajectory.t, "stationary")
    s_stationary = singular_values_planar_imu(trajectory.state, trajectory.ideal_imu, auxiliary_stationary, trajectory.dt)
    s_moving = singular_values_planar_imu(trajectory.state, trajectory.ideal_imu, auxiliary_moving, trajectory.dt)
    metrics = {
        "ideal_imu_dr": summary_metrics(ideal_dr, trajectory.state),
        "noisy_imu_dr": summary_metrics(noisy_dr, trajectory.state),
        "pf_uwb": summary_metrics(pf.estimate, trajectory.state),
        "observability": {
            "state": ["x", "y", "v_x", "v_y", "psi"],
            "stationary_singular_values": [float(x) for x in s_stationary],
            "moving_singular_values": [float(x) for x in s_moving],
            "stationary_condition_ratio": float(s_stationary[-1] / s_stationary[0]),
            "moving_condition_ratio": float(s_moving[-1] / s_moving[0]),
        },
    }
    out = ROOT / "results/phase1"
    out.mkdir(parents=True, exist_ok=True)
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2))
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
