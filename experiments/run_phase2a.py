from __future__ import annotations

import csv
import json
import math
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.fleet_simulation import all_pairs, generate_four_node_fleet, pairwise_distances
from mt_ag.imu import simulate_imu_measurements
from mt_ag.joint_ekf import run_joint_ekf

ROOT = Path(__file__).resolve().parents[1]


def _fleet_metrics(estimate: np.ndarray, truth: np.ndarray) -> dict:
    position_error_vector = estimate[:, :, :2] - truth[:, :, :2]
    position_error = np.linalg.norm(position_error_vector, axis=2)
    yaw_error = (estimate[:, :, 4] - truth[:, :, 4] + np.pi) % (2.0 * np.pi) - np.pi
    node_rmse = np.sqrt(np.mean(position_error**2, axis=0))

    centroid_error_vector = np.mean(position_error_vector, axis=1)
    centroid_error = np.linalg.norm(centroid_error_vector, axis=1)
    estimate_centered = estimate[:, :, :2] - np.mean(estimate[:, :, :2], axis=1, keepdims=True)
    truth_centered = truth[:, :, :2] - np.mean(truth[:, :, :2], axis=1, keepdims=True)
    shape_error = np.linalg.norm(estimate_centered - truth_centered, axis=2)

    return {
        "fleet_position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
        "fleet_position_p95_m": float(np.quantile(position_error, 0.95)),
        "worst_node_position_rmse_m": float(np.max(node_rmse)),
        "centroid_position_rmse_m": float(np.sqrt(np.mean(centroid_error**2))),
        "relative_shape_rmse_m": float(np.sqrt(np.mean(shape_error**2))),
        "fleet_yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error**2)))),
        "node_position_rmse_m": [float(value) for value in node_rmse],
    }


def _initial_covariance(cfg: dict, n_nodes: int) -> np.ndarray:
    pos = float(cfg["covariance_position_std_m"])
    vel = float(cfg["covariance_velocity_std_mps"])
    yaw = np.deg2rad(float(cfg["covariance_yaw_std_deg"]))
    node_var = np.array([pos**2, pos**2, vel**2, vel**2, yaw**2])
    return np.diag(np.tile(node_var, n_nodes))


def _range_data(truth: np.ndarray, sigma: float, rng: np.random.Generator):
    n_time, n_nodes, _ = truth.shape
    ranges = np.zeros((n_time, n_nodes, n_nodes), dtype=float)
    for i, j in all_pairs(n_nodes):
        distance = np.linalg.norm(truth[:, i, :2] - truth[:, j, :2], axis=1)
        noisy = distance + rng.normal(0.0, sigma, n_time)
        ranges[:, i, j] = noisy
        ranges[:, j, i] = noisy
    return ranges


def _full_mask(n_time: int, n_nodes: int, stride: int) -> np.ndarray:
    mask = np.zeros((n_time, n_nodes, n_nodes), dtype=bool)
    for i, j in all_pairs(n_nodes):
        mask[stride::stride, i, j] = True
        mask[stride::stride, j, i] = True
    return mask


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase2a.yaml").read_text())
    sim = cfg["simulation"]
    imu_cfg = cfg["imu"]
    uwb_cfg = cfg["uwb"]
    init_cfg = cfg["initialization"]
    ekf_cfg = cfg["ekf"]
    campaign = cfg["campaign"]

    fleet = generate_four_node_fleet(
        dt=float(sim["dt_s"]),
        duration=float(sim["duration_s"]),
    )
    if fleet.n_nodes != int(sim["n_nodes"]):
        raise RuntimeError("P2A config and trajectory generator disagree on node count")

    uwb_stride = round(1.0 / (float(uwb_cfg["max_rate_hz"]) * fleet.dt))
    report_stride = round(1.0 / (float(campaign["report_sample_rate_hz"]) * fleet.dt))
    if uwb_stride < 1 or report_stride < 1:
        raise RuntimeError("Invalid sampling-rate configuration")

    n_time = len(fleet.t)
    n_nodes = fleet.n_nodes
    full_mask = _full_mask(n_time, n_nodes, uwb_stride)
    no_mask = np.zeros_like(full_mask)
    p0 = _initial_covariance(init_cfg, n_nodes)

    seed_rows = []
    representative_rows = []
    time_imu = []
    time_full = []
    time_imu_centroid = []
    time_full_centroid = []
    time_imu_shape = []
    time_full_shape = []
    representative_seed = int(campaign["representative_seed"])
    base_seed = int(campaign["seed_offset"])
    start_campaign = time.perf_counter()

    for seed in range(int(campaign["n_seeds"])):
        measured_imu = np.zeros_like(fleet.ideal_imu)
        for node in range(n_nodes):
            rng_imu = np.random.default_rng(base_seed + 10000 * seed + 100 * node + 1)
            imu = simulate_imu_measurements(
                fleet.ideal_imu[:, node],
                fleet.dt,
                rng_imu,
                sigma_accel_mps2=float(imu_cfg["sigma_accel_mps2"]),
                sigma_gyro_rps=float(imu_cfg["sigma_gyro_rps"]),
                sigma_accel_bias_rw_mps2_sqrt_s=float(
                    imu_cfg["sigma_accel_bias_rw_mps2_sqrt_s"]
                ),
                sigma_gyro_bias_rw_rps_sqrt_s=float(
                    imu_cfg["sigma_gyro_bias_rw_rps_sqrt_s"]
                ),
            )
            measured_imu[:, node] = imu.measured

        rng_range = np.random.default_rng(base_seed + 10000 * seed + 9001)
        ranges = _range_data(
            fleet.state,
            float(uwb_cfg["sigma_range_m"]),
            rng_range,
        )

        initial = fleet.state[0].copy()
        rng_init = np.random.default_rng(base_seed + 10000 * seed + 8001)
        initial[:, 2:4] += rng_init.normal(
            0.0,
            float(init_cfg["velocity_error_std_mps"]),
            size=(n_nodes, 2),
        )

        ekf_args = {
            "imu_measurements": measured_imu,
            "pairwise_ranges": ranges,
            "initial_state": initial,
            "initial_covariance": p0,
            "dt": fleet.dt,
            "sigma_range_m": float(uwb_cfg["sigma_range_m"]),
            "sigma_accel_process_mps2": float(ekf_cfg["sigma_accel_process_mps2"]),
            "sigma_gyro_process_rps": float(ekf_cfg["sigma_gyro_process_rps"]),
        }

        start = time.perf_counter()
        imu_only = run_joint_ekf(range_mask=no_mask, **ekf_args)
        imu_runtime = time.perf_counter() - start

        start = time.perf_counter()
        full = run_joint_ekf(range_mask=full_mask, **ekf_args)
        full_runtime = time.perf_counter() - start

        imu_metrics = _fleet_metrics(imu_only.state, fleet.state)
        full_metrics = _fleet_metrics(full.state, fleet.state)
        imu_error = np.linalg.norm(imu_only.state[:, :, :2] - fleet.state[:, :, :2], axis=2)
        full_error = np.linalg.norm(full.state[:, :, :2] - fleet.state[:, :, :2], axis=2)
        time_imu.append(np.sqrt(np.mean(imu_error**2, axis=1)))
        time_full.append(np.sqrt(np.mean(full_error**2, axis=1)))

        imu_error_vector = imu_only.state[:, :, :2] - fleet.state[:, :, :2]
        full_error_vector = full.state[:, :, :2] - fleet.state[:, :, :2]
        time_imu_centroid.append(np.linalg.norm(np.mean(imu_error_vector, axis=1), axis=1))
        time_full_centroid.append(np.linalg.norm(np.mean(full_error_vector, axis=1), axis=1))

        truth_centered = fleet.state[:, :, :2] - np.mean(
            fleet.state[:, :, :2], axis=1, keepdims=True
        )
        imu_centered = imu_only.state[:, :, :2] - np.mean(
            imu_only.state[:, :, :2], axis=1, keepdims=True
        )
        full_centered = full.state[:, :, :2] - np.mean(
            full.state[:, :, :2], axis=1, keepdims=True
        )
        time_imu_shape.append(
            np.sqrt(np.mean(np.linalg.norm(imu_centered - truth_centered, axis=2) ** 2, axis=1))
        )
        time_full_shape.append(
            np.sqrt(np.mean(np.linalg.norm(full_centered - truth_centered, axis=2) ** 2, axis=1))
        )

        row = {
            "seed": seed,
            "imu_fleet_rmse_m": imu_metrics["fleet_position_rmse_m"],
            "full_fleet_rmse_m": full_metrics["fleet_position_rmse_m"],
            "imu_worst_node_rmse_m": imu_metrics["worst_node_position_rmse_m"],
            "full_worst_node_rmse_m": full_metrics["worst_node_position_rmse_m"],
            "imu_centroid_rmse_m": imu_metrics["centroid_position_rmse_m"],
            "full_centroid_rmse_m": full_metrics["centroid_position_rmse_m"],
            "imu_relative_shape_rmse_m": imu_metrics["relative_shape_rmse_m"],
            "full_relative_shape_rmse_m": full_metrics["relative_shape_rmse_m"],
            "imu_yaw_rmse_deg": imu_metrics["fleet_yaw_rmse_deg"],
            "full_yaw_rmse_deg": full_metrics["fleet_yaw_rmse_deg"],
            "imu_runtime_s": imu_runtime,
            "full_runtime_s": full_runtime,
        }
        for node in range(n_nodes):
            row[f"imu_node{node + 1}_rmse_m"] = imu_metrics["node_position_rmse_m"][node]
            row[f"full_node{node + 1}_rmse_m"] = full_metrics["node_position_rmse_m"][node]
        seed_rows.append(row)

        if seed == representative_seed:
            sample_indices = np.arange(0, n_time, report_stride)
            for k in sample_indices:
                for node in range(n_nodes):
                    sl = slice(5 * node, 5 * node + 2)
                    p_pos = full.covariance[k, sl, sl]
                    representative_rows.append(
                        {
                            "time_s": fleet.t[k],
                            "node": node + 1,
                            "truth_x_m": fleet.state[k, node, 0],
                            "truth_y_m": fleet.state[k, node, 1],
                            "imu_x_m": imu_only.state[k, node, 0],
                            "imu_y_m": imu_only.state[k, node, 1],
                            "full_x_m": full.state[k, node, 0],
                            "full_y_m": full.state[k, node, 1],
                            "imu_position_error_m": imu_error[k, node],
                            "full_position_error_m": full_error[k, node],
                            "full_position_sigma_radial_m": math.sqrt(float(np.trace(p_pos))),
                        }
                    )

    time_imu = np.asarray(time_imu)
    time_full = np.asarray(time_full)
    time_imu_centroid = np.asarray(time_imu_centroid)
    time_full_centroid = np.asarray(time_full_centroid)
    time_imu_shape = np.asarray(time_imu_shape)
    time_full_shape = np.asarray(time_full_shape)
    sample_indices = np.arange(0, n_time, report_stride)
    profile_rows = []
    for k in sample_indices:
        profile_rows.append(
            {
                "time_s": fleet.t[k],
                "imu_mean_fleet_error_m": float(np.mean(time_imu[:, k])),
                "imu_p95_fleet_error_m": float(np.quantile(time_imu[:, k], 0.95)),
                "full_mean_fleet_error_m": float(np.mean(time_full[:, k])),
                "full_p95_fleet_error_m": float(np.quantile(time_full[:, k], 0.95)),
                "imu_mean_centroid_error_m": float(np.mean(time_imu_centroid[:, k])),
                "full_mean_centroid_error_m": float(np.mean(time_full_centroid[:, k])),
                "imu_mean_relative_shape_error_m": float(np.mean(time_imu_shape[:, k])),
                "full_mean_relative_shape_error_m": float(np.mean(time_full_shape[:, k])),
            }
        )

    distances = pairwise_distances(fleet.state[:, :, :2])
    speed = np.linalg.norm(fleet.state[:, :, 2:4], axis=2)
    accel_body = np.linalg.norm(fleet.ideal_imu[:, :, :2], axis=2)
    yaw_rate = np.abs(fleet.ideal_imu[:, :, 2])

    imu_rmse = np.array([row["imu_fleet_rmse_m"] for row in seed_rows])
    full_rmse = np.array([row["full_fleet_rmse_m"] for row in seed_rows])
    imu_worst = np.array([row["imu_worst_node_rmse_m"] for row in seed_rows])
    full_worst = np.array([row["full_worst_node_rmse_m"] for row in seed_rows])
    imu_centroid = np.array([row["imu_centroid_rmse_m"] for row in seed_rows])
    full_centroid = np.array([row["full_centroid_rmse_m"] for row in seed_rows])
    imu_shape = np.array([row["imu_relative_shape_rmse_m"] for row in seed_rows])
    full_shape = np.array([row["full_relative_shape_rmse_m"] for row in seed_rows])
    n_opportunities = (n_time - 1) // uwb_stride
    full_communications = n_opportunities * len(all_pairs(n_nodes))

    summary = {
        "design": {
            "n_nodes": n_nodes,
            "n_links": len(all_pairs(n_nodes)),
            "duration_s": float(sim["duration_s"]),
            "imu_rate_hz": float(imu_cfg["rate_hz"]),
            "uwb_rate_hz": float(uwb_cfg["max_rate_hz"]),
            "uwb_sigma_m": float(uwb_cfg["sigma_range_m"]),
            "n_seeds": int(campaign["n_seeds"]),
            "known_initial_position_yaw": True,
            "initial_velocity_error_std_mps": float(init_cfg["velocity_error_std_mps"]),
            "full_communication_count": int(full_communications),
            "full_communication_ratio": 1.0,
        },
        "trajectory": {
            "speed_min_mps": float(np.min(speed)),
            "speed_max_mps": float(np.max(speed)),
            "body_accel_max_mps2": float(np.max(accel_body)),
            "yaw_rate_max_deg_s": float(np.rad2deg(np.max(yaw_rate))),
            "pair_distance_min_m": float(
                min(np.min(values) for values in distances.values())
            ),
            "pair_distance_max_m": float(
                max(np.max(values) for values in distances.values())
            ),
        },
        "imu_only": {
            "fleet_rmse_mean_m": float(np.mean(imu_rmse)),
            "fleet_rmse_std_m": float(np.std(imu_rmse, ddof=1)),
            "fleet_rmse_p95_m": float(np.quantile(imu_rmse, 0.95)),
            "worst_node_rmse_mean_m": float(np.mean(imu_worst)),
            "centroid_rmse_mean_m": float(np.mean(imu_centroid)),
            "relative_shape_rmse_mean_m": float(np.mean(imu_shape)),
        },
        "full_communication": {
            "fleet_rmse_mean_m": float(np.mean(full_rmse)),
            "fleet_rmse_std_m": float(np.std(full_rmse, ddof=1)),
            "fleet_rmse_p95_m": float(np.quantile(full_rmse, 0.95)),
            "worst_node_rmse_mean_m": float(np.mean(full_worst)),
            "centroid_rmse_mean_m": float(np.mean(full_centroid)),
            "relative_shape_rmse_mean_m": float(np.mean(full_shape)),
            "better_than_imu_fraction": float(np.mean(full_rmse < imu_rmse)),
        },
        "improvement": {
            "mean_fleet_rmse_reduction_fraction": float(
                1.0 - np.mean(full_rmse) / np.mean(imu_rmse)
            ),
            "mean_worst_node_rmse_reduction_fraction": float(
                1.0 - np.mean(full_worst) / np.mean(imu_worst)
            ),
            "mean_centroid_rmse_reduction_fraction": float(
                1.0 - np.mean(full_centroid) / np.mean(imu_centroid)
            ),
            "mean_relative_shape_rmse_reduction_fraction": float(
                1.0 - np.mean(full_shape) / np.mean(imu_shape)
            ),
        },
        "campaign_runtime_s": float(time.perf_counter() - start_campaign),
    }

    out = ROOT / "results" / "phase2" / "p2a"
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "summary.md").write_text(
        "# P2A four-node benchmark\n\n"
        f"- IMU-only fleet RMSE: {summary['imu_only']['fleet_rmse_mean_m']:.3f} "
        f"+/- {summary['imu_only']['fleet_rmse_std_m']:.3f} m\n"
        f"- Full-communication fleet RMSE: "
        f"{summary['full_communication']['fleet_rmse_mean_m']:.3f} +/- "
        f"{summary['full_communication']['fleet_rmse_std_m']:.3f} m\n"
        f"- RMSE reduction: "
        f"{100 * summary['improvement']['mean_fleet_rmse_reduction_fraction']:.1f}%\n"
        f"- Full UWB exchanges/run: {full_communications}\n"
        f"- Full communication better than IMU-only: "
        f"{100 * summary['full_communication']['better_than_imu_fraction']:.0f}% of seeds\n"
    )

    _write_csv(out / "seed_metrics.csv", list(seed_rows[0].keys()), seed_rows)
    _write_csv(
        out / "representative.csv",
        list(representative_rows[0].keys()),
        representative_rows,
    )
    _write_csv(out / "time_profile.csv", list(profile_rows[0].keys()), profile_rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
