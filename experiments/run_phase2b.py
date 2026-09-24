"""P2B: paired uniform all-pairs UWB rate sweep on the frozen P2A fleet."""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from experiments.run_phase2a import (
    _fleet_metrics,
    _full_mask,
    _initial_covariance,
    _range_data,
)
from mt_ag.fleet_simulation import all_pairs, generate_four_node_fleet
from mt_ag.imu import simulate_imu_measurements
from mt_ag.joint_ekf import run_joint_ekf

def uniform_schedule(n_time: int, n_nodes: int, dt: float, rate_hz: float) -> np.ndarray:
    """Activate all links at the first positive multiple of the rate interval."""
    if rate_hz == 0:
        return np.zeros((n_time, n_nodes, n_nodes), dtype=bool)
    if rate_hz < 0:
        raise ValueError("UWB rate cannot be negative")
    stride = round(1.0 / (rate_hz * dt))
    if stride < 1 or not np.isclose(stride * rate_hz * dt, 1.0, atol=1e-12):
        raise ValueError("UWB rate must correspond to an integer number of IMU steps")
    return _full_mask(n_time, n_nodes, stride)


def _time_errors(estimate: np.ndarray, truth: np.ndarray) -> dict[str, np.ndarray]:
    error = estimate[:, :, :2] - truth[:, :, :2]
    centroid = np.mean(error, axis=1)
    fleet_sq = np.mean(np.sum(error**2, axis=2), axis=1)
    centroid_sq = np.sum(centroid**2, axis=1)
    shape_sq = np.mean(np.sum((error - centroid[:, None, :]) ** 2, axis=2), axis=1)
    if not np.allclose(fleet_sq, centroid_sq + shape_sq, rtol=1e-10, atol=1e-12):
        raise RuntimeError("Fleet error decomposition failed")
    return {
        "fleet_rms_error_m": np.sqrt(fleet_sq),
        "centroid_error_m": np.sqrt(centroid_sq),
        "shape_rms_error_m": np.sqrt(shape_sq),
    }


def _write_csv(path: Path, rows: list[dict], *, significant_digits: int | None = None) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        if significant_digits is None:
            writer.writerows(rows)
        else:
            writer.writerows({key: format(value, f".{significant_digits}g")
                              if isinstance(value, float) else value
                              for key, value in row.items()} for row in rows)


def _verify_p2a_endpoints(rows: list[dict], rates: list[float]) -> None:
    with (ROOT / "results/phase2/p2a/seed_metrics.csv").open(newline="") as handle:
        p2a = {int(row["seed"]): row for row in csv.DictReader(handle)}
    if len(p2a) != len({row["seed"] for row in rows}):
        raise RuntimeError("P2A and P2B seed sets differ")
    mapping = {
        "fleet_rmse_m": "fleet_rmse_m",
        "worst_node_rmse_m": "worst_node_rmse_m",
        "centroid_rmse_m": "centroid_rmse_m",
        "relative_shape_rmse_m": "relative_shape_rmse_m",
        "yaw_rmse_deg": "yaw_rmse_deg",
    }
    for row in rows:
        rate = row["rate_hz"]
        if rate not in (0.0, max(rates)):
            continue
        prefix = "imu" if rate == 0 else "full"
        reference = p2a[row["seed"]]
        for key, suffix in mapping.items():
            expected = float(reference[f"{prefix}_{suffix}"])
            if not np.isclose(row[key], expected, atol=1e-10, rtol=1e-10):
                raise RuntimeError(f"P2B {rate:g} Hz does not reproduce P2A seed {row['seed']} {key}")
        for node in range(1, 5):
            expected = float(reference[f"{prefix}_node{node}_rmse_m"])
            if not np.isclose(row[f"node{node}_rmse_m"], expected, atol=1e-10, rtol=1e-10):
                raise RuntimeError(f"P2B {rate:g} Hz disagrees with P2A node {node}")


def main() -> None:
    sweep = yaml.safe_load((ROOT / "configs/phase2b.yaml").read_text())
    base = yaml.safe_load((ROOT / sweep["base_config"]).read_text())
    sim, imu_cfg, uwb_cfg = base["simulation"], base["imu"], base["uwb"]
    init_cfg, ekf_cfg, campaign = base["initialization"], base["ekf"], base["campaign"]
    rates = [float(rate) for rate in sweep["uwb_rates_hz"]]
    if sorted(set(rates)) != rates or rates[0] != 0 or rates[-1] != float(uwb_cfg["max_rate_hz"]):
        raise ValueError("Rates must be unique, ascending, and include the P2A endpoints")

    fleet = generate_four_node_fleet(dt=float(sim["dt_s"]), duration=float(sim["duration_s"]))
    if fleet.n_nodes != int(sim["n_nodes"]):
        raise RuntimeError("P2B and P2A disagree on node count")
    n_time, n_nodes = len(fleet.t), fleet.n_nodes
    masks = {rate: uniform_schedule(n_time, n_nodes, fleet.dt, rate) for rate in rates}
    counts = {rate: int(np.count_nonzero(mask) // 2) for rate, mask in masks.items()}
    full_count = counts[rates[-1]]
    expected_count = ((n_time - 1) // round(1.0 / (rates[-1] * fleet.dt))) * len(
        all_pairs(n_nodes)
    )
    if full_count != expected_count or full_count != 3600:
        raise RuntimeError("P2B full communication count disagrees with P2A")
    for rate in rates[1:]:
        if np.any(masks[rate] & ~masks[rates[-1]]):
            raise RuntimeError("A P2B range opportunity is outside the P2A schedule")

    report_stride = round(1.0 / (float(sweep["profile_sample_rate_hz"]) * fleet.dt))
    if report_stride < 1 or not np.isclose(
        report_stride * float(sweep["profile_sample_rate_hz"]) * fleet.dt, 1.0
    ):
        raise ValueError("Invalid profile sample rate")
    sampled = np.arange(0, n_time, report_stride)
    p0 = _initial_covariance(init_cfg, n_nodes)
    n_seeds = int(campaign["n_seeds"])
    base_seed = int(campaign["seed_offset"])
    rows: list[dict] = []
    profiles = {rate: {key: [] for key in (
        "fleet_rms_error_m", "centroid_error_m", "shape_rms_error_m"
    )} for rate in rates}
    start_campaign = time.perf_counter()

    for seed in range(n_seeds):
        measured_imu = np.zeros_like(fleet.ideal_imu)
        for node in range(n_nodes):
            rng = np.random.default_rng(base_seed + 10000 * seed + 100 * node + 1)
            measured_imu[:, node] = simulate_imu_measurements(
                fleet.ideal_imu[:, node],
                fleet.dt,
                rng,
                sigma_accel_mps2=float(imu_cfg["sigma_accel_mps2"]),
                sigma_gyro_rps=float(imu_cfg["sigma_gyro_rps"]),
                sigma_accel_bias_rw_mps2_sqrt_s=float(
                    imu_cfg["sigma_accel_bias_rw_mps2_sqrt_s"]
                ),
                sigma_gyro_bias_rw_rps_sqrt_s=float(
                    imu_cfg["sigma_gyro_bias_rw_rps_sqrt_s"]
                ),
            ).measured

        rng_range = np.random.default_rng(base_seed + 10000 * seed + 9001)
        ranges = _range_data(fleet.state, float(uwb_cfg["sigma_range_m"]), rng_range)
        initial = fleet.state[0].copy()
        rng_init = np.random.default_rng(base_seed + 10000 * seed + 8001)
        initial[:, 2:4] += rng_init.normal(
            0.0, float(init_cfg["velocity_error_std_mps"]), size=(n_nodes, 2)
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

        for rate in rates:
            start = time.perf_counter()
            estimate = run_joint_ekf(range_mask=masks[rate], **ekf_args)
            runtime = time.perf_counter() - start
            if estimate.range_update_count != counts[rate]:
                raise RuntimeError("EKF range-update count disagrees with schedule")
            metric = _fleet_metrics(estimate.state, fleet.state)
            row = {
                "seed": seed,
                "rate_hz": rate,
                "uwb_exchanges": counts[rate],
                "communication_ratio": counts[rate] / full_count,
                "fleet_rmse_m": metric["fleet_position_rmse_m"],
                "worst_node_rmse_m": metric["worst_node_position_rmse_m"],
                "centroid_rmse_m": metric["centroid_position_rmse_m"],
                "relative_shape_rmse_m": metric["relative_shape_rmse_m"],
                "yaw_rmse_deg": metric["fleet_yaw_rmse_deg"],
                "runtime_s": runtime,
            }
            for node, rmse in enumerate(metric["node_position_rmse_m"], start=1):
                row[f"node{node}_rmse_m"] = rmse
            rows.append(row)
            time_error = _time_errors(estimate.state, fleet.state)
            for key in profiles[rate]:
                profiles[rate][key].append(time_error[key][sampled])
        print(f"P2B completed seed {seed + 1}/{n_seeds}", flush=True)

    _verify_p2a_endpoints(rows, rates)
    out = ROOT / "results/phase2/p2b"
    out.mkdir(parents=True, exist_ok=True)
    _write_csv(out / "seed_metrics.csv", rows)

    profile_rows = []
    for rate in rates:
        arrays = {key: np.asarray(samples) for key, samples in profiles[rate].items()}
        for idx, k in enumerate(sampled):
            row = {"rate_hz": rate, "time_s": float(fleet.t[k])}
            for key, array in arrays.items():
                observations = array[:, idx]
                row[f"mean_{key}"] = float(np.mean(observations))
                row[f"std_{key}"] = float(np.std(observations, ddof=1))
                row[f"p95_{key}"] = float(np.quantile(observations, 0.95))
            profile_rows.append(row)
    _write_csv(out / "time_profile.csv", profile_rows, significant_digits=7)

    by_rate = {rate: [row for row in rows if row["rate_hz"] == rate] for rate in rates}
    baseline = {row["seed"]: row for row in by_rate[0.0]}
    summary = {
        "design": {
            "base_config": sweep["base_config"],
            "rates_hz": rates,
            "n_seeds": n_seeds,
            "duration_s": float(sim["duration_s"]),
            "n_nodes": n_nodes,
            "n_links": len(all_pairs(n_nodes)),
            "full_communication_count": full_count,
            "endpoints_reproduce_p2a": True,
        },
        "rates": [],
        "campaign_runtime_s": float(time.perf_counter() - start_campaign),
    }
    metric_keys = (
        "fleet_rmse_m", "worst_node_rmse_m", "centroid_rmse_m",
        "relative_shape_rmse_m", "yaw_rmse_deg"
    )
    for rate in rates:
        group = by_rate[rate]
        baseline_fleet = np.array([baseline[row["seed"]]["fleet_rmse_m"] for row in group])
        fleet_rmse = np.array([row["fleet_rmse_m"] for row in group])
        item = {
            "rate_hz": rate,
            "uwb_exchanges": counts[rate],
            "communication_ratio": counts[rate] / full_count,
            "paired_fleet_improvement_fraction": float(np.mean(fleet_rmse < baseline_fleet)),
            "mean_fleet_reduction_vs_imu_fraction": float(
                1.0 - np.mean(fleet_rmse) / np.mean(baseline_fleet)
            ),
            "metrics": {},
        }
        for key in metric_keys:
            values = np.array([row[key] for row in group])
            item["metrics"][key] = {
                "mean": float(np.mean(values)),
                "std": float(np.std(values, ddof=1)),
                "p95": float(np.quantile(values, 0.95)),
            }
        summary["rates"].append(item)

    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = ["# P2B uniform all-pairs UWB rate sweep", "", "20 paired P2A seeds; endpoints verified.", ""]
    lines += ["| Rate [Hz] | Exchanges | Fleet RMSE [m] | Centroid [m] | Shape [m] |",
              "|---:|---:|---:|---:|---:|"]
    for item in summary["rates"]:
        metric = item["metrics"]
        lines.append(
            f"| {item['rate_hz']:g} | {item['uwb_exchanges']} | "
            f"{metric['fleet_rmse_m']['mean']:.3f} +/- {metric['fleet_rmse_m']['std']:.3f} | "
            f"{metric['centroid_rmse_m']['mean']:.3f} | "
            f"{metric['relative_shape_rmse_m']['mean']:.3f} |"
        )
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
