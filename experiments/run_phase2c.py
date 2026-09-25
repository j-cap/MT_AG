"""P2C: matched-budget link selection on the frozen four-node benchmark."""

from __future__ import annotations

import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from experiments.run_phase2a import _fleet_metrics, _initial_covariance, _range_data
from experiments.run_phase2b import _time_errors, _write_csv
from mt_ag.fleet_simulation import all_pairs, generate_four_node_fleet
from mt_ag.imu import simulate_imu_measurements
from mt_ag.joint_ekf import run_joint_ekf
from mt_ag.link_selection import (
    PredictiveLinkSelector,
    baseline_mask,
    eligible_mask,
    exchange_capacity,
)


def _reference_rows(budgets: list[int]) -> dict[tuple[int, int], dict]:
    with (ROOT / "results/phase2/p2b/seed_metrics.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    result = {(int(row["seed"]), int(row["uwb_exchanges"])): row for row in rows
              if int(row["uwb_exchanges"]) in budgets}
    return result


def main() -> None:
    cfg = yaml.safe_load((ROOT / "configs/phase2c.yaml").read_text())
    base = yaml.safe_load((ROOT / cfg["base_config"]).read_text())
    sim, imu_cfg, uwb_cfg = base["simulation"], base["imu"], base["uwb"]
    init_cfg, ekf_cfg, campaign = base["initialization"], base["ekf"], base["campaign"]
    budgets = [int(value) for value in cfg["budgets"]]
    policies = cfg["policies"]
    if budgets != [180, 360, 720] or policies != [
        "cyclic", "random", "geometry", "information"
    ]:
        raise ValueError("P2C protocol requires its frozen three budgets and four policies")

    fleet = generate_four_node_fleet(dt=float(sim["dt_s"]), duration=float(sim["duration_s"]))
    n_time, n_nodes = len(fleet.t), fleet.n_nodes
    if n_nodes != 4 or n_time != 6001:
        raise ValueError("P2C expects the frozen P2A fleet")
    stride = round(1 / (float(uwb_cfg["max_rate_hz"]) * fleet.dt))
    capacities = {b: exchange_capacity(n_time, stride, b, len(all_pairs(n_nodes)))
                  for b in budgets}
    eligibility = {b: eligible_mask(capacities[b], n_nodes) for b in budgets}
    if any(int(np.sum(capacities[b])) != b for b in budgets):
        raise RuntimeError("Invalid per-tick budget")
    report_stride = round(1 / (float(cfg["profile_sample_rate_hz"]) * fleet.dt))
    if not np.isclose(report_stride * float(cfg["profile_sample_rate_hz"]) * fleet.dt, 1):
        raise ValueError("Profile output must land on simulation steps")
    sample_indices = np.arange(0, n_time, report_stride)
    p0 = _initial_covariance(init_cfg, n_nodes)
    n_seeds = int(campaign["n_seeds"])
    base_seed = int(campaign["seed_offset"])
    p2b = _reference_rows(budgets)
    if len(p2b) != n_seeds * len(budgets):
        raise RuntimeError("P2B comparison lacks a seed/budget endpoint")
    rows, link_rows, example = [], [], []
    histories = {(b, policy): [] for b in cfg["profile_budgets"] for policy in policies}
    start_campaign = time.perf_counter()

    for seed in range(n_seeds):
        measured_imu = np.zeros_like(fleet.ideal_imu)
        for node in range(n_nodes):
            rng = np.random.default_rng(base_seed + 10000 * seed + 100 * node + 1)
            measured_imu[:, node] = simulate_imu_measurements(
                fleet.ideal_imu[:, node], fleet.dt, rng,
                sigma_accel_mps2=float(imu_cfg["sigma_accel_mps2"]),
                sigma_gyro_rps=float(imu_cfg["sigma_gyro_rps"]),
                sigma_accel_bias_rw_mps2_sqrt_s=float(
                    imu_cfg["sigma_accel_bias_rw_mps2_sqrt_s"]),
                sigma_gyro_bias_rw_rps_sqrt_s=float(
                    imu_cfg["sigma_gyro_bias_rw_rps_sqrt_s"]),
            ).measured
        ranges = _range_data(
            fleet.state, float(uwb_cfg["sigma_range_m"]),
            np.random.default_rng(base_seed + 10000 * seed + 9001),
        )
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
        for budget in budgets:
            capacity = capacities[budget]
            for policy in policies:
                selector = None
                if policy in ("cyclic", "random"):
                    rng = (np.random.default_rng(base_seed + 10000 * seed + 7301 + budget)
                           if policy == "random" else None)
                    mask = baseline_mask(capacity, n_nodes, policy, rng)
                    if not np.array_equal(np.count_nonzero(mask, axis=(1, 2)) // 2,
                                          capacity):
                        raise RuntimeError("Fixed policy violates per-tick capacity")
                    selections = [(int(t), i, j) for t in np.flatnonzero(capacity)
                                  for i, j in all_pairs(n_nodes) if mask[t, i, j]]
                else:
                    mask = eligibility[budget]
                    selector = PredictiveLinkSelector(
                        policy, capacity, float(uwb_cfg["sigma_range_m"]), fleet.dt,
                        n_nodes, float(cfg["geometry_memory_s"]),
                        float(cfg["geometry_prior"]),
                    )
                start = time.perf_counter()
                estimate = run_joint_ekf(range_mask=mask, range_selector=selector, **ekf_args)
                runtime = time.perf_counter() - start
                if selector is not None:
                    selections = [(t, i, j) for t, i, j, _ in selector.records]
                    selected_at = Counter(t for t, _, _ in selections)
                    if any(selected_at.get(t, 0) != int(count)
                           for t, count in enumerate(capacity) if count):
                        raise RuntimeError("Adaptive policy violates per-tick capacity")
                if estimate.range_update_count != budget or len(selections) != budget:
                    raise RuntimeError("Range count disagrees with P2C budget")

                metrics = _fleet_metrics(estimate.state, fleet.state)
                row = {
                    "seed": seed, "budget": budget, "communication_ratio": budget / 3600,
                    "policy": policy,
                    "fleet_rmse_m": metrics["fleet_position_rmse_m"],
                    "worst_node_rmse_m": metrics["worst_node_position_rmse_m"],
                    "centroid_rmse_m": metrics["centroid_position_rmse_m"],
                    "relative_shape_rmse_m": metrics["relative_shape_rmse_m"],
                    "yaw_rmse_deg": metrics["fleet_yaw_rmse_deg"],
                    "runtime_s": runtime,
                }
                for node, value in enumerate(metrics["node_position_rmse_m"], start=1):
                    row[f"node{node}_rmse_m"] = value
                rows.append(row)
                link_counts = Counter((i, j) for _, i, j in selections)
                for i, j in all_pairs(n_nodes):
                    link_rows.append({"seed": seed, "budget": budget, "policy": policy,
                                      "pair": f"{i+1}-{j+1}", "exchanges": link_counts[i, j]})
                if seed == 0 and budget == int(cfg["representative_budget"]):
                    example.extend({"policy": policy, "time_s": round(t * fleet.dt, 2),
                                    "pair": f"{i+1}-{j+1}"} for t, i, j in selections)
                if budget in cfg["profile_budgets"]:
                    errors = _time_errors(estimate.state, fleet.state)
                    histories[budget, policy].append(
                        {key: values[sample_indices] for key, values in errors.items()}
                    )
        print(f"P2C completed seed {seed + 1}/{n_seeds}", flush=True)

    out = ROOT / "results/phase2/p2c"
    out.mkdir(parents=True, exist_ok=True)
    _write_csv(out / "seed_metrics.csv", rows)
    _write_csv(out / "link_counts.csv", link_rows)
    _write_csv(out / "representative_schedule.csv", example)
    profile_rows = []
    for budget, policy in histories:
        history = histories[budget, policy]
        for idx, time_index in enumerate(sample_indices):
            result = {"budget": budget, "policy": policy,
                      "time_s": float(fleet.t[time_index])}
            for key in history[0]:
                observations = np.array([run[key][idx] for run in history])
                result[f"mean_{key}"] = float(np.mean(observations))
                result[f"std_{key}"] = float(np.std(observations, ddof=1))
            profile_rows.append(result)
    _write_csv(out / "time_profile.csv", profile_rows, significant_digits=7)

    summary = {
        "protocol": {
            "base_config": cfg["base_config"], "n_seeds": n_seeds,
            "budgets": budgets, "policies": policies,
            "max_opportunities": 600, "pairwise_links": 6,
            "geometry_memory_s": cfg["geometry_memory_s"],
            "geometry_prior": cfg["geometry_prior"],
            "per_tick_capacity_shared": True,
            "same_seed_sensor_draws_as_p2a_p2b": True,
        },
        "results": [], "campaign_runtime_s": float(time.perf_counter() - start_campaign),
    }
    keys = ("fleet_rmse_m", "worst_node_rmse_m", "centroid_rmse_m",
            "relative_shape_rmse_m", "yaw_rmse_deg")
    lookup = {(row["seed"], row["budget"], row["policy"]): row for row in rows}
    for budget in budgets:
        for policy in policies:
            group = [lookup[seed, budget, policy] for seed in range(n_seeds)]
            item = {"budget": budget, "policy": policy, "metrics": {}}
            for key in keys:
                vals = np.array([row[key] for row in group])
                item["metrics"][key] = {"mean": float(np.mean(vals)),
                                        "std": float(np.std(vals, ddof=1)),
                                        "p95": float(np.quantile(vals, .95))}
            item["paired_differences_m"] = {}
            for reference in ("cyclic", "random", "p2b_uniform"):
                for key in ("fleet_rmse_m", "worst_node_rmse_m",
                            "relative_shape_rmse_m"):
                    benchmark = [float(p2b[seed, budget][key]) if reference == "p2b_uniform"
                                 else lookup[seed, budget, reference][key]
                                 for seed in range(n_seeds)]
                    diff = np.array([row[key] for row in group]) - benchmark
                    item["paired_differences_m"][f"{key}_vs_{reference}"] = {
                        "mean": float(np.mean(diff)),
                        "std": float(np.std(diff, ddof=1)),
                        "wins": int(np.count_nonzero(diff < 0)),
                    }
            summary["results"].append(item)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = ["# P2C matched-budget link selection", "",
             "All policies have identical exchange counts at every 10 Hz opportunity.", "",
             "| Budget | Policy | Fleet RMSE [m] | Worst node [m] | Shape [m] |",
             "|---:|---|---:|---:|---:|"]
    for item in summary["results"]:
        metrics = item["metrics"]
        lines.append(f"| {item['budget']} | {item['policy']} | "
                     f"{metrics['fleet_rmse_m']['mean']:.3f} +/- "
                     f"{metrics['fleet_rmse_m']['std']:.3f} | "
                     f"{metrics['worst_node_rmse_m']['mean']:.3f} | "
                     f"{metrics['relative_shape_rmse_m']['mean']:.3f} |")
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
