from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mt_ag.imu import mechanize_planar, simulate_imu_measurements  # noqa: E402
from mt_ag.metrics import summary_metrics  # noqa: E402
from mt_ag.particle_filter import run_imu_bootstrap_pf  # noqa: E402
from mt_ag.sensors import generate_uwb_ranges  # noqa: E402
from mt_ag.simulation import auxiliary_trajectory, generate_curved_trajectory  # noqa: E402


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


def aggregate(runs):
    return {
        "n": len(runs),
        "dr_rmse_m": distribution([r["dr"]["position_rmse_m"] for r in runs]),
        "pf_rmse_m": distribution([r["pf"]["position_rmse_m"] for r in runs]),
        "pf_p95_m": distribution([r["pf"]["position_p95_m"] for r in runs]),
        "pf_final_m": distribution([r["pf"]["position_final_m"] for r in runs]),
        "dr_heading_rmse_deg": distribution(
            [r["dr"]["heading_rmse_deg"] for r in runs]
        ),
        "pf_heading_rmse_deg": distribution(
            [r["pf"]["heading_rmse_deg"] for r in runs]
        ),
        "rmse_ratio": distribution([r["rmse_ratio_pf_over_dr"] for r in runs]),
        "runtime_s": distribution([r["runtime_s"] for r in runs]),
        "pf_better_fraction": float(np.mean([r["pf_better"] for r in runs])),
        "mean_neff_fraction": distribution([r["mean_neff_fraction"] for r in runs]),
    }


def compact_aggregate(a):
    return {
        "n": a["n"],
        "dr_rmse_mean_m": a["dr_rmse_m"]["mean"],
        "dr_rmse_std_m": a["dr_rmse_m"]["std"],
        "pf_rmse_mean_m": a["pf_rmse_m"]["mean"],
        "pf_rmse_std_m": a["pf_rmse_m"]["std"],
        "pf_rmse_median_m": a["pf_rmse_m"]["median"],
        "pf_rmse_q95_m": a["pf_rmse_m"]["q95"],
        "pf_rmse_max_m": a["pf_rmse_m"]["max"],
        "dr_heading_rmse_mean_deg": a["dr_heading_rmse_deg"]["mean"],
        "pf_heading_rmse_mean_deg": a["pf_heading_rmse_deg"]["mean"],
        "rmse_ratio_mean": a["rmse_ratio"]["mean"],
        "pf_better_fraction": a["pf_better_fraction"],
        "runtime_mean_s": a["runtime_s"]["mean"],
    }


def count(cfg, section, quick):
    if quick:
        return min(3, int(cfg[section]["n_seeds"]))
    return int(cfg[section]["n_seeds"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run three seeds per condition for a fast implementation check.",
    )
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1b.yaml").read_text())
    trajectory = generate_curved_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    uwb_sigma = float(cfg["uwb"]["sigma_range_m"])
    base_imu = dict(cfg["baseline_imu"])
    base_pf = dict(cfg["baseline_pf"])
    base_initial_std = np.asarray(base_pf["initial_std"], dtype=float)
    seed_cfg = cfg["randomness"]

    def sensor_data(seed, imu_overrides=None):
        imu_cfg = dict(base_imu)
        if imu_overrides:
            imu_cfg.update(imu_overrides)
        rng_sensor = np.random.default_rng(seed_cfg["sensor_seed_offset"] + seed)
        imu = simulate_imu_measurements(
            trajectory.ideal_imu,
            trajectory.dt,
            rng_sensor,
            **imu_cfg,
        )
        ranges = generate_uwb_ranges(
            trajectory.state[:, :2], auxiliary, rng_sensor, uwb_sigma
        )
        return imu, ranges

    def run_case(
        seed,
        *,
        imu_overrides=None,
        initial_mean=None,
        initial_std=None,
        n_particles=None,
        process_accel=None,
        process_gyro=None,
    ):
        imu, ranges = sensor_data(seed, imu_overrides)
        if initial_mean is None:
            initial_mean = trajectory.state[0]
        if initial_std is None:
            initial_std = base_initial_std
        n_particles = int(n_particles or base_pf["n_particles"])

        noisy_dr = mechanize_planar(initial_mean, imu.measured, trajectory.dt)
        rng_pf = np.random.default_rng(seed_cfg["pf_seed_offset"] + seed)
        start = time.perf_counter()
        pf = run_imu_bootstrap_pf(
            imu.measured,
            ranges,
            auxiliary,
            initial_mean,
            initial_std,
            trajectory.dt,
            rng_pf,
            n_particles=n_particles,
            sigma_process_accel_mps2=(
                base_pf["sigma_process_accel_mps2"]
                if process_accel is None
                else process_accel
            ),
            sigma_process_gyro_rps=(
                base_pf["sigma_process_gyro_rps"]
                if process_gyro is None
                else process_gyro
            ),
            sigma_uwb_m=uwb_sigma,
            resample_fraction=base_pf["resample_fraction"],
        )
        runtime = time.perf_counter() - start
        dr_metrics = summary_metrics(noisy_dr, trajectory.state)
        pf_metrics = summary_metrics(pf.estimate, trajectory.state)
        dr_rmse = dr_metrics["position_rmse_m"]
        pf_rmse = pf_metrics["position_rmse_m"]
        return {
            "seed": seed,
            "dr": dr_metrics,
            "pf": pf_metrics,
            "runtime_s": runtime,
            "rmse_ratio_pf_over_dr": float(pf_rmse / dr_rmse),
            "pf_better": bool(pf_rmse < dr_rmse),
            "mean_neff_fraction": float(np.mean(pf.neff) / n_particles),
        }

    raw = {
        "design": {
            "trajectory": "P1A deterministic 60 s curved trajectory",
            "uwb": "moving auxiliary; Gaussian sigma 0.12 m; synchronous 10 Hz",
            "quick_mode": bool(args.quick),
        },
        "multi_seed": {},
        "initial_state_uncertainty": {},
        "imu_quality": {},
        "initial_bias": {},
        "particle_count": {},
        "process_noise": {},
    }

    n_seed = count(cfg, "multi_seed", args.quick)
    baseline_runs = [run_case(seed) for seed in range(n_seed)]
    raw["multi_seed"] = {"aggregate": aggregate(baseline_runs), "runs": baseline_runs}

    n_seed = count(cfg, "initial_state_uncertainty", args.quick)
    for name, values in cfg["initial_state_uncertainty"]["levels"].items():
        error_std = np.asarray(values, dtype=float)
        runs = []
        for seed in range(n_seed):
            rng_initial = np.random.default_rng(
                seed_cfg["initial_state_seed_offset"] + seed
            )
            initial_error = rng_initial.normal(0.0, error_std)
            initial_mean = trajectory.state[0] + initial_error
            runs.append(
                run_case(
                    seed,
                    initial_mean=initial_mean,
                    initial_std=error_std,
                )
            )
        raw["initial_state_uncertainty"][name] = {
            "error_std": [float(x) for x in error_std],
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    n_seed = count(cfg, "imu_quality", args.quick)
    for scale in cfg["imu_quality"]["scales"]:
        scale = float(scale)
        overrides = {
            "sigma_accel_mps2": base_imu["sigma_accel_mps2"] * scale,
            "sigma_gyro_rps": base_imu["sigma_gyro_rps"] * scale,
            "sigma_accel_bias_rw_mps2_sqrt_s": (
                base_imu["sigma_accel_bias_rw_mps2_sqrt_s"] * scale
            ),
            "sigma_gyro_bias_rw_rps_sqrt_s": (
                base_imu["sigma_gyro_bias_rw_rps_sqrt_s"] * scale
            ),
        }
        runs = [run_case(seed, imu_overrides=overrides) for seed in range(n_seed)]
        raw["imu_quality"][str(scale)] = {
            "scale": scale,
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    n_seed = count(cfg, "initial_bias", args.quick)
    for name, values in cfg["initial_bias"]["levels"].items():
        overrides = {
            "accel_bias_initial_mps2": values["accel_bias_mps2"],
            "gyro_bias_initial_rps": values["gyro_bias_rps"],
        }
        runs = [run_case(seed, imu_overrides=overrides) for seed in range(n_seed)]
        raw["initial_bias"][name] = {
            **values,
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    n_seed = count(cfg, "particle_count", args.quick)
    for n_particles in cfg["particle_count"]["values"]:
        n_particles = int(n_particles)
        runs = [run_case(seed, n_particles=n_particles) for seed in range(n_seed)]
        raw["particle_count"][str(n_particles)] = {
            "n_particles": n_particles,
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    n_seed = count(cfg, "process_noise", args.quick)
    for scale in cfg["process_noise"]["scales"]:
        scale = float(scale)
        runs = [
            run_case(
                seed,
                process_accel=base_pf["sigma_process_accel_mps2"] * scale,
                process_gyro=base_pf["sigma_process_gyro_rps"] * scale,
            )
            for seed in range(n_seed)
        ]
        raw["process_noise"][str(scale)] = {
            "scale": scale,
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    compact = {
        "multi_seed": compact_aggregate(raw["multi_seed"]["aggregate"]),
        "initial_state_uncertainty": {},
        "imu_quality": {},
        "initial_bias": {},
        "particle_count": {},
        "process_noise": {},
    }
    for group in compact:
        if group == "multi_seed":
            continue
        for name, result in raw[group].items():
            metadata = {
                key: value
                for key, value in result.items()
                if key not in {"aggregate", "runs"}
            }
            compact[group][name] = {
                **metadata,
                "metrics": compact_aggregate(result["aggregate"]),
            }

    out = ROOT / "results/phase1/p1b"
    out.mkdir(parents=True, exist_ok=True)
    (out / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out / "summary.json").write_text(json.dumps(compact, indent=2))
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
