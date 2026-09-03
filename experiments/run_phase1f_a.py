from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.geometry import wrap_angle
from mt_ag.paper_pf import (
    generate_paper_trajectory,
    initialize_random_annulus_particles,
    initialize_structured_annulus_particles,
    propagate_paper_state,
    run_paper_bootstrap_pf,
)
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory

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


def terminal_convergence_time(mask, dt, hold_steps):
    mask = np.asarray(mask, dtype=bool)
    false_indices = np.flatnonzero(~mask)
    start = 0 if len(false_indices) == 0 else int(false_indices[-1] + 1)
    if len(mask) - start < hold_steps:
        return None
    return float(start * dt)


def aggregate(runs):
    pose_times = [r["pose_convergence_s"] for r in runs if r["pose_convergence_s"] is not None]
    return {
        "n": len(runs),
        "n_particles": runs[0]["n_particles"],
        "pose_success_fraction": float(np.mean([r["pose_convergence_s"] is not None for r in runs])),
        "pose_convergence_s": distribution(pose_times) if pose_times else None,
        "position_rmse_m": distribution([r["position_rmse_m"] for r in runs]),
        "late_position_rmse_m": distribution([r["late_position_rmse_m"] for r in runs]),
        "final_position_error_m": distribution([r["final_position_error_m"] for r in runs]),
        "late_yaw_rmse_deg": distribution([r["late_yaw_rmse_deg"] for r in runs]),
        "final_yaw_error_deg": distribution([r["final_yaw_error_deg"] for r in runs]),
        "final_initial_yaw_lineage_error_deg": distribution(
            [r["final_initial_yaw_lineage_error_deg"] for r in runs]
        ),
        "final_correct_mode_mass": distribution([r["final_correct_mode_mass"] for r in runs]),
        "late_correct_mode_mass": distribution([r["late_correct_mode_mass"] for r in runs]),
        "mean_neff_fraction": distribution([r["mean_neff_fraction"] for r in runs]),
        "resampling_events": distribution([r["resampling_events"] for r in runs]),
        "mean_unique_fraction_at_resampling": distribution(
            [r["mean_unique_fraction_at_resampling"] for r in runs]
        ),
        "runtime_s": distribution([r["runtime_s"] for r in runs]),
    }


def compact(value):
    return {
        "n": value["n"],
        "n_particles": value["n_particles"],
        "pose_success_fraction": value["pose_success_fraction"],
        "pose_convergence_median_s": (
            None if value["pose_convergence_s"] is None else value["pose_convergence_s"]["median"]
        ),
        "late_position_rmse_mean_m": value["late_position_rmse_m"]["mean"],
        "late_yaw_rmse_mean_deg": value["late_yaw_rmse_deg"]["mean"],
        "final_initial_yaw_lineage_error_mean_deg": value[
            "final_initial_yaw_lineage_error_deg"
        ]["mean"],
        "final_correct_mode_mass_mean": value["final_correct_mode_mass"]["mean"],
        "mean_neff_fraction": value["mean_neff_fraction"]["mean"],
        "resampling_events_mean": value["resampling_events"]["mean"],
        "runtime_mean_s": value["runtime_s"]["mean"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1f_a.yaml").read_text())
    trajectory = generate_paper_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    uwb_cfg = cfg["uwb"]
    init_cfg = cfg["initialization"]
    pf_cfg = cfg["particle_filter"]
    conv_cfg = cfg["convergence"]
    seed_cfg = cfg["randomness"]
    delta_d = init_cfg["delta_d_sigma_factor"] * uwb_cfg["sigma_range_m"]
    hold_steps = round(conv_cfg["hold_time_s"] / trajectory.dt)
    late_steps = round(conv_cfg["late_window_s"] / trajectory.dt)
    yaw_threshold = np.deg2rad(conv_cfg["yaw_threshold_deg"])

    # Exact model-consistency check for the primary Equation-(3) convention.
    reconstruction = np.zeros_like(trajectory.state)
    reconstruction[0] = trajectory.state[0]
    for k, increment in enumerate(trajectory.increments):
        reconstruction[k + 1] = propagate_paper_state(
            reconstruction[k], increment, convention="pre_turn"
        )
    reconstruction_max_abs_error = float(np.max(np.abs(reconstruction - trajectory.state)))

    range_cache = {}

    def ranges_for_seed(seed):
        if seed not in range_cache:
            rng = np.random.default_rng(seed_cfg["range_seed_offset"] + seed)
            range_cache[seed] = generate_uwb_ranges(
                trajectory.state[:, :2],
                auxiliary,
                rng,
                uwb_cfg["sigma_range_m"],
            )
        return range_cache[seed]

    def run_case(seed, *, convention, init_mode, n_particles=None, n_bearing=None, n_yaw=None):
        ranges = ranges_for_seed(seed)
        rng_init = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        if init_mode == "random":
            particles = initialize_random_annulus_particles(
                ranges[0], auxiliary[0], n_particles, rng_init, delta_d
            )
        elif init_mode == "structured":
            particles = initialize_structured_annulus_particles(
                ranges[0], auxiliary[0], n_bearing, n_yaw, rng_init, delta_d
            )
        else:
            raise ValueError(init_mode)

        rng_pf = np.random.default_rng(seed_cfg["pf_seed_offset"] + seed)
        start = time.perf_counter()
        result = run_paper_bootstrap_pf(
            trajectory.increments,
            ranges,
            auxiliary,
            particles,
            rng_pf,
            sigma_uwb_m=uwb_cfg["sigma_range_m"],
            resample_fraction=pf_cfg["resample_fraction"],
            propagation_convention=convention,
            sigma_delta_l_m=pf_cfg["sigma_delta_l_m"],
            sigma_delta_phi_rad=pf_cfg["sigma_delta_phi_rad"],
            initial_particles_conditioned_on_z0=True,
            truth_state=trajectory.state,
            correct_mode_position_threshold_m=conv_cfg["position_threshold_m"],
            correct_mode_yaw_threshold_rad=yaw_threshold,
        )
        runtime = time.perf_counter() - start

        position_error = np.linalg.norm(result.estimate[:, :2] - trajectory.state[:, :2], axis=1)
        yaw_error = np.abs(wrap_angle(result.estimate[:, 2] - trajectory.state[:, 2]))
        pose_ok = (position_error < conv_cfg["position_threshold_m"]) & (yaw_error < yaw_threshold)
        late = slice(max(0, len(position_error) - late_steps), None)
        resampled = result.unique_fraction_post_transition < 1.0 - 1e-15
        mean_unique = (
            float(np.mean(result.unique_fraction_post_transition[resampled])) if np.any(resampled) else 1.0
        )
        initial_yaw_error = abs(wrap_angle(result.map_initial_yaw[-1] - trajectory.state[0, 2]))
        return {
            "seed": seed,
            "n_particles": len(particles),
            "propagation_convention": convention,
            "initialization": init_mode,
            "pose_convergence_s": terminal_convergence_time(pose_ok, trajectory.dt, hold_steps),
            "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
            "late_position_rmse_m": float(np.sqrt(np.mean(position_error[late] ** 2))),
            "final_position_error_m": float(position_error[-1]),
            "late_yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error[late] ** 2)))),
            "final_yaw_error_deg": float(np.rad2deg(yaw_error[-1])),
            "final_initial_yaw_lineage_error_deg": float(np.rad2deg(initial_yaw_error)),
            "final_correct_mode_mass": float(result.correct_mode_mass[-1]),
            "late_correct_mode_mass": float(np.mean(result.correct_mode_mass[late])),
            "mean_neff_fraction": float(np.mean(result.neff) / len(particles)),
            "resampling_events": int(np.sum(resampled)),
            "mean_unique_fraction_at_resampling": mean_unique,
            "runtime_s": float(runtime),
        }

    def n_seeds(section):
        n = cfg[section]["n_seeds"]
        return min(3, n) if args.quick else n

    raw = {
        "design": {
            "state": "[x,y,phi]",
            "inputs": "[delta_L,delta_phi]",
            "truth_convention": "pre_turn Equation (3)",
            "auxiliary": "moving and globally known",
            "uwb_sigma_m": uwb_cfg["sigma_range_m"],
            "delta_d_m": delta_d,
            "dr_noise": "none in P1F-A primary campaign",
            "first_range_handling": "initialization conditioned on z0; likelihood begins at z1",
            "quick": args.quick,
        },
        "reconstruction_max_abs_error": reconstruction_max_abs_error,
        "core_random": {},
        "propagation_sensitivity": {},
        "structured_control": {},
        "particle_count": {},
    }

    core = cfg["core_random"]
    core_particles = 3000 if args.quick else core["n_particles"]
    runs = [
        run_case(
            seed,
            convention=core["propagation_convention"],
            init_mode="random",
            n_particles=core_particles,
        )
        for seed in range(n_seeds("core_random"))
    ]
    raw["core_random"] = {"aggregate": aggregate(runs), "runs": runs}

    sensitivity = cfg["propagation_sensitivity"]
    sens_particles = 3000 if args.quick else sensitivity["n_particles"]
    for convention in sensitivity["conventions"]:
        runs = [
            run_case(
                seed,
                convention=convention,
                init_mode="random",
                n_particles=sens_particles,
            )
            for seed in range(n_seeds("propagation_sensitivity"))
        ]
        raw["propagation_sensitivity"][convention] = {
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    structured = cfg["structured_control"]
    grid = 40 if args.quick else structured["n_bearing"]
    runs = [
        run_case(
            seed,
            convention=structured["propagation_convention"],
            init_mode="structured",
            n_bearing=grid,
            n_yaw=grid,
        )
        for seed in range(n_seeds("structured_control"))
    ]
    raw["structured_control"] = {"aggregate": aggregate(runs), "runs": runs}

    count_values = [2500, 10000] if args.quick else cfg["particle_count"]["values"]
    for n_particles in count_values:
        runs = [
            run_case(
                seed,
                convention=cfg["particle_count"]["propagation_convention"],
                init_mode="random",
                n_particles=n_particles,
            )
            for seed in range(n_seeds("particle_count"))
        ]
        raw["particle_count"][str(n_particles)] = {
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    summary = {
        "reconstruction_max_abs_error": reconstruction_max_abs_error,
        "core_random": compact(raw["core_random"]["aggregate"]),
        "propagation_sensitivity": {
            key: compact(value["aggregate"])
            for key, value in raw["propagation_sensitivity"].items()
        },
        "structured_control": compact(raw["structured_control"]["aggregate"]),
        "particle_count": {
            key: compact(value["aggregate"]) for key, value in raw["particle_count"].items()
        },
    }

    out = ROOT / "results/phase1/p1f_a"
    out.mkdir(parents=True, exist_ok=True)
    (out / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
