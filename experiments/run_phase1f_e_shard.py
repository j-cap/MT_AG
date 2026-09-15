from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.aacopf import run_literal_small_aacopf
from mt_ag.paper_pf import (
    generate_paper_trajectory,
    initialize_random_annulus_particles,
    run_paper_bootstrap_pf,
)
from mt_ag.simulation import auxiliary_trajectory
from run_phase1f_e import initial_support_count, method_metrics, observability_summary

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", required=True, choices=["stationary", "constant_bearing"])
    parser.add_argument("--seed-start", required=True, type=int)
    parser.add_argument("--n-seeds", required=True, type=int)
    args = parser.parse_args()

    cfg = yaml.safe_load((ROOT / "configs/phase1f_e.yaml").read_text())
    trajectory = generate_paper_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(
        trajectory.t,
        args.geometry,
        target_positions=(
            trajectory.state[:, :2] if args.geometry == "constant_bearing" else None
        ),
    )
    n_particles = cfg["comparison"]["n_particles"]
    sigma_uwb = cfg["uwb"]["sigma_range_m"]
    delta_d = cfg["initialization"]["delta_d_sigma_factor"] * sigma_uwb
    support_yaw = np.deg2rad(cfg["support_diagnostics"]["correct_region_yaw_threshold_deg"])
    conv_yaw = np.deg2rad(cfg["convergence"]["yaw_threshold_deg"])
    seed_cfg = cfg["randomness"]
    aco_cfg = cfg["frozen_aacopf"]
    geometry_index = cfg["geometry_cases"].index(args.geometry)

    runs = []
    for seed in range(args.seed_start, args.seed_start + args.n_seeds):
        noise_rng = np.random.default_rng(seed_cfg["range_seed_offset"] + seed)
        shared_noise = noise_rng.normal(0.0, sigma_uwb, len(trajectory.t))
        ideal_ranges = np.linalg.norm(trajectory.state[:, :2] - auxiliary, axis=1)
        ranges = ideal_ranges + shared_noise

        init_rng = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        initial_particles = initialize_random_annulus_particles(
            ranges[0], auxiliary[0], n_particles, init_rng, delta_d
        )
        support_count = initial_support_count(
            initial_particles,
            trajectory.state[0],
            cfg["support_diagnostics"]["correct_region_position_threshold_m"],
            support_yaw,
        )

        pf_rng = np.random.default_rng(
            seed_cfg["pf_seed_offset"] + 1000 * geometry_index + seed
        )
        start = time.perf_counter()
        pf_result = run_paper_bootstrap_pf(
            trajectory.increments,
            ranges,
            auxiliary,
            initial_particles.copy(),
            pf_rng,
            sigma_uwb_m=sigma_uwb,
            resample_fraction=cfg["conventional_pf"]["resample_fraction"],
            propagation_convention=cfg["paper_model"]["propagation_convention"],
            sigma_delta_l_m=cfg["paper_model"]["sigma_delta_l_m"],
            sigma_delta_phi_rad=cfg["paper_model"]["sigma_delta_phi_rad"],
            initial_particles_conditioned_on_z0=True,
            truth_state=trajectory.state,
            correct_mode_position_threshold_m=cfg["convergence"]["position_threshold_m"],
            correct_mode_yaw_threshold_rad=conv_yaw,
        )
        pf_runtime = time.perf_counter() - start
        pf_metrics = method_metrics(
            pf_result.estimate,
            pf_result.map_state,
            trajectory.state,
            auxiliary,
            cfg["convergence"],
            trajectory.dt,
        )
        pf_metrics.update(
            {
                "runtime_s": float(pf_runtime),
                "final_correct_mode_mass": float(pf_result.correct_mode_mass[-1]),
                "resampling_events": int(
                    np.sum(pf_result.unique_fraction_post_transition < 1.0 - 1e-15)
                ),
                "minimum_unique_fraction": float(
                    np.min(pf_result.unique_fraction_post_transition)
                ),
                "final_position_spread_m": float(pf_result.position_spread[-1]),
                "final_yaw_resultant": float(pf_result.yaw_resultant[-1]),
            }
        )

        start = time.perf_counter()
        aco_result = run_literal_small_aacopf(
            trajectory.increments,
            ranges,
            auxiliary,
            initial_particles.copy(),
            sigma_uwb_m=sigma_uwb,
            propagation_convention=cfg["paper_model"]["propagation_convention"],
            alpha=aco_cfg["alpha"],
            beta=aco_cfg["beta"],
            c_lambda=aco_cfg["c_lambda"],
            epsilon_distance=aco_cfg["epsilon_distance"],
            epsilon_weight=aco_cfg["epsilon_weight"],
            initial_particles_conditioned_on_z0=True,
            truth_state=trajectory.state,
            correct_mode_position_threshold_m=cfg["convergence"]["position_threshold_m"],
            correct_mode_yaw_threshold_rad=conv_yaw,
        )
        aco_runtime = time.perf_counter() - start
        aco_metrics = method_metrics(
            aco_result.estimate,
            aco_result.map_state,
            trajectory.state,
            auxiliary,
            cfg["convergence"],
            trajectory.dt,
        )
        update_slice = slice(1, None)
        min_parent = float(np.min(aco_result.unique_parent_fraction[update_slice]))
        max_multiplicity = int(np.max(aco_result.max_destination_multiplicity[update_slice]))
        aco_metrics.update(
            {
                "runtime_s": float(aco_runtime),
                "transition_runtime_s": float(
                    np.sum(aco_result.transition_runtime_s[update_slice])
                ),
                "final_correct_mode_mass_pre_aco": float(
                    aco_result.correct_mode_mass_pre_aco[-1]
                ),
                "final_correct_mode_fraction_post_aco": float(
                    aco_result.correct_mode_fraction_post_aco[-1]
                ),
                "minimum_unique_parent_fraction": min_parent,
                "catastrophic_collapse": bool(
                    min_parent
                    < cfg["ancestry_diagnostics"]["catastrophic_unique_parent_fraction"]
                ),
                "dominant_clone": bool(
                    max_multiplicity
                    >= cfg["ancestry_diagnostics"]["dominant_clone_fraction"] * n_particles
                ),
                "maximum_destination_multiplicity": max_multiplicity,
                "mean_moved_fraction": float(
                    np.mean(aco_result.moved_fraction[update_slice])
                ),
                "final_position_spread_m": float(aco_result.position_spread[-1]),
                "final_yaw_resultant": float(aco_result.yaw_resultant[-1]),
            }
        )
        runs.append(
            {
                "seed": seed,
                "initial_support_count": support_count,
                "pf": pf_metrics,
                "aacopf": aco_metrics,
            }
        )

    payload = {
        "geometry": args.geometry,
        "seed_start": args.seed_start,
        "n_seeds": args.n_seeds,
        "n_particles": n_particles,
        "observability": observability_summary(trajectory, auxiliary, cfg),
        "runs": runs,
    }
    out_dir = ROOT / "results/phase1/p1f_e/shards"
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / f"{args.geometry}_{args.seed_start:03d}.json"
    output.write_text(json.dumps(payload, indent=2))
    print(json.dumps({"output": str(output), "geometry": args.geometry, "n": len(runs)}, indent=2))


if __name__ == "__main__":
    main()
