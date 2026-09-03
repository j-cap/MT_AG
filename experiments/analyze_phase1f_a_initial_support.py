from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml

from mt_ag.geometry import wrap_angle
from mt_ag.paper_pf import (
    generate_paper_trajectory,
    initialize_random_annulus_particles,
    initialize_structured_annulus_particles,
)
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory

ROOT = Path(__file__).resolve().parents[1]


def summarize(values):
    x = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1f_a.yaml").read_text())
    trajectory = generate_paper_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    sigma_uwb = cfg["uwb"]["sigma_range_m"]
    delta_d = cfg["initialization"]["delta_d_sigma_factor"] * sigma_uwb
    seed_cfg = cfg["randomness"]
    pos_threshold = cfg["convergence"]["position_threshold_m"]
    yaw_threshold = np.deg2rad(cfg["convergence"]["yaw_threshold_deg"])

    def noisy_z0(seed):
        rng = np.random.default_rng(seed_cfg["range_seed_offset"] + seed)
        z = generate_uwb_ranges(trajectory.state[:, :2], auxiliary, rng, sigma_uwb)
        return float(z[0])

    def evaluate(particles):
        pos_error = np.linalg.norm(particles[:, :2] - trajectory.state[0, :2], axis=1)
        yaw_error = np.abs(wrap_angle(particles[:, 2] - trajectory.state[0, 2]))
        correct = (pos_error < pos_threshold) & (yaw_error < yaw_threshold)
        normalized_score = np.sqrt(
            (pos_error / pos_threshold) ** 2 + (yaw_error / yaw_threshold) ** 2
        )
        yaw_compatible = yaw_error < yaw_threshold
        pos_compatible = pos_error < pos_threshold
        return {
            "n_particles": len(particles),
            "correct_mode_count": int(np.sum(correct)),
            "correct_mode_fraction": float(np.mean(correct)),
            "best_normalized_joint_score": float(np.min(normalized_score)),
            "min_position_error_m": float(np.min(pos_error)),
            "min_yaw_error_deg": float(np.rad2deg(np.min(yaw_error))),
            "min_position_error_given_yaw_compatible_m": (
                float(np.min(pos_error[yaw_compatible])) if np.any(yaw_compatible) else None
            ),
            "min_yaw_error_given_position_compatible_deg": (
                float(np.rad2deg(np.min(yaw_error[pos_compatible])))
                if np.any(pos_compatible)
                else None
            ),
        }

    def aggregate(runs):
        return {
            "n": len(runs),
            "fraction_with_at_least_one_correct_mode_particle": float(
                np.mean([r["correct_mode_count"] > 0 for r in runs])
            ),
            "correct_mode_count": summarize([r["correct_mode_count"] for r in runs]),
            "correct_mode_fraction": summarize([r["correct_mode_fraction"] for r in runs]),
            "best_normalized_joint_score": summarize(
                [r["best_normalized_joint_score"] for r in runs]
            ),
        }

    output = {"random_10000": {}, "structured_100x100": {}, "particle_count": {}}

    core_cfg = cfg["core_random"]
    core_runs = []
    for seed in range(core_cfg["n_seeds"]):
        rng = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        particles = initialize_random_annulus_particles(
            noisy_z0(seed), auxiliary[0], core_cfg["n_particles"], rng, delta_d
        )
        run = {"seed": seed, **evaluate(particles)}
        core_runs.append(run)
    output["random_10000"] = {"aggregate": aggregate(core_runs), "runs": core_runs}

    structured_cfg = cfg["structured_control"]
    structured_runs = []
    for seed in range(structured_cfg["n_seeds"]):
        rng = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        particles = initialize_structured_annulus_particles(
            noisy_z0(seed),
            auxiliary[0],
            structured_cfg["n_bearing"],
            structured_cfg["n_yaw"],
            rng,
            delta_d,
        )
        run = {"seed": seed, **evaluate(particles)}
        structured_runs.append(run)
    output["structured_100x100"] = {
        "aggregate": aggregate(structured_runs),
        "runs": structured_runs,
    }

    count_cfg = cfg["particle_count"]
    for n_particles in count_cfg["values"]:
        runs = []
        for seed in range(count_cfg["n_seeds"]):
            rng = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
            particles = initialize_random_annulus_particles(
                noisy_z0(seed), auxiliary[0], n_particles, rng, delta_d
            )
            runs.append({"seed": seed, **evaluate(particles)})
        output["particle_count"][str(n_particles)] = {
            "aggregate": aggregate(runs),
            "runs": runs,
        }

    out_dir = ROOT / "results/phase1/p1f_a"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "initial_support.json").write_text(json.dumps(output, indent=2))
    print(
        json.dumps(
            {
                "random_10000": output["random_10000"]["aggregate"],
                "structured_100x100": output["structured_100x100"]["aggregate"],
                "particle_count": {
                    key: value["aggregate"]
                    for key, value in output["particle_count"].items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
