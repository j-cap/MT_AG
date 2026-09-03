from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml

from mt_ag.aacopf import run_literal_small_aacopf
from mt_ag.geometry import wrap_angle
from mt_ag.paper_pf import generate_paper_trajectory
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory

ROOT = Path(__file__).resolve().parents[1]


def controlled_particles(trajectory, auxiliary, cfg, rng):
    n_particles = cfg["n_particles"]
    n_correct = round(cfg["correct_fraction"] * n_particles)
    n_wrong = n_particles - n_correct
    position_std = cfg["position_std_m"]
    yaw_std = np.deg2rad(cfg["yaw_std_deg"])
    truth_position = trajectory.state[0, :2]
    truth_yaw = trajectory.state[0, 2]
    vector = truth_position - auxiliary[0]
    range0 = np.linalg.norm(vector)
    truth_bearing = np.arctan2(vector[1], vector[0])
    wrong_bearing = truth_bearing + np.deg2rad(cfg["wrong_bearing_offset_deg"])
    wrong_center = auxiliary[0] + range0 * np.array(
        [np.cos(wrong_bearing), np.sin(wrong_bearing)]
    )
    wrong_yaw = wrap_angle(truth_yaw + np.deg2rad(cfg["wrong_yaw_offset_deg"]))
    correct = np.column_stack(
        [
            rng.normal(truth_position[0], position_std, n_correct),
            rng.normal(truth_position[1], position_std, n_correct),
            wrap_angle(rng.normal(truth_yaw, yaw_std, n_correct)),
        ]
    )
    wrong = np.column_stack(
        [
            rng.normal(wrong_center[0], position_std, n_wrong),
            rng.normal(wrong_center[1], position_std, n_wrong),
            wrap_angle(rng.normal(wrong_yaw, yaw_std, n_wrong)),
        ]
    )
    particles = np.vstack([correct, wrong])
    rng.shuffle(particles, axis=0)
    return particles


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1f_b.yaml").read_text())
    trajectory = generate_paper_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    case_cfg = cfg["controlled_bimodal"]
    aco_cfg = cfg["aco"]
    seed_cfg = cfg["randomness"]
    yaw_threshold = np.deg2rad(case_cfg["correct_mode_yaw_threshold_deg"])

    output = {}
    for seed in [0, 4]:
        range_rng = np.random.default_rng(seed_cfg["range_seed_offset"] + seed)
        ranges = generate_uwb_ranges(
            trajectory.state[:, :2],
            auxiliary,
            range_rng,
            cfg["uwb"]["sigma_range_m"],
        )
        init_rng = np.random.default_rng(seed_cfg["initializer_seed_offset"] + seed)
        particles = controlled_particles(trajectory, auxiliary, case_cfg, init_rng)
        result = run_literal_small_aacopf(
            trajectory.increments,
            ranges,
            auxiliary,
            particles,
            sigma_uwb_m=cfg["uwb"]["sigma_range_m"],
            alpha=aco_cfg["alpha"],
            beta=aco_cfg["beta"],
            c_lambda=aco_cfg["c_lambda_mechanism_reference"],
            epsilon_distance=aco_cfg["epsilon_distance"],
            epsilon_weight=aco_cfg["epsilon_weight"],
            truth_state=trajectory.state,
            correct_mode_position_threshold_m=case_cfg[
                "correct_mode_position_threshold_m"
            ],
            correct_mode_yaw_threshold_rad=yaw_threshold,
        )
        position_error = np.linalg.norm(
            result.estimate[:, :2] - trajectory.state[:, :2], axis=1
        )
        yaw_error = np.abs(wrap_angle(result.estimate[:, 2] - trajectory.state[:, 2]))
        output[str(seed)] = {
            "t_s": trajectory.t.tolist(),
            "position_error_m": position_error.tolist(),
            "yaw_error_deg": np.rad2deg(yaw_error).tolist(),
            "correct_mode_mass_pre_aco": result.correct_mode_mass_pre_aco.tolist(),
            "correct_mode_fraction_post_aco": (
                result.correct_mode_fraction_post_aco.tolist()
            ),
            "moved_fraction": result.moved_fraction.tolist(),
            "unique_parent_fraction": result.unique_parent_fraction.tolist(),
            "max_destination_multiplicity": (
                result.max_destination_multiplicity.tolist()
            ),
            "mean_max_transition_probability": result.mean_max_probability.tolist(),
            "neff_fraction_pre_aco": (
                result.neff_pre_aco / case_cfg["n_particles"]
            ).tolist(),
        }

    out_dir = ROOT / "results/phase1/p1f_b"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "diagnostic_traces.json").write_text(json.dumps(output, indent=2))

    for seed, trace in output.items():
        unique = np.asarray(trace["unique_parent_fraction"])
        moved = np.asarray(trace["moved_fraction"])
        idx = int(np.argmin(unique[1:]) + 1)
        print(
            f"seed={seed}: strongest collapse at t={trace['t_s'][idx]:.2f}s, "
            f"unique_parent_fraction={unique[idx]:.4f}, "
            f"moved_fraction={moved[idx]:.4f}, "
            f"correct_post={trace['correct_mode_fraction_post_aco'][idx]:.4f}"
        )


if __name__ == "__main__":
    main()
