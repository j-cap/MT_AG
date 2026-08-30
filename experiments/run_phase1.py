from pathlib import Path
import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mt_ag.metrics import position_error, summary_metrics  # noqa: E402
from mt_ag.observability import singular_values  # noqa: E402
from mt_ag.particle_filter import run_bootstrap_pf  # noqa: E402
from mt_ag.sensors import generate_uwb_ranges  # noqa: E402
from mt_ag.simulation import (  # noqa: E402
    auxiliary_trajectory,
    generate_curved_trajectory,
    integrate_dr,
    noisy_dr_increments,
)


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1.yaml").read_text())
    rng = np.random.default_rng(cfg["seed"])
    tr = generate_curved_trajectory(**cfg["simulation"])
    ideal_dr = integrate_dr(tr.state[0], tr.increments)
    dr_cfg = cfg["dead_reckoning"]
    noisy_inc = noisy_dr_increments(
        tr.increments,
        rng,
        dr_cfg["sigma_distance_m"],
        dr_cfg["sigma_heading_rad"],
        dr_cfg["heading_bias_rw_rad"],
    )
    noisy_dr = integrate_dr(tr.state[0], noisy_inc)

    aux_moving = auxiliary_trajectory(tr.t, "moving")
    uwb_sigma = cfg["uwb"]["sigma_range_m"]
    z = generate_uwb_ranges(tr.state[:, :2], aux_moving, rng, uwb_sigma)
    pf_cfg = cfg["particle_filter"]
    pf = run_bootstrap_pf(
        noisy_inc,
        z,
        aux_moving,
        tr.state[0],
        pf_cfg["initial_std"],
        rng,
        pf_cfg["n_particles"],
        pf_cfg["sigma_process_distance_m"],
        pf_cfg["sigma_process_heading_rad"],
        uwb_sigma,
        pf_cfg["resample_fraction"],
    )

    metrics = {
        "ideal_dr": summary_metrics(ideal_dr, tr.state),
        "noisy_dr": summary_metrics(noisy_dr, tr.state),
        "pf_uwb": summary_metrics(pf.estimate, tr.state),
    }

    aux_stationary = auxiliary_trajectory(tr.t, "stationary")
    s_stat = singular_values(tr.state, tr.increments, aux_stationary)
    s_move = singular_values(tr.state, tr.increments, aux_moving)
    metrics["observability"] = {
        "stationary_singular_values": [float(x) for x in s_stat],
        "moving_singular_values": [float(x) for x in s_move],
        "stationary_condition_ratio": float(s_stat[-1] / s_stat[0]),
        "moving_condition_ratio": float(s_move[-1] / s_move[0]),
    }

    out = ROOT / "results/phase1"
    figdir = out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2))

    data = np.column_stack((tr.t, tr.state, noisy_dr, pf.estimate))
    header = "t,x_true,y_true,psi_true,x_dr,y_dr,psi_dr,x_pf,y_pf,psi_pf"
    np.savetxt(out / "trajectory.csv", data, delimiter=",", header=header, comments="")
    np.savetxt(
        out / "trajectory_report.csv",
        data[::10],
        delimiter=",",
        header=header,
        comments="",
    )

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.plot(tr.state[:, 0], tr.state[:, 1], label="Ground truth")
    ax.plot(noisy_dr[:, 0], noisy_dr[:, 1], label="Noisy DR")
    ax.plot(pf.estimate[:, 0], pf.estimate[:, 1], label="PF + UWB")
    ax.plot(aux_moving[:, 0], aux_moving[:, 1], "--", label="Auxiliary node")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.axis("equal")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figdir / "trajectory.pdf")
    fig.savefig(figdir / "trajectory.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(tr.t, position_error(noisy_dr, tr.state), label="Noisy DR")
    ax.plot(tr.t, position_error(pf.estimate, tr.state), label="PF + UWB")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("position error [m]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figdir / "position_error.pdf")
    fig.savefig(figdir / "position_error.png", dpi=180)
    plt.close(fig)

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
