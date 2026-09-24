from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

from mt_ag.fleet_simulation import all_pairs, generate_four_node_fleet, pairwise_distances

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "report" / "figures" / "phase2"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9,
        "legend.fontsize": 7.3,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
    }
)


def _read_csv(path: Path) -> list[dict[str, float]]:
    with path.open() as handle:
        reader = csv.DictReader(handle)
        return [{key: float(value) for key, value in row.items()} for row in reader]


def _panel(ax, label: str):
    ax.text(
        0.01,
        0.99,
        label,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontweight="bold",
    )


def figure_truth_and_geometry(fleet):
    distances = pairwise_distances(fleet.state[:, :, :2])
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.9))

    ax = axes[0]
    for node in range(fleet.n_nodes):
        ax.plot(
            fleet.state[:, node, 0],
            fleet.state[:, node, 1],
            linewidth=1.4,
            label=f"Node {node + 1}",
        )
        ax.scatter(
            fleet.state[0, node, 0],
            fleet.state[0, node, 1],
            s=22,
            marker="o",
        )
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    _panel(ax, "(a)")

    ax = axes[1]
    for i, j in all_pairs(fleet.n_nodes):
        ax.plot(fleet.t, distances[(i, j)], linewidth=1.1, label=f"{i + 1}-{j + 1}")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("True pairwise distance [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, loc="best")
    _panel(ax, "(b)")

    fig.tight_layout(w_pad=1.3)
    fig.savefig(OUT / "p2a_truth_and_geometry.pdf")
    plt.close(fig)


def figure_imu_truth(fleet):
    speed = np.linalg.norm(fleet.state[:, :, 2:4], axis=2)
    accel = np.linalg.norm(fleet.ideal_imu[:, :, :2], axis=2)
    yaw_rate = np.rad2deg(fleet.ideal_imu[:, :, 2])

    fig, axes = plt.subplots(3, 1, figsize=(6.8, 5.1), sharex=True)
    for node in range(fleet.n_nodes):
        axes[0].plot(fleet.t, speed[:, node], linewidth=1.0, label=f"Node {node + 1}")
        axes[1].plot(fleet.t[:-1], accel[:, node], linewidth=1.0)
        axes[2].plot(fleet.t[:-1], yaw_rate[:, node], linewidth=1.0)

    axes[0].set_ylabel("Speed [m/s]")
    axes[1].set_ylabel("Body accel. norm [m/s$^2$]")
    axes[2].set_ylabel("Yaw rate [deg/s]")
    axes[2].set_xlabel("Time [s]")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    axes[0].legend(ncol=4, loc="upper center")
    _panel(axes[0], "(a)")
    _panel(axes[1], "(b)")
    _panel(axes[2], "(c)")
    fig.tight_layout()
    fig.savefig(OUT / "p2a_truth_imu_signals.pdf")
    plt.close(fig)


def figure_representative_errors(rows):
    by_node = defaultdict(list)
    for row in rows:
        by_node[int(row["node"])].append(row)

    fig, axes = plt.subplots(2, 2, figsize=(6.8, 4.8), sharex=True, sharey=True)
    for node, ax in enumerate(axes.flat, start=1):
        data = by_node[node]
        t = np.array([row["time_s"] for row in data])
        imu = np.array([row["imu_position_error_m"] for row in data])
        full = np.array([row["full_position_error_m"] for row in data])
        ax.plot(t, imu, linestyle="--", linewidth=1.1, label="IMU only")
        ax.plot(t, full, linewidth=1.2, label="All-pair UWB")
        ax.set_title(f"Node {node}")
        ax.grid(True, alpha=0.25)
        _panel(ax, f"({chr(96 + node)})")
    axes[1, 0].set_xlabel("Time [s]")
    axes[1, 1].set_xlabel("Time [s]")
    axes[0, 0].set_ylabel("Position error [m]")
    axes[1, 0].set_ylabel("Position error [m]")
    axes[0, 0].legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "p2a_representative_node_errors.pdf")
    plt.close(fig)


def figure_aggregate(profile, seeds):
    t = np.array([row["time_s"] for row in profile])
    imu_mean = np.array([row["imu_mean_fleet_error_m"] for row in profile])
    imu_p95 = np.array([row["imu_p95_fleet_error_m"] for row in profile])
    full_mean = np.array([row["full_mean_fleet_error_m"] for row in profile])
    full_p95 = np.array([row["full_p95_fleet_error_m"] for row in profile])

    imu_seed = np.array([row["imu_fleet_rmse_m"] for row in seeds])
    full_seed = np.array([row["full_fleet_rmse_m"] for row in seeds])
    maxv = 1.05 * max(np.max(imu_seed), np.max(full_seed))

    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.9))
    ax = axes[0]
    ax.plot(t, imu_mean, linestyle="--", linewidth=1.2, label="IMU mean")
    ax.plot(t, full_mean, linewidth=1.3, label="All-pair mean")
    ax.plot(t, imu_p95, linestyle=":", linewidth=1.0, label="IMU p95")
    ax.plot(t, full_p95, linestyle="-.", linewidth=1.0, label="All-pair p95")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Fleet RMS position error [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    _panel(ax, "(a)")

    ax = axes[1]
    ax.plot([0, maxv], [0, maxv], linestyle="--", linewidth=0.9, label="Equal RMSE")
    ax.scatter(imu_seed, full_seed, s=28)
    ax.set_xlim(0, maxv)
    ax.set_ylim(0, maxv)
    ax.set_xlabel("IMU-only fleet RMSE [m]")
    ax.set_ylabel("All-pair fleet RMSE [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    _panel(ax, "(b)")

    fig.tight_layout(w_pad=1.4)
    fig.savefig(OUT / "p2a_aggregate_performance.pdf")
    plt.close(fig)


def figure_error_decomposition(profile):
    t = np.array([row["time_s"] for row in profile])
    imu_centroid = np.array([row["imu_mean_centroid_error_m"] for row in profile])
    full_centroid = np.array([row["full_mean_centroid_error_m"] for row in profile])
    imu_shape = np.array([row["imu_mean_relative_shape_error_m"] for row in profile])
    full_shape = np.array([row["full_mean_relative_shape_error_m"] for row in profile])

    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.9))
    ax = axes[0]
    ax.plot(t, imu_centroid, linestyle="--", linewidth=1.1, label="IMU only")
    ax.plot(t, full_centroid, linewidth=1.2, label="All-pair UWB")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mean centroid error [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    _panel(ax, "(a)")

    ax = axes[1]
    ax.plot(t, imu_shape, linestyle="--", linewidth=1.1, label="IMU only")
    ax.plot(t, full_shape, linewidth=1.2, label="All-pair UWB")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mean relative-shape RMS error [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    _panel(ax, "(b)")

    fig.tight_layout(w_pad=1.4)
    fig.savefig(OUT / "p2a_error_decomposition.pdf")
    plt.close(fig)


def figure_consistency(rows):
    by_node = defaultdict(list)
    for row in rows:
        by_node[int(row["node"])].append(row)

    fig, axes = plt.subplots(2, 2, figsize=(6.8, 4.8), sharex=True, sharey=True)
    for node, ax in enumerate(axes.flat, start=1):
        data = by_node[node]
        t = np.array([row["time_s"] for row in data])
        error = np.array([row["full_position_error_m"] for row in data])
        sigma = np.array([row["full_position_sigma_radial_m"] for row in data])
        ax.plot(t, error, linewidth=1.2, label="Position error")
        ax.plot(t, 2.0 * sigma, linestyle="--", linewidth=1.0, label="$2\\sqrt{\\mathrm{tr}(P_p)}$")
        ax.set_title(f"Node {node}")
        ax.grid(True, alpha=0.25)
        _panel(ax, f"({chr(96 + node)})")
    axes[1, 0].set_xlabel("Time [s]")
    axes[1, 1].set_xlabel("Time [s]")
    axes[0, 0].set_ylabel("Magnitude [m]")
    axes[1, 0].set_ylabel("Magnitude [m]")
    axes[0, 0].legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "p2a_ekf_consistency.pdf")
    plt.close(fig)


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase2a.yaml").read_text())
    fleet = generate_four_node_fleet(
        dt=float(cfg["simulation"]["dt_s"]),
        duration=float(cfg["simulation"]["duration_s"]),
    )
    result_dir = ROOT / "results" / "phase2" / "p2a"
    representative = _read_csv(result_dir / "representative.csv")
    profile = _read_csv(result_dir / "time_profile.csv")
    seeds = _read_csv(result_dir / "seed_metrics.csv")

    figure_truth_and_geometry(fleet)
    figure_imu_truth(fleet)
    figure_representative_errors(representative)
    figure_aggregate(profile, seeds)
    figure_error_decomposition(profile)
    figure_consistency(representative)

    summary = json.loads((result_dir / "summary.json").read_text())
    print("Generated P2A figures; full-communication fleet RMSE:",
          summary["full_communication"]["fleet_rmse_mean_m"])


if __name__ == "__main__":
    main()
