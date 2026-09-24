"""Build P2B report figures from the committed compact campaign outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/phase2/p2b"
OUT = ROOT / "report/figures/phase2"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 8.5,
    "axes.labelsize": 8.5,
    "axes.titlesize": 9,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
})


def _save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / name)
    plt.close(fig)


def communication_curve(summary: dict) -> None:
    rates = summary["rates"]
    x = np.array([entry["uwb_exchanges"] for entry in rates])
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.85))
    for ax, components in zip(axes, [
        [("fleet_rmse_m", "Fleet position", "tab:blue")],
        [("centroid_rmse_m", "Centroid", "tab:orange"),
         ("relative_shape_rmse_m", "Relative shape", "tab:green")],
    ]):
        for key, label, color in components:
            mean = np.array([entry["metrics"][key]["mean"] for entry in rates])
            sd = np.array([entry["metrics"][key]["std"] for entry in rates])
            ax.errorbar(x, mean, yerr=sd, marker="o", markersize=3.5,
                        linewidth=1.2, capsize=2.2, color=color, label=label)
        ax.set_xlabel("Pairwise UWB exchanges per 60 s run")
        ax.set_ylabel("Across-seed RMSE [m]")
        ax.set_xticks(x)
        ax.set_xticklabels([str(int(value)) for value in x], rotation=35, ha="right")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best")
    axes[0].set_title("(a) Absolute fleet error")
    axes[1].set_title("(b) Centroid and fleet shape")
    _save(fig, "p2b_rate_sweep.pdf")


def temporal_profiles(rows: list[dict[str, float]]) -> None:
    selected = [0.0, 0.5, 2.0, 10.0]
    colors = ["#555555", "#9467bd", "#1f77b4", "#2ca02c"]
    labels = ["0 Hz", "0.5 Hz", "2 Hz", "10 Hz"]
    metrics = [
        ("fleet_rms_error_m", "Fleet RMS error [m]"),
        ("centroid_error_m", "Centroid error [m]"),
        ("shape_rms_error_m", "Relative-shape RMS error [m]"),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(6.8, 4.8), sharex=True)
    for rate, color, label in zip(selected, colors, labels):
        subset = [row for row in rows if row["rate_hz"] == rate]
        t = np.array([row["time_s"] for row in subset])
        for ax, (metric, ylabel) in zip(axes, metrics):
            if rate == 0.0 and metric == "shape_rms_error_m":
                continue  # Use a close scale to resolve the positive-rate shape traces.
            mean = np.array([row[f"mean_{metric}"] for row in subset])
            std = np.array([row[f"std_{metric}"] for row in subset])
            ax.fill_between(t, np.maximum(mean - std, 0), mean + std,
                            color=color, alpha=0.10, linewidth=0)
            ax.plot(t, mean, linewidth=1.15, color=color, label=label)
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.25)
    axes[0].legend(loc="upper left", ncol=4)
    axes[-1].set_ylim(0, 0.32)
    axes[-1].set_title("Positive UWB rates only; 0 Hz omitted to show detail", fontsize=8)
    axes[-1].set_xlabel("Time [s]")
    _save(fig, "p2b_time_profiles.pdf")


def main() -> None:
    summary = json.loads((DATA / "summary.json").read_text())
    with (DATA / "time_profile.csv").open(newline="") as handle:
        rows = [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]
    communication_curve(summary)
    temporal_profiles(rows)


if __name__ == "__main__":
    main()
