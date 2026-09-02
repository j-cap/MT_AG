import numpy as np

from .geometry import wrap_angle


def position_error(estimate, truth):
    return np.linalg.norm(np.asarray(estimate)[:, :2] - np.asarray(truth)[:, :2], axis=1)


def summary_metrics(estimate, truth, heading_index=4):
    estimate = np.asarray(estimate)
    truth = np.asarray(truth)
    e = position_error(estimate, truth)
    psi_e = wrap_angle(estimate[:, heading_index] - truth[:, heading_index])
    return {
        "position_rmse_m": float(np.sqrt(np.mean(e**2))),
        "position_mae_m": float(np.mean(e)),
        "position_p95_m": float(np.quantile(e, 0.95)),
        "position_final_m": float(e[-1]),
        "heading_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(psi_e**2)))),
    }
