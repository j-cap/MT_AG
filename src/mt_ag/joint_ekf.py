from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fleet_simulation import all_pairs
from .geometry import wrap_angle


@dataclass(frozen=True)
class JointEKFResult:
    state: np.ndarray  # [time,node,5]
    covariance: np.ndarray  # [time,5N,5N]
    range_update_count: int


def _node_predict(
    state: np.ndarray,
    covariance: np.ndarray,
    imu: np.ndarray,
    dt: float,
    sigma_accel_process_mps2: float,
    sigma_gyro_process_rps: float,
) -> tuple[np.ndarray, np.ndarray]:
    ax_b, ay_b, omega = np.asarray(imu, dtype=float)
    x, y, vx, vy, psi = np.asarray(state, dtype=float)

    psi_mid = psi + 0.5 * omega * dt
    c = np.cos(psi_mid)
    s = np.sin(psi_mid)
    a_nav = np.array([c * ax_b - s * ay_b, s * ax_b + c * ay_b])
    da_dpsi = np.array([-s * ax_b - c * ay_b, c * ax_b - s * ay_b])

    predicted = np.array(
        [
            x + vx * dt + 0.5 * a_nav[0] * dt**2,
            y + vy * dt + 0.5 * a_nav[1] * dt**2,
            vx + a_nav[0] * dt,
            vy + a_nav[1] * dt,
            wrap_angle(psi + omega * dt),
        ],
        dtype=float,
    )

    f = np.eye(5)
    f[0, 2] = dt
    f[1, 3] = dt
    f[0, 4] = 0.5 * da_dpsi[0] * dt**2
    f[1, 4] = 0.5 * da_dpsi[1] * dt**2
    f[2, 4] = da_dpsi[0] * dt
    f[3, 4] = da_dpsi[1] * dt

    sa2 = sigma_accel_process_mps2**2
    sg2 = sigma_gyro_process_rps**2
    q = np.zeros((5, 5), dtype=float)
    q[0, 0] = 0.25 * dt**4 * sa2
    q[0, 2] = q[2, 0] = 0.5 * dt**3 * sa2
    q[2, 2] = dt**2 * sa2
    q[1, 1] = 0.25 * dt**4 * sa2
    q[1, 3] = q[3, 1] = 0.5 * dt**3 * sa2
    q[3, 3] = dt**2 * sa2
    q[4, 4] = dt**2 * sg2

    predicted_covariance = f @ covariance @ f.T + q
    return predicted, predicted_covariance


def _range_update(
    state_flat: np.ndarray,
    covariance: np.ndarray,
    measurement_m: float,
    pair: tuple[int, int],
    sigma_range_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    n_nodes = len(state_flat) // 5
    i, j = pair
    pi = state_flat[5 * i : 5 * i + 2]
    pj = state_flat[5 * j : 5 * j + 2]
    delta = pi - pj
    distance = float(np.linalg.norm(delta))
    if distance < 1e-9:
        return state_flat, covariance

    unit = delta / distance
    h = np.zeros(5 * n_nodes, dtype=float)
    h[5 * i : 5 * i + 2] = unit
    h[5 * j : 5 * j + 2] = -unit

    innovation = float(measurement_m - distance)
    s = float(h @ covariance @ h + sigma_range_m**2)
    gain = covariance @ h / s
    updated = state_flat + gain * innovation
    for node in range(n_nodes):
        updated[5 * node + 4] = wrap_angle(updated[5 * node + 4])

    identity = np.eye(5 * n_nodes)
    kh = np.outer(gain, h)
    updated_covariance = (
        (identity - kh) @ covariance @ (identity - kh).T
        + np.outer(gain, gain) * sigma_range_m**2
    )
    updated_covariance = 0.5 * (updated_covariance + updated_covariance.T)
    return updated, updated_covariance


def run_joint_ekf(
    imu_measurements: np.ndarray,
    pairwise_ranges: np.ndarray,
    range_mask: np.ndarray,
    initial_state: np.ndarray,
    initial_covariance: np.ndarray,
    dt: float,
    sigma_range_m: float,
    sigma_accel_process_mps2: float,
    sigma_gyro_process_rps: float,
) -> JointEKFResult:
    """Run a centralized joint EKF for symmetric mobile IMU+UWB nodes."""
    imu_measurements = np.asarray(imu_measurements, dtype=float)
    pairwise_ranges = np.asarray(pairwise_ranges, dtype=float)
    range_mask = np.asarray(range_mask, dtype=bool)
    initial_state = np.asarray(initial_state, dtype=float)

    n_steps, n_nodes, _ = imu_measurements.shape
    if pairwise_ranges.shape != (n_steps + 1, n_nodes, n_nodes):
        raise ValueError("pairwise_ranges must have shape [time,node,node]")
    if range_mask.shape != pairwise_ranges.shape:
        raise ValueError("range_mask must match pairwise_ranges")

    dim = 5 * n_nodes
    covariance = np.asarray(initial_covariance, dtype=float).copy()
    if covariance.shape != (dim, dim):
        raise ValueError("initial_covariance has incompatible shape")

    history = np.zeros((n_steps + 1, n_nodes, 5), dtype=float)
    covariance_history = np.zeros((n_steps + 1, dim, dim), dtype=float)
    history[0] = initial_state
    covariance_history[0] = covariance
    current = initial_state.reshape(-1).copy()
    update_count = 0

    for k in range(n_steps):
        predicted = np.zeros_like(current)
        f_global = np.zeros((dim, dim), dtype=float)
        q_global = np.zeros((dim, dim), dtype=float)

        for node in range(n_nodes):
            sl = slice(5 * node, 5 * node + 5)
            _, p_node = _node_predict(
                current[sl],
                covariance[sl, sl],
                imu_measurements[k, node],
                dt,
                sigma_accel_process_mps2,
                sigma_gyro_process_rps,
            )

            # Recompute the node transition matrix through covariance propagation identity.
            # Cross-node covariance needs an explicit F; use a numerical-free analytic build here.
            ax_b, ay_b, omega = imu_measurements[k, node]
            psi = current[5 * node + 4]
            psi_mid = psi + 0.5 * omega * dt
            c = np.cos(psi_mid)
            s = np.sin(psi_mid)
            da = np.array([-s * ax_b - c * ay_b, c * ax_b - s * ay_b])
            f = np.eye(5)
            f[0, 2] = dt
            f[1, 3] = dt
            f[0, 4] = 0.5 * da[0] * dt**2
            f[1, 4] = 0.5 * da[1] * dt**2
            f[2, 4] = da[0] * dt
            f[3, 4] = da[1] * dt
            f_global[sl, sl] = f

            predicted_node, _ = _node_predict(
                current[sl],
                covariance[sl, sl],
                imu_measurements[k, node],
                dt,
                sigma_accel_process_mps2,
                sigma_gyro_process_rps,
            )
            predicted[sl] = predicted_node
            q_global[sl, sl] = p_node - f @ covariance[sl, sl] @ f.T

        covariance = f_global @ covariance @ f_global.T + q_global
        current = predicted

        time_index = k + 1
        for pair in all_pairs(n_nodes):
            i, j = pair
            if range_mask[time_index, i, j]:
                current, covariance = _range_update(
                    current,
                    covariance,
                    pairwise_ranges[time_index, i, j],
                    pair,
                    sigma_range_m,
                )
                update_count += 1

        history[time_index] = current.reshape(n_nodes, 5)
        covariance_history[time_index] = covariance

    return JointEKFResult(
        state=history,
        covariance=covariance_history,
        range_update_count=update_count,
    )
