from dataclasses import dataclass
import numpy as np
from .geometry import wrap_angle


@dataclass(frozen=True)
class Trajectory:
    t: np.ndarray
    state: np.ndarray  # [x, y, psi]
    increments: np.ndarray  # [dL, dpsi], length len(t)-1


def generate_curved_trajectory(dt=0.1, duration=60.0):
    """Generate a persistently curved deterministic planar trajectory."""
    t = np.arange(0.0, duration + 0.5 * dt, dt)
    n = len(t)
    state = np.zeros((n, 3), dtype=float)
    inc = np.zeros((n - 1, 2), dtype=float)
    for k in range(n - 1):
        v = 0.75 + 0.12 * np.sin(0.11 * t[k])
        omega = 0.055 + 0.025 * np.sin(0.07 * t[k])
        dL = v * dt
        dpsi = omega * dt
        inc[k] = (dL, dpsi)
        x, y, psi = state[k]
        state[k + 1, 0] = x + dL * np.cos(psi)
        state[k + 1, 1] = y + dL * np.sin(psi)
        state[k + 1, 2] = wrap_angle(psi + dpsi)
    return Trajectory(t=t, state=state, increments=inc)


def auxiliary_trajectory(t, mode="moving"):
    t = np.asarray(t, dtype=float)
    if mode == "stationary":
        return np.column_stack((np.full_like(t, 8.0), np.full_like(t, -4.0)))
    if mode == "moving":
        return np.column_stack((8.0 - 0.10 * t, -4.0 + 0.22 * t + 1.2 * np.sin(0.08 * t)))
    raise ValueError(f"Unknown auxiliary mode: {mode}")


def integrate_dr(initial_state, increments):
    increments = np.asarray(increments, dtype=float)
    state = np.zeros((len(increments) + 1, 3), dtype=float)
    state[0] = np.asarray(initial_state, dtype=float)
    for k, (dL, dpsi) in enumerate(increments):
        x, y, psi = state[k]
        state[k + 1, 0] = x + dL * np.cos(psi)
        state[k + 1, 1] = y + dL * np.sin(psi)
        state[k + 1, 2] = wrap_angle(psi + dpsi)
    return state


def noisy_dr_increments(true_increments, rng, sigma_distance=0.004, sigma_heading=0.0015,
                        heading_bias_rw=4e-5):
    true_increments = np.asarray(true_increments, dtype=float)
    out = true_increments.copy()
    out[:, 0] += rng.normal(0.0, sigma_distance, len(out))
    bias = 0.0
    for k in range(len(out)):
        bias += rng.normal(0.0, heading_bias_rw)
        out[k, 1] += bias + rng.normal(0.0, sigma_heading)
    return out
