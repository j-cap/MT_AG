from dataclasses import dataclass
import numpy as np

from .imu import mechanize_planar


@dataclass(frozen=True)
class Trajectory:
    t: np.ndarray
    state: np.ndarray      # [x, y, v_x, v_y, psi]
    ideal_imu: np.ndarray  # [a_x^b, a_y^b, omega_z], length len(t)-1
    dt: float


def generate_curved_trajectory(dt=0.1, duration=60.0):
    """Generate a deterministic planar trajectory with consistent ideal IMU data."""
    t = np.arange(0.0, duration + 0.5 * dt, dt)
    state = np.zeros((len(t), 5), dtype=float)
    state[0] = np.array([0.0, 0.0, 0.75, 0.0, 0.0])
    ideal_imu = np.zeros((len(t) - 1, 3), dtype=float)

    for k in range(len(t) - 1):
        speed = np.linalg.norm(state[k, 2:4])
        a_long = 0.018 * np.sin(0.14 * t[k]) + 0.006 * np.cos(0.05 * t[k])
        omega_z = 0.055 + 0.025 * np.sin(0.07 * t[k])
        a_lat = speed * omega_z
        ideal_imu[k] = np.array([a_long, a_lat, omega_z])
        state[k + 1] = mechanize_planar(state[k], ideal_imu[k:k + 1], dt)[1]
    return Trajectory(t=t, state=state, ideal_imu=ideal_imu, dt=dt)


def auxiliary_trajectory(t, mode="moving"):
    t = np.asarray(t, dtype=float)
    if mode == "stationary":
        return np.column_stack((np.full_like(t, 8.0), np.full_like(t, -4.0)))
    if mode == "moving":
        return np.column_stack((8.0 - 0.10 * t, -4.0 + 0.22 * t + 1.2 * np.sin(0.08 * t)))
    raise ValueError(f"Unknown auxiliary mode: {mode}")
