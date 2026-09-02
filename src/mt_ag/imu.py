from dataclasses import dataclass

import numpy as np

from .geometry import rotation_2d, wrap_angle


@dataclass(frozen=True)
class IMUData:
    measured: np.ndarray  # [a_x^b, a_y^b, omega_z]
    bias: np.ndarray      # corresponding additive sensor bias


def mechanize_planar(initial_state, imu_measurements, dt):
    """Integrate a level planar IMU.

    State convention: [x, y, v_x, v_y, psi] in the navigation frame.
    IMU convention: [a_x^b, a_y^b, omega_z] in the body frame.
    Horizontal accelerations are assumed gravity-compensated. A midpoint yaw
    is used to rotate body acceleration into the navigation frame.
    """
    imu_measurements = np.asarray(imu_measurements, dtype=float)
    state = np.zeros((len(imu_measurements) + 1, 5), dtype=float)
    state[0] = np.asarray(initial_state, dtype=float)

    for k, (a_x_b, a_y_b, omega_z) in enumerate(imu_measurements):
        x, y, v_x, v_y, psi = state[k]
        psi_mid = psi + 0.5 * omega_z * dt
        a_nav = rotation_2d(psi_mid) @ np.array([a_x_b, a_y_b])
        state[k + 1, 0] = x + v_x * dt + 0.5 * a_nav[0] * dt**2
        state[k + 1, 1] = y + v_y * dt + 0.5 * a_nav[1] * dt**2
        state[k + 1, 2] = v_x + a_nav[0] * dt
        state[k + 1, 3] = v_y + a_nav[1] * dt
        state[k + 1, 4] = wrap_angle(psi + omega_z * dt)
    return state


def simulate_imu_measurements(
    ideal_imu,
    dt,
    rng,
    sigma_accel_mps2=0.02,
    sigma_gyro_rps=0.001,
    accel_bias_initial_mps2=(0.0, 0.0),
    gyro_bias_initial_rps=0.0,
    sigma_accel_bias_rw_mps2_sqrt_s=5e-5,
    sigma_gyro_bias_rw_rps_sqrt_s=5e-6,
):
    """Add white sensor noise and slowly varying bias to ideal planar IMU data."""
    ideal_imu = np.asarray(ideal_imu, dtype=float)
    measured = np.empty_like(ideal_imu)
    bias = np.zeros_like(ideal_imu)
    accel_bias = np.asarray(accel_bias_initial_mps2, dtype=float).copy()
    gyro_bias = float(gyro_bias_initial_rps)
    sqrt_dt = np.sqrt(dt)

    for k in range(len(ideal_imu)):
        if k > 0:
            accel_bias += rng.normal(0.0, sigma_accel_bias_rw_mps2_sqrt_s * sqrt_dt, size=2)
            gyro_bias += rng.normal(0.0, sigma_gyro_bias_rw_rps_sqrt_s * sqrt_dt)
        bias[k, :2] = accel_bias
        bias[k, 2] = gyro_bias
        measured[k, :2] = ideal_imu[k, :2] + accel_bias + rng.normal(0.0, sigma_accel_mps2, size=2)
        measured[k, 2] = ideal_imu[k, 2] + gyro_bias + rng.normal(0.0, sigma_gyro_rps)
    return IMUData(measured=measured, bias=bias)
