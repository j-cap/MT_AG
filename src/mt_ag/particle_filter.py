from dataclasses import dataclass

import numpy as np

from .geometry import wrap_angle


@dataclass
class PFResult:
    estimate: np.ndarray  # [x, y, v_x, v_y, psi]
    neff: np.ndarray
    particles_final: np.ndarray


def systematic_resample(weights, rng):
    n = len(weights)
    positions = (rng.random() + np.arange(n)) / n
    cumulative = np.cumsum(weights)
    cumulative[-1] = 1.0
    return np.searchsorted(cumulative, positions)


def circular_mean(theta, weights):
    return np.arctan2(np.sum(weights * np.sin(theta)), np.sum(weights * np.cos(theta)))


def run_imu_bootstrap_pf(
    imu_measurements,
    ranges_meas,
    auxiliary_positions,
    initial_mean,
    initial_std,
    dt,
    rng,
    n_particles=5000,
    sigma_process_accel_mps2=0.025,
    sigma_process_gyro_rps=0.0015,
    sigma_uwb_m=0.12,
    resample_fraction=0.5,
):
    """Bootstrap PF using planar accelerometer/gyro measurements for propagation."""
    imu_measurements = np.asarray(imu_measurements, dtype=float)
    z = np.asarray(ranges_meas, dtype=float)
    aux = np.asarray(auxiliary_positions, dtype=float)
    n_steps = len(imu_measurements) + 1
    particles = rng.normal(
        np.asarray(initial_mean, dtype=float),
        np.asarray(initial_std, dtype=float),
        size=(n_particles, 5),
    )
    particles[:, 4] = wrap_angle(particles[:, 4])
    weights = np.full(n_particles, 1.0 / n_particles)
    estimate = np.zeros((n_steps, 5), dtype=float)
    neff = np.zeros(n_steps, dtype=float)

    def measurement_update(k):
        nonlocal particles, weights
        predicted_range = np.linalg.norm(particles[:, :2] - aux[k], axis=1)
        residual = z[k] - predicted_range
        loglik = -0.5 * (residual / sigma_uwb_m) ** 2
        logw = np.log(weights + 1e-300) + loglik
        logw -= np.max(logw)
        weights = np.exp(logw)
        weights /= np.sum(weights)
        neff[k] = 1.0 / np.sum(weights**2)
        estimate[k, :4] = np.sum(particles[:, :4] * weights[:, None], axis=0)
        estimate[k, 4] = circular_mean(particles[:, 4], weights)
        if neff[k] < resample_fraction * n_particles:
            indices = systematic_resample(weights, rng)
            particles = particles[indices]
            weights = np.full(n_particles, 1.0 / n_particles)

    measurement_update(0)
    for k, measurement in enumerate(imu_measurements):
        a_x_b = measurement[0] + rng.normal(0.0, sigma_process_accel_mps2, n_particles)
        a_y_b = measurement[1] + rng.normal(0.0, sigma_process_accel_mps2, n_particles)
        omega_z = measurement[2] + rng.normal(0.0, sigma_process_gyro_rps, n_particles)
        psi = particles[:, 4].copy()
        psi_mid = psi + 0.5 * omega_z * dt
        c = np.cos(psi_mid)
        s = np.sin(psi_mid)
        a_x_n = c * a_x_b - s * a_y_b
        a_y_n = s * a_x_b + c * a_y_b
        particles[:, 0] += particles[:, 2] * dt + 0.5 * a_x_n * dt**2
        particles[:, 1] += particles[:, 3] * dt + 0.5 * a_y_n * dt**2
        particles[:, 2] += a_x_n * dt
        particles[:, 3] += a_y_n * dt
        particles[:, 4] = wrap_angle(psi + omega_z * dt)
        measurement_update(k + 1)
    return PFResult(estimate=estimate, neff=neff, particles_final=particles.copy())
