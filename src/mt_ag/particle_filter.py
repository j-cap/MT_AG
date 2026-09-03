from dataclasses import dataclass

import numpy as np

from .geometry import wrap_angle


@dataclass
class PFResult:
    estimate: np.ndarray  # [x, y, v_x, v_y, psi]
    neff: np.ndarray
    particles_final: np.ndarray
    position_spread: np.ndarray
    yaw_resultant: np.ndarray
    snapshots: dict | None = None


def systematic_resample(weights, rng):
    n = len(weights)
    positions = (rng.random() + np.arange(n)) / n
    cumulative = np.cumsum(weights)
    cumulative[-1] = 1.0
    return np.searchsorted(cumulative, positions)


def circular_mean(theta, weights):
    return np.arctan2(np.sum(weights * np.sin(theta)), np.sum(weights * np.cos(theta)))


def _run_imu_bootstrap_pf_core(
    imu_measurements,
    ranges_meas,
    auxiliary_positions,
    particles,
    dt,
    rng,
    *,
    sigma_process_accel_mps2,
    sigma_process_gyro_rps,
    sigma_uwb_m,
    resample_fraction,
    use_initial_measurement,
    snapshot_steps=None,
):
    imu_measurements = np.asarray(imu_measurements, dtype=float)
    z = np.asarray(ranges_meas, dtype=float)
    aux = np.asarray(auxiliary_positions, dtype=float)
    particles = np.asarray(particles, dtype=float).copy()
    if particles.ndim != 2 or particles.shape[1] != 5:
        raise ValueError("particles must have shape (N, 5)")

    n_particles = len(particles)
    n_steps = len(imu_measurements) + 1
    if len(z) != n_steps or len(aux) != n_steps:
        raise ValueError("range and auxiliary arrays must have len(imu)+1 samples")

    requested_snapshots = set() if snapshot_steps is None else set(snapshot_steps)
    if any(k < 0 or k >= n_steps for k in requested_snapshots):
        raise ValueError("snapshot_steps must be valid state indices")

    particles[:, 4] = wrap_angle(particles[:, 4])
    weights = np.full(n_particles, 1.0 / n_particles)
    estimate = np.zeros((n_steps, 5), dtype=float)
    neff = np.zeros(n_steps, dtype=float)
    position_spread = np.zeros(n_steps, dtype=float)
    yaw_resultant = np.zeros(n_steps, dtype=float)
    snapshots = {}

    def record_estimate(k):
        estimate[k, :4] = np.sum(particles[:, :4] * weights[:, None], axis=0)
        estimate[k, 4] = circular_mean(particles[:, 4], weights)
        delta = particles[:, :2] - estimate[k, :2]
        position_spread[k] = np.sqrt(np.sum(weights * np.sum(delta**2, axis=1)))
        c = np.sum(weights * np.cos(particles[:, 4]))
        s = np.sum(weights * np.sin(particles[:, 4]))
        yaw_resultant[k] = np.hypot(c, s)
        neff[k] = 1.0 / np.sum(weights**2)
        if k in requested_snapshots:
            snapshots[k] = {
                "particles": particles.copy(),
                "weights": weights.copy(),
            }

    def measurement_update(k):
        nonlocal particles, weights
        predicted_range = np.linalg.norm(particles[:, :2] - aux[k], axis=1)
        residual = z[k] - predicted_range
        loglik = -0.5 * (residual / sigma_uwb_m) ** 2
        logw = np.log(weights + 1e-300) + loglik
        logw -= np.max(logw)
        weights = np.exp(logw)
        weights /= np.sum(weights)
        record_estimate(k)
        if neff[k] < resample_fraction * n_particles:
            indices = systematic_resample(weights, rng)
            particles = particles[indices]
            weights = np.full(n_particles, 1.0 / n_particles)

    if use_initial_measurement:
        measurement_update(0)
    else:
        record_estimate(0)

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

    return PFResult(
        estimate=estimate,
        neff=neff,
        particles_final=particles.copy(),
        position_spread=position_spread,
        yaw_resultant=yaw_resultant,
        snapshots=snapshots,
    )


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
    snapshot_steps=None,
):
    """Bootstrap PF with a local Gaussian initial state distribution."""
    particles = rng.normal(
        np.asarray(initial_mean, dtype=float),
        np.asarray(initial_std, dtype=float),
        size=(n_particles, 5),
    )
    return _run_imu_bootstrap_pf_core(
        imu_measurements,
        ranges_meas,
        auxiliary_positions,
        particles,
        dt,
        rng,
        sigma_process_accel_mps2=sigma_process_accel_mps2,
        sigma_process_gyro_rps=sigma_process_gyro_rps,
        sigma_uwb_m=sigma_uwb_m,
        resample_fraction=resample_fraction,
        use_initial_measurement=True,
        snapshot_steps=snapshot_steps,
    )


def run_imu_bootstrap_pf_from_particles(
    imu_measurements,
    ranges_meas,
    auxiliary_positions,
    initial_particles,
    dt,
    rng,
    *,
    sigma_process_accel_mps2=0.025,
    sigma_process_gyro_rps=0.0015,
    sigma_uwb_m=0.12,
    resample_fraction=0.5,
    initial_particles_conditioned_on_z0=False,
    snapshot_steps=None,
):
    """Bootstrap PF starting from an explicitly constructed particle cloud.

    Set ``initial_particles_conditioned_on_z0=True`` when the particles were
    created from the first range measurement. The filter then begins applying
    UWB likelihoods at sample 1 to avoid counting z_0 twice.
    """
    return _run_imu_bootstrap_pf_core(
        imu_measurements,
        ranges_meas,
        auxiliary_positions,
        initial_particles,
        dt,
        rng,
        sigma_process_accel_mps2=sigma_process_accel_mps2,
        sigma_process_gyro_rps=sigma_process_gyro_rps,
        sigma_uwb_m=sigma_uwb_m,
        resample_fraction=resample_fraction,
        use_initial_measurement=not initial_particles_conditioned_on_z0,
        snapshot_steps=snapshot_steps,
    )
