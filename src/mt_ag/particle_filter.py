from dataclasses import dataclass
import numpy as np
from .geometry import wrap_angle


@dataclass
class PFResult:
    estimate: np.ndarray
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


def run_bootstrap_pf(increments, ranges_meas, auxiliary_positions, initial_mean, initial_std,
                     rng, n_particles=3000, sigma_process_distance=0.006,
                     sigma_process_heading=0.002, sigma_uwb=0.12, resample_fraction=0.5):
    increments = np.asarray(increments, dtype=float)
    z = np.asarray(ranges_meas, dtype=float)
    aux = np.asarray(auxiliary_positions, dtype=float)
    n_steps = len(increments) + 1
    particles = rng.normal(np.asarray(initial_mean), np.asarray(initial_std), size=(n_particles, 3))
    particles[:, 2] = wrap_angle(particles[:, 2])
    weights = np.full(n_particles, 1.0 / n_particles)
    est = np.zeros((n_steps, 3))
    neff = np.zeros(n_steps)

    def update(k):
        nonlocal weights, particles
        pred = np.linalg.norm(particles[:, :2] - aux[k], axis=1)
        residual = z[k] - pred
        loglik = -0.5 * (residual / sigma_uwb) ** 2
        logw = np.log(weights + 1e-300) + loglik
        logw -= np.max(logw)
        weights = np.exp(logw)
        weights /= np.sum(weights)
        neff[k] = 1.0 / np.sum(weights**2)
        est[k, 0:2] = np.sum(particles[:, 0:2] * weights[:, None], axis=0)
        est[k, 2] = circular_mean(particles[:, 2], weights)
        if neff[k] < resample_fraction * n_particles:
            idx = systematic_resample(weights, rng)
            particles = particles[idx]
            weights = np.full(n_particles, 1.0 / n_particles)

    update(0)
    for k, (dL, dpsi) in enumerate(increments):
        dl_i = dL + rng.normal(0.0, sigma_process_distance, n_particles)
        dp_i = dpsi + rng.normal(0.0, sigma_process_heading, n_particles)
        psi = particles[:, 2].copy()
        particles[:, 0] += dl_i * np.cos(psi)
        particles[:, 1] += dl_i * np.sin(psi)
        particles[:, 2] = wrap_angle(psi + dp_i)
        update(k + 1)
    return PFResult(est, neff, particles.copy())
