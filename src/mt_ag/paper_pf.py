from dataclasses import dataclass

import numpy as np

from .geometry import wrap_angle
from .particle_filter import circular_mean, systematic_resample


@dataclass(frozen=True)
class PaperTrajectory:
    t: np.ndarray
    state: np.ndarray  # [x, y, phi]
    increments: np.ndarray  # [delta_L, delta_phi], length len(t)-1
    dt: float


@dataclass
class PaperPFResult:
    estimate: np.ndarray
    map_state: np.ndarray
    map_initial_yaw: np.ndarray
    neff: np.ndarray
    position_spread: np.ndarray
    yaw_resultant: np.ndarray
    unique_fraction_post_transition: np.ndarray
    correct_mode_mass: np.ndarray | None
    particles_final: np.ndarray
    lineage_initial_yaw_final: np.ndarray


def propagate_paper_state(state, increments, convention="pre_turn"):
    """Propagate [x, y, phi] states using Han et al.'s paper-level DR model.

    ``pre_turn`` follows Equation (3): translation uses the previous azimuth.
    ``post_turn`` follows the conflicting Algorithm-1 description: translation
    uses the azimuth after applying the current yaw increment.
    """
    state = np.asarray(state, dtype=float)
    increments = np.asarray(increments, dtype=float)
    if state.shape[-1] != 3:
        raise ValueError("state must end in three components [x, y, phi]")
    if increments.shape[-1] != 2:
        raise ValueError("increments must end in [delta_L, delta_phi]")
    if convention not in {"pre_turn", "post_turn"}:
        raise ValueError("convention must be 'pre_turn' or 'post_turn'")

    phi_old = state[..., 2]
    delta_l = increments[..., 0]
    delta_phi = increments[..., 1]
    heading_for_translation = phi_old if convention == "pre_turn" else phi_old + delta_phi

    out = np.array(state, dtype=float, copy=True)
    out[..., 0] += delta_l * np.cos(heading_for_translation)
    out[..., 1] += delta_l * np.sin(heading_for_translation)
    out[..., 2] = wrap_angle(phi_old + delta_phi)
    return out


def generate_paper_trajectory(dt=0.1, duration=60.0):
    """Generate a deterministic trajectory exactly consistent with Equation (3)."""
    t = np.arange(0.0, duration + 0.5 * dt, dt)
    increments = np.zeros((len(t) - 1, 2), dtype=float)
    state = np.zeros((len(t), 3), dtype=float)
    state[0] = np.array([0.0, 0.0, 0.0])

    for k in range(len(t) - 1):
        speed = 0.75 + 0.08 * np.sin(0.11 * t[k]) + 0.03 * np.cos(0.05 * t[k])
        yaw_rate = 0.055 + 0.025 * np.sin(0.07 * t[k])
        increments[k] = np.array([speed * dt, yaw_rate * dt])
        state[k + 1] = propagate_paper_state(state[k], increments[k], convention="pre_turn")
    return PaperTrajectory(t=t, state=state, increments=increments, dt=dt)


def initialize_random_annulus_particles(z0, auxiliary0, n_particles, rng, delta_d_m):
    """Random first-range annulus and random yaw used as the P1E primary interpretation."""
    auxiliary0 = np.asarray(auxiliary0, dtype=float)
    if auxiliary0.shape != (2,):
        raise ValueError("auxiliary0 must have shape (2,)")
    theta = rng.uniform(-np.pi, np.pi, n_particles)
    radius = z0 + rng.uniform(-delta_d_m, delta_d_m, n_particles)
    radius = np.maximum(radius, 0.0)
    yaw = rng.uniform(-np.pi, np.pi, n_particles)
    particles = np.empty((n_particles, 3), dtype=float)
    particles[:, 0] = auxiliary0[0] + radius * np.cos(theta)
    particles[:, 1] = auxiliary0[1] + radius * np.sin(theta)
    particles[:, 2] = yaw
    return particles


def initialize_structured_annulus_particles(
    z0,
    auxiliary0,
    n_bearing,
    n_yaw,
    rng,
    delta_d_m,
):
    """Structured angular coverage retained only as a P1F-A control variant."""
    auxiliary0 = np.asarray(auxiliary0, dtype=float)
    bearings = -np.pi + 2.0 * np.pi * np.arange(n_bearing) / n_bearing
    yaws = -np.pi + 2.0 * np.pi * np.arange(n_yaw) / n_yaw
    bearing_grid, yaw_grid = np.meshgrid(bearings, yaws, indexing="ij")
    bearing_flat = bearing_grid.ravel()
    yaw_flat = yaw_grid.ravel()
    n_particles = len(bearing_flat)
    radius = z0 + rng.uniform(-delta_d_m, delta_d_m, n_particles)
    radius = np.maximum(radius, 0.0)
    particles = np.empty((n_particles, 3), dtype=float)
    particles[:, 0] = auxiliary0[0] + radius * np.cos(bearing_flat)
    particles[:, 1] = auxiliary0[1] + radius * np.sin(bearing_flat)
    particles[:, 2] = yaw_flat
    return particles


def run_paper_bootstrap_pf(
    increments,
    ranges_meas,
    auxiliary_positions,
    initial_particles,
    rng,
    *,
    sigma_uwb_m=0.12,
    resample_fraction=0.5,
    propagation_convention="pre_turn",
    sigma_delta_l_m=0.0,
    sigma_delta_phi_rad=0.0,
    initial_particles_conditioned_on_z0=True,
    truth_state=None,
    correct_mode_position_threshold_m=1.0,
    correct_mode_yaw_threshold_rad=np.deg2rad(10.0),
):
    """Conventional bootstrap PF for the audited three-state paper model."""
    increments = np.asarray(increments, dtype=float)
    z = np.asarray(ranges_meas, dtype=float)
    aux = np.asarray(auxiliary_positions, dtype=float)
    particles = np.asarray(initial_particles, dtype=float).copy()
    if particles.ndim != 2 or particles.shape[1] != 3:
        raise ValueError("initial_particles must have shape (N, 3)")
    n_particles = len(particles)
    n_steps = len(increments) + 1
    if len(z) != n_steps or aux.shape != (n_steps, 2):
        raise ValueError("ranges and auxiliary positions must have len(increments)+1 samples")
    if truth_state is not None:
        truth_state = np.asarray(truth_state, dtype=float)
        if truth_state.shape != (n_steps, 3):
            raise ValueError("truth_state must have shape (n_steps, 3)")

    particles[:, 2] = wrap_angle(particles[:, 2])
    lineage_initial_yaw = particles[:, 2].copy()
    weights = np.full(n_particles, 1.0 / n_particles)

    estimate = np.zeros((n_steps, 3), dtype=float)
    map_state = np.zeros((n_steps, 3), dtype=float)
    map_initial_yaw = np.zeros(n_steps, dtype=float)
    neff = np.zeros(n_steps, dtype=float)
    position_spread = np.zeros(n_steps, dtype=float)
    yaw_resultant = np.zeros(n_steps, dtype=float)
    unique_fraction = np.ones(n_steps, dtype=float)
    correct_mode_mass = None if truth_state is None else np.zeros(n_steps, dtype=float)

    def record_posterior(k):
        estimate[k, :2] = np.sum(particles[:, :2] * weights[:, None], axis=0)
        estimate[k, 2] = circular_mean(particles[:, 2], weights)
        map_idx = int(np.argmax(weights))
        map_state[k] = particles[map_idx]
        map_initial_yaw[k] = lineage_initial_yaw[map_idx]
        delta = particles[:, :2] - estimate[k, :2]
        position_spread[k] = np.sqrt(np.sum(weights * np.sum(delta**2, axis=1)))
        c = np.sum(weights * np.cos(particles[:, 2]))
        s = np.sum(weights * np.sin(particles[:, 2]))
        yaw_resultant[k] = np.hypot(c, s)
        neff[k] = 1.0 / np.sum(weights**2)
        if correct_mode_mass is not None:
            pos_error = np.linalg.norm(particles[:, :2] - truth_state[k, :2], axis=1)
            yaw_error = np.abs(wrap_angle(particles[:, 2] - truth_state[k, 2]))
            mask = (pos_error < correct_mode_position_threshold_m) & (
                yaw_error < correct_mode_yaw_threshold_rad
            )
            correct_mode_mass[k] = float(np.sum(weights[mask]))

    def update(k):
        nonlocal particles, lineage_initial_yaw, weights
        predicted_range = np.linalg.norm(particles[:, :2] - aux[k], axis=1)
        residual = z[k] - predicted_range
        logw = np.log(weights + 1e-300) - 0.5 * (residual / sigma_uwb_m) ** 2
        logw -= np.max(logw)
        weights = np.exp(logw)
        weights /= np.sum(weights)
        record_posterior(k)
        if neff[k] < resample_fraction * n_particles:
            indices = systematic_resample(weights, rng)
            particles = particles[indices]
            lineage_initial_yaw = lineage_initial_yaw[indices]
            unique_fraction[k] = len(np.unique(indices)) / n_particles
            weights = np.full(n_particles, 1.0 / n_particles)

    if initial_particles_conditioned_on_z0:
        record_posterior(0)
    else:
        update(0)

    for k, nominal_increment in enumerate(increments):
        if sigma_delta_l_m > 0.0 or sigma_delta_phi_rad > 0.0:
            particle_increments = np.empty((n_particles, 2), dtype=float)
            particle_increments[:, 0] = nominal_increment[0] + rng.normal(
                0.0, sigma_delta_l_m, n_particles
            )
            particle_increments[:, 1] = nominal_increment[1] + rng.normal(
                0.0, sigma_delta_phi_rad, n_particles
            )
        else:
            particle_increments = np.broadcast_to(nominal_increment, (n_particles, 2))
        particles = propagate_paper_state(
            particles,
            particle_increments,
            convention=propagation_convention,
        )
        update(k + 1)

    return PaperPFResult(
        estimate=estimate,
        map_state=map_state,
        map_initial_yaw=map_initial_yaw,
        neff=neff,
        position_spread=position_spread,
        yaw_resultant=yaw_resultant,
        unique_fraction_post_transition=unique_fraction,
        correct_mode_mass=correct_mode_mass,
        particles_final=particles.copy(),
        lineage_initial_yaw_final=lineage_initial_yaw.copy(),
    )
