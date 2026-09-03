from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np

from .geometry import wrap_angle
from .particle_filter import circular_mean


@dataclass(frozen=True)
class ACOTransitionDiagnostics:
    moved_count: int
    moved_fraction: float
    unique_parent_count: int
    unique_parent_fraction: float
    unique_destination_count: int
    mean_destination_multiplicity: float
    max_destination_multiplicity: int
    mean_candidate_count: float
    mean_max_probability: float
    median_max_probability: float
    max_probability: float
    candidate_score_count: int
    dense_pair_count: int
    runtime_s: float


@dataclass
class LiteralAACOPFResult:
    estimate: np.ndarray
    map_state: np.ndarray
    map_initial_yaw: np.ndarray
    neff_pre_aco: np.ndarray
    position_spread: np.ndarray
    yaw_resultant: np.ndarray
    correct_mode_mass_pre_aco: np.ndarray | None
    correct_mode_fraction_post_aco: np.ndarray | None
    moved_fraction: np.ndarray
    unique_parent_fraction: np.ndarray
    unique_destination_count: np.ndarray
    max_destination_multiplicity: np.ndarray
    mean_destination_multiplicity: np.ndarray
    mean_candidate_count: np.ndarray
    mean_max_probability: np.ndarray
    max_probability: np.ndarray
    transition_runtime_s: np.ndarray
    particles_final: np.ndarray
    lineage_initial_yaw_final: np.ndarray


def literal_all_pairs_aco_transition(
    particles,
    weights,
    lineage_initial_yaw=None,
    *,
    alpha=1.0,
    beta=1.0,
    c_lambda=0.5,
    epsilon_distance=1.0e-12,
    epsilon_weight=1.0e-15,
):
    """Apply the audited P1F literal-small all-pairs ACO transition.

    Each source particle considers every strictly higher-weight particle as a
    destination. The transition score is

        (w_j - w_i + epsilon_weight)**alpha
        * (1 / (||p_i - p_j|| + epsilon_distance))**beta.

    Scores are normalized over the candidate set of source particle ``i``.
    The best destination is copied if its probability is strictly larger than
    ``lambda_i = c_lambda / K_i``. All destinations are chosen from an
    immutable pre-transition cloud, so the update is synchronous.

    This is an explicit repository interpretation of Han et al.; alpha, beta
    and the threshold parameterization are not source-reported values.
    """
    start = perf_counter()
    particles = np.asarray(particles, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if particles.ndim != 2 or particles.shape[1] != 3:
        raise ValueError("particles must have shape (N, 3)")
    if weights.shape != (len(particles),):
        raise ValueError("weights must have shape (N,)")
    if len(particles) == 0:
        raise ValueError("at least one particle is required")
    if np.any(~np.isfinite(particles)) or np.any(~np.isfinite(weights)):
        raise ValueError("particles and weights must be finite")
    if np.any(weights < 0.0) or np.sum(weights) <= 0.0:
        raise ValueError("weights must be non-negative and have positive sum")
    if alpha < 0.0 or beta < 0.0 or c_lambda < 0.0:
        raise ValueError("alpha, beta and c_lambda must be non-negative")
    if epsilon_distance <= 0.0 or epsilon_weight <= 0.0:
        raise ValueError("epsilon values must be strictly positive")

    n_particles = len(particles)
    normalized_weights = weights / np.sum(weights)
    if lineage_initial_yaw is None:
        lineage = particles[:, 2].copy()
    else:
        lineage = np.asarray(lineage_initial_yaw, dtype=float)
        if lineage.shape != (n_particles,):
            raise ValueError("lineage_initial_yaw must have shape (N,)")
        lineage = lineage.copy()

    positions = particles[:, :2]
    delta = positions[:, None, :] - positions[None, :, :]
    distances = np.linalg.norm(delta, axis=2)

    weight_difference = normalized_weights[None, :] - normalized_weights[:, None]
    candidate_mask = weight_difference > 0.0
    np.fill_diagonal(candidate_mask, False)

    candidate_counts = np.sum(candidate_mask, axis=1)
    candidate_score_count = int(np.sum(candidate_counts))

    scores = np.zeros((n_particles, n_particles), dtype=float)
    if candidate_score_count > 0:
        weight_factor = np.power(
            np.maximum(weight_difference + epsilon_weight, epsilon_weight), alpha
        )
        distance_factor = np.power(1.0 / (distances + epsilon_distance), beta)
        scores[candidate_mask] = (weight_factor * distance_factor)[candidate_mask]

    score_sums = np.sum(scores, axis=1)
    probabilities = np.divide(
        scores,
        score_sums[:, None],
        out=np.zeros_like(scores),
        where=score_sums[:, None] > 0.0,
    )

    best_destination = np.argmax(probabilities, axis=1)
    max_probability = probabilities[np.arange(n_particles), best_destination]
    thresholds = np.full(n_particles, np.inf, dtype=float)
    has_candidates = candidate_counts > 0
    thresholds[has_candidates] = c_lambda / candidate_counts[has_candidates]
    move_mask = has_candidates & (max_probability > thresholds)

    parent_index = np.arange(n_particles)
    parent_index[move_mask] = best_destination[move_mask]
    transitioned_particles = particles[parent_index].copy()
    transitioned_lineage = lineage[parent_index].copy()

    moved_destinations = best_destination[move_mask]
    moved_count = int(np.sum(move_mask))
    if moved_count > 0:
        destination_ids, destination_counts = np.unique(moved_destinations, return_counts=True)
        unique_destination_count = len(destination_ids)
        mean_destination_multiplicity = float(np.mean(destination_counts))
        max_destination_multiplicity = int(np.max(destination_counts))
    else:
        unique_destination_count = 0
        mean_destination_multiplicity = 0.0
        max_destination_multiplicity = 0

    parent_counts = np.bincount(parent_index, minlength=n_particles)
    unique_parent_count = int(np.count_nonzero(parent_counts))
    candidate_max_probabilities = max_probability[has_candidates]
    if len(candidate_max_probabilities) > 0:
        mean_max_probability = float(np.mean(candidate_max_probabilities))
        median_max_probability = float(np.median(candidate_max_probabilities))
        global_max_probability = float(np.max(candidate_max_probabilities))
    else:
        mean_max_probability = 0.0
        median_max_probability = 0.0
        global_max_probability = 0.0

    diagnostics = ACOTransitionDiagnostics(
        moved_count=moved_count,
        moved_fraction=moved_count / n_particles,
        unique_parent_count=unique_parent_count,
        unique_parent_fraction=unique_parent_count / n_particles,
        unique_destination_count=unique_destination_count,
        mean_destination_multiplicity=mean_destination_multiplicity,
        max_destination_multiplicity=max_destination_multiplicity,
        mean_candidate_count=float(np.mean(candidate_counts)),
        mean_max_probability=mean_max_probability,
        median_max_probability=median_max_probability,
        max_probability=global_max_probability,
        candidate_score_count=candidate_score_count,
        dense_pair_count=n_particles * n_particles,
        runtime_s=float(perf_counter() - start),
    )
    return transitioned_particles, transitioned_lineage, parent_index, diagnostics


def run_literal_small_aacopf(
    increments,
    ranges_meas,
    auxiliary_positions,
    initial_particles,
    *,
    sigma_uwb_m=0.12,
    propagation_convention="pre_turn",
    alpha=1.0,
    beta=1.0,
    c_lambda=0.5,
    epsilon_distance=1.0e-12,
    epsilon_weight=1.0e-15,
    initial_particles_conditioned_on_z0=True,
    truth_state=None,
    correct_mode_position_threshold_m=1.0,
    correct_mode_yaw_threshold_rad=np.deg2rad(10.0),
):
    """Run the audited literal-small AACOPF sequence on the paper state.

    The posterior estimate and MAP lineage are recorded before the ACO
    transition. The ACO transition then creates the next generation and the
    weights are reset uniformly, matching the explicit P1E implementation
    contract. No conventional systematic resampling is mixed into this path.
    """
    from .paper_pf import propagate_paper_state

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
    lineage = particles[:, 2].copy()
    weights = np.full(n_particles, 1.0 / n_particles)

    estimate = np.zeros((n_steps, 3), dtype=float)
    map_state = np.zeros((n_steps, 3), dtype=float)
    map_initial_yaw = np.zeros(n_steps, dtype=float)
    neff = np.zeros(n_steps, dtype=float)
    position_spread = np.zeros(n_steps, dtype=float)
    yaw_resultant = np.zeros(n_steps, dtype=float)
    correct_mode_mass = None if truth_state is None else np.zeros(n_steps, dtype=float)
    correct_mode_fraction_post = None if truth_state is None else np.zeros(n_steps, dtype=float)
    moved_fraction = np.zeros(n_steps, dtype=float)
    unique_parent_fraction = np.ones(n_steps, dtype=float)
    unique_destination_count = np.zeros(n_steps, dtype=int)
    max_destination_multiplicity = np.zeros(n_steps, dtype=int)
    mean_destination_multiplicity = np.zeros(n_steps, dtype=float)
    mean_candidate_count = np.zeros(n_steps, dtype=float)
    mean_max_probability = np.zeros(n_steps, dtype=float)
    max_probability = np.zeros(n_steps, dtype=float)
    transition_runtime_s = np.zeros(n_steps, dtype=float)

    def record_posterior(k):
        estimate[k, :2] = np.sum(particles[:, :2] * weights[:, None], axis=0)
        estimate[k, 2] = circular_mean(particles[:, 2], weights)
        map_idx = int(np.argmax(weights))
        map_state[k] = particles[map_idx]
        map_initial_yaw[k] = lineage[map_idx]
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

    def measurement_update(k):
        nonlocal weights
        predicted_range = np.linalg.norm(particles[:, :2] - aux[k], axis=1)
        residual = z[k] - predicted_range
        logw = -0.5 * (residual / sigma_uwb_m) ** 2
        logw -= np.max(logw)
        weights = np.exp(logw)
        weights /= np.sum(weights)
        record_posterior(k)

    def transition(k):
        nonlocal particles, lineage, weights
        particles, lineage, _, diag = literal_all_pairs_aco_transition(
            particles,
            weights,
            lineage,
            alpha=alpha,
            beta=beta,
            c_lambda=c_lambda,
            epsilon_distance=epsilon_distance,
            epsilon_weight=epsilon_weight,
        )
        moved_fraction[k] = diag.moved_fraction
        unique_parent_fraction[k] = diag.unique_parent_fraction
        unique_destination_count[k] = diag.unique_destination_count
        max_destination_multiplicity[k] = diag.max_destination_multiplicity
        mean_destination_multiplicity[k] = diag.mean_destination_multiplicity
        mean_candidate_count[k] = diag.mean_candidate_count
        mean_max_probability[k] = diag.mean_max_probability
        max_probability[k] = diag.max_probability
        transition_runtime_s[k] = diag.runtime_s
        if correct_mode_fraction_post is not None:
            pos_error = np.linalg.norm(particles[:, :2] - truth_state[k, :2], axis=1)
            yaw_error = np.abs(wrap_angle(particles[:, 2] - truth_state[k, 2]))
            correct = (pos_error < correct_mode_position_threshold_m) & (
                yaw_error < correct_mode_yaw_threshold_rad
            )
            correct_mode_fraction_post[k] = float(np.mean(correct))
        weights = np.full(n_particles, 1.0 / n_particles)

    if initial_particles_conditioned_on_z0:
        record_posterior(0)
        if correct_mode_fraction_post is not None:
            pos_error = np.linalg.norm(particles[:, :2] - truth_state[0, :2], axis=1)
            yaw_error = np.abs(wrap_angle(particles[:, 2] - truth_state[0, 2]))
            correct_mode_fraction_post[0] = float(
                np.mean(
                    (pos_error < correct_mode_position_threshold_m)
                    & (yaw_error < correct_mode_yaw_threshold_rad)
                )
            )
    else:
        measurement_update(0)
        transition(0)

    for k, nominal_increment in enumerate(increments):
        particles = propagate_paper_state(
            particles,
            np.broadcast_to(nominal_increment, (n_particles, 2)),
            convention=propagation_convention,
        )
        measurement_update(k + 1)
        transition(k + 1)

    return LiteralAACOPFResult(
        estimate=estimate,
        map_state=map_state,
        map_initial_yaw=map_initial_yaw,
        neff_pre_aco=neff,
        position_spread=position_spread,
        yaw_resultant=yaw_resultant,
        correct_mode_mass_pre_aco=correct_mode_mass,
        correct_mode_fraction_post_aco=correct_mode_fraction_post,
        moved_fraction=moved_fraction,
        unique_parent_fraction=unique_parent_fraction,
        unique_destination_count=unique_destination_count,
        max_destination_multiplicity=max_destination_multiplicity,
        mean_destination_multiplicity=mean_destination_multiplicity,
        mean_candidate_count=mean_candidate_count,
        mean_max_probability=mean_max_probability,
        max_probability=max_probability,
        transition_runtime_s=transition_runtime_s,
        particles_final=particles.copy(),
        lineage_initial_yaw_final=lineage.copy(),
    )
