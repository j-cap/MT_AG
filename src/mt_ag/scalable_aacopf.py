from __future__ import annotations

from time import perf_counter

import numpy as np

from .aacopf import ACOTransitionDiagnostics, LiteralAACOPFResult
from .geometry import wrap_angle
from .particle_filter import circular_mean


def bounded_candidate_aco_transition(
    particles,
    weights,
    lineage_initial_yaw=None,
    *,
    candidate_count=16,
    max_move_fraction=0.25,
    destination_capacity_fraction=0.025,
    alpha=0.5,
    beta=0.0,
    c_lambda=2.0,
    epsilon_distance=1.0e-12,
    epsilon_weight=1.0e-15,
):
    """Apply the P1F-G bounded-candidate AACOPF adaptation.

    Each source particle scores only a global pool of the ``candidate_count``
    highest-weight particles and retains strictly higher-weight members of that
    pool as possible destinations. The movement threshold keeps the literal
    normalization ``c_lambda / K_i`` with ``K_i`` equal to the total number of
    strictly higher-weight particles in the full cloud. That count is obtained
    from sorted weights in O(N log N) without constructing a dense pair matrix.
    This avoids the unintended threshold inflation that would occur if K_i were
    replaced by the truncated candidate-pool size.

    Two explicit diversity guards are applied after the score/threshold test:

    1. at most ``max_move_fraction`` of the cloud may move in one transition;
    2. one destination may receive at most ``destination_capacity_fraction`` of
       the cloud as incoming copies.

    Candidate and guard decisions use the immutable pre-transition cloud. The
    candidate scoring cost is O(N K), plus O(N log N) sorting. This is a
    repository adaptation, not a literal claim about Han et al.'s implementation.
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
    if int(candidate_count) != candidate_count or candidate_count < 1:
        raise ValueError("candidate_count must be a positive integer")
    if not 0.0 <= max_move_fraction <= 1.0:
        raise ValueError("max_move_fraction must lie in [0, 1]")
    if not 0.0 < destination_capacity_fraction <= 1.0:
        raise ValueError("destination_capacity_fraction must lie in (0, 1]")
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

    pool_size = min(int(candidate_count), n_particles)
    elite_index = np.argsort(-normalized_weights, kind="stable")[:pool_size]
    elite_weights = normalized_weights[elite_index]
    weight_difference = elite_weights[None, :] - normalized_weights[:, None]
    candidate_mask = weight_difference > 0.0
    candidate_counts = np.sum(candidate_mask, axis=1)
    candidate_score_count = int(np.sum(candidate_counts))

    scores = np.zeros((n_particles, pool_size), dtype=float)
    if candidate_score_count > 0:
        weight_factor = np.power(
            np.maximum(weight_difference + epsilon_weight, epsilon_weight), alpha
        )
        if beta == 0.0:
            score_matrix = weight_factor
        else:
            source_positions = particles[:, None, :2]
            destination_positions = particles[elite_index][None, :, :2]
            distances = np.linalg.norm(source_positions - destination_positions, axis=2)
            distance_factor = np.power(1.0 / (distances + epsilon_distance), beta)
            score_matrix = weight_factor * distance_factor
        scores[candidate_mask] = score_matrix[candidate_mask]

    score_sums = np.sum(scores, axis=1)
    probabilities = np.divide(
        scores,
        score_sums[:, None],
        out=np.zeros_like(scores),
        where=score_sums[:, None] > 0.0,
    )
    best_column = np.argmax(probabilities, axis=1)
    best_destination = elite_index[best_column]
    max_probability = probabilities[np.arange(n_particles), best_column]

    sorted_weights = np.sort(normalized_weights)
    full_candidate_counts = n_particles - np.searchsorted(
        sorted_weights, normalized_weights, side="right"
    )
    thresholds = np.full(n_particles, np.inf, dtype=float)
    has_candidates = (candidate_counts > 0) & (full_candidate_counts > 0)
    thresholds[has_candidates] = c_lambda / full_candidate_counts[has_candidates]
    eligible = has_candidates & (max_probability > thresholds)
    margin = max_probability - thresholds
    margin[~eligible] = -np.inf

    eligible_index = np.flatnonzero(eligible)
    if len(eligible_index) > 0:
        priority = eligible_index[
            np.lexsort((eligible_index, -margin[eligible_index]))
        ]
    else:
        priority = eligible_index

    max_moves = int(np.floor(max_move_fraction * n_particles))
    if max_move_fraction > 0.0 and max_moves == 0:
        max_moves = 1
    priority = priority[:max_moves]

    destination_capacity = max(
        1, int(np.ceil(destination_capacity_fraction * n_particles))
    )
    accepted = []
    incoming_count = {}
    for source_index in priority:
        destination = int(best_destination[source_index])
        count = incoming_count.get(destination, 0)
        if count >= destination_capacity:
            continue
        incoming_count[destination] = count + 1
        accepted.append(int(source_index))

    move_mask = np.zeros(n_particles, dtype=bool)
    if accepted:
        move_mask[np.asarray(accepted, dtype=int)] = True

    parent_index = np.arange(n_particles)
    parent_index[move_mask] = best_destination[move_mask]
    transitioned_particles = particles[parent_index].copy()
    transitioned_lineage = lineage[parent_index].copy()

    moved_destinations = best_destination[move_mask]
    moved_count = int(np.sum(move_mask))
    if moved_count > 0:
        destination_ids, destination_counts = np.unique(
            moved_destinations, return_counts=True
        )
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


def run_bounded_candidate_aacopf(
    increments,
    ranges_meas,
    auxiliary_positions,
    initial_particles,
    *,
    sigma_uwb_m=0.12,
    propagation_convention="pre_turn",
    candidate_count=16,
    max_move_fraction=0.25,
    destination_capacity_fraction=0.025,
    alpha=0.5,
    beta=0.0,
    c_lambda=2.0,
    epsilon_distance=1.0e-12,
    epsilon_weight=1.0e-15,
    initial_particles_conditioned_on_z0=True,
    truth_state=None,
    correct_mode_position_threshold_m=1.0,
    correct_mode_yaw_threshold_rad=0.17453292519943295,
):
    """Run the P1F-G bounded-candidate AACOPF adaptation on the paper state."""
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
        particles, lineage, _, diag = bounded_candidate_aco_transition(
            particles,
            weights,
            lineage,
            candidate_count=candidate_count,
            max_move_fraction=max_move_fraction,
            destination_capacity_fraction=destination_capacity_fraction,
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
