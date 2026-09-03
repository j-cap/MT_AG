import numpy as np

from mt_ag.aacopf import literal_all_pairs_aco_transition, run_literal_small_aacopf
from mt_ag.paper_pf import generate_paper_trajectory
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory


def test_equal_weights_produce_no_aco_moves():
    particles = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.1],
            [2.0, 0.0, 0.2],
        ]
    )
    weights = np.full(3, 1.0 / 3.0)
    moved, lineage, parent, diag = literal_all_pairs_aco_transition(
        particles, weights, c_lambda=0.0
    )
    assert np.allclose(moved, particles)
    assert np.allclose(lineage, particles[:, 2])
    assert np.array_equal(parent, np.arange(3))
    assert diag.moved_count == 0
    assert diag.candidate_score_count == 0


def test_transition_prefers_weight_distance_score_and_is_synchronous():
    particles = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.1],
            [4.0, 0.0, 0.2],
        ]
    )
    weights = np.array([0.1, 0.3, 0.6])
    lineage = np.array([-0.4, 0.5, 1.2])
    moved, moved_lineage, parent, diag = literal_all_pairs_aco_transition(
        particles,
        weights,
        lineage,
        alpha=1.0,
        beta=1.0,
        c_lambda=0.5,
    )
    # Particle 0 chooses particle 1: 0.2/1 > 0.5/4.
    # Particle 1 has only particle 2 as a higher-weight destination.
    assert np.array_equal(parent, np.array([1, 2, 2]))
    # Synchronous update: particle 0 receives the pre-transition particle 1,
    # not particle 2 after particle 1 itself moves.
    assert np.allclose(moved[0], particles[1])
    assert np.allclose(moved[1], particles[2])
    assert moved_lineage[0] == lineage[1]
    assert moved_lineage[1] == lineage[2]
    assert diag.moved_count == 2
    assert diag.unique_parent_count == 2
    assert diag.max_destination_multiplicity == 1


def test_normalized_threshold_can_block_all_moves():
    particles = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.1],
            [4.0, 0.0, 0.2],
        ]
    )
    weights = np.array([0.1, 0.3, 0.6])
    moved, _, parent, diag = literal_all_pairs_aco_transition(
        particles,
        weights,
        c_lambda=2.0,
    )
    assert np.allclose(moved, particles)
    assert np.array_equal(parent, np.arange(3))
    assert diag.moved_count == 0


def test_only_strictly_higher_weight_destinations_are_allowed():
    particles = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.1],
            [10.0, 0.0, 0.2],
        ]
    )
    weights = np.array([0.2, 0.2, 0.6])
    _, _, parent, _ = literal_all_pairs_aco_transition(
        particles,
        weights,
        c_lambda=0.0,
    )
    # Equal-weight particles 0 and 1 cannot select one another even though
    # they are much closer than particle 2.
    assert parent[0] == 2
    assert parent[1] == 2
    assert parent[2] == 2


def test_literal_small_filter_runs_and_records_transition_diagnostics():
    trajectory = generate_paper_trajectory(dt=0.1, duration=1.0)
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    ranges = generate_uwb_ranges(
        trajectory.state[:, :2], auxiliary, np.random.default_rng(10), 0.05
    )
    rng = np.random.default_rng(11)
    n_particles = 80
    initial = np.column_stack(
        [
            rng.normal(trajectory.state[0, 0], 0.4, n_particles),
            rng.normal(trajectory.state[0, 1], 0.4, n_particles),
            rng.normal(trajectory.state[0, 2], 0.2, n_particles),
        ]
    )
    result = run_literal_small_aacopf(
        trajectory.increments,
        ranges,
        auxiliary,
        initial,
        sigma_uwb_m=0.05,
        alpha=1.0,
        beta=1.0,
        c_lambda=0.5,
        truth_state=trajectory.state,
    )
    assert result.estimate.shape == trajectory.state.shape
    assert np.all(np.isfinite(result.estimate))
    assert np.all((result.moved_fraction >= 0.0) & (result.moved_fraction <= 1.0))
    assert np.all(
        (result.unique_parent_fraction >= 0.0) & (result.unique_parent_fraction <= 1.0)
    )
    assert np.all(result.transition_runtime_s >= 0.0)
    assert np.any(result.moved_fraction[1:] > 0.0)
    assert np.all(
        (result.correct_mode_mass_pre_aco >= 0.0)
        & (result.correct_mode_mass_pre_aco <= 1.0)
    )
