import numpy as np

from mt_ag.aacopf import literal_all_pairs_aco_transition
from mt_ag.paper_pf import generate_paper_trajectory
from mt_ag.scalable_aacopf import (
    bounded_candidate_aco_transition,
    run_bounded_candidate_aacopf,
)
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory


def test_full_candidate_unconstrained_matches_literal_for_unique_weights():
    particles = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.1],
            [2.0, 0.0, 0.2],
            [3.0, 0.0, 0.3],
        ]
    )
    weights = np.array([0.05, 0.15, 0.30, 0.50])
    literal = literal_all_pairs_aco_transition(
        particles,
        weights,
        alpha=0.5,
        beta=0.0,
        c_lambda=2.0,
    )
    bounded = bounded_candidate_aco_transition(
        particles,
        weights,
        candidate_count=len(particles),
        max_move_fraction=1.0,
        destination_capacity_fraction=1.0,
        alpha=0.5,
        beta=0.0,
        c_lambda=2.0,
    )
    assert np.array_equal(bounded[2], literal[2])
    assert np.allclose(bounded[0], literal[0])


def test_move_and_destination_caps_are_enforced():
    rng = np.random.default_rng(3)
    n_particles = 100
    particles = np.column_stack(
        [
            rng.normal(size=n_particles),
            rng.normal(size=n_particles),
            rng.normal(size=n_particles),
        ]
    )
    weights = np.linspace(1.0, 100.0, n_particles)
    moved, _, parent, diag = bounded_candidate_aco_transition(
        particles,
        weights,
        candidate_count=16,
        max_move_fraction=0.25,
        destination_capacity_fraction=0.03,
        c_lambda=0.0,
    )
    assert diag.moved_count <= 25
    assert diag.max_destination_multiplicity <= 3
    assert diag.candidate_score_count <= n_particles * 16
    assert diag.dense_pair_count == n_particles**2
    assert np.sum(np.any(moved != particles, axis=1)) <= 25
    assert parent.shape == (n_particles,)


def test_equal_weights_produce_no_bounded_moves():
    particles = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.1], [2.0, 0.0, 0.2]]
    )
    weights = np.full(3, 1.0 / 3.0)
    moved, _, parent, diag = bounded_candidate_aco_transition(
        particles,
        weights,
        candidate_count=2,
        c_lambda=0.0,
    )
    assert np.allclose(moved, particles)
    assert np.array_equal(parent, np.arange(3))
    assert diag.moved_count == 0


def test_bounded_filter_runs_and_records_guarded_diagnostics():
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
    result = run_bounded_candidate_aacopf(
        trajectory.increments,
        ranges,
        auxiliary,
        initial,
        sigma_uwb_m=0.05,
        candidate_count=8,
        max_move_fraction=0.25,
        destination_capacity_fraction=0.05,
        c_lambda=0.5,
        truth_state=trajectory.state,
    )
    assert result.estimate.shape == trajectory.state.shape
    assert np.all(np.isfinite(result.estimate))
    assert np.all(result.moved_fraction <= 0.25 + 1e-15)
    assert np.all(result.max_destination_multiplicity <= 4)
    assert np.all(
        (result.unique_parent_fraction >= 0.0) & (result.unique_parent_fraction <= 1.0)
    )
