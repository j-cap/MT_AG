import numpy as np

from mt_ag.paper_pf import (
    generate_paper_trajectory,
    initialize_random_annulus_particles,
    propagate_paper_state,
    run_paper_bootstrap_pf,
)
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory


def test_paper_trajectory_reconstructs_under_equation3_convention():
    trajectory = generate_paper_trajectory(dt=0.1, duration=5.0)
    state = trajectory.state[0].copy()
    reconstructed = [state.copy()]
    for increment in trajectory.increments:
        state = propagate_paper_state(state, increment, convention="pre_turn")
        reconstructed.append(state.copy())
    assert np.max(np.abs(np.asarray(reconstructed) - trajectory.state)) < 1e-12


def test_post_turn_convention_is_distinct_for_nonzero_yaw_increment():
    state = np.array([0.0, 0.0, 0.2])
    increment = np.array([1.0, 0.3])
    pre = propagate_paper_state(state, increment, convention="pre_turn")
    post = propagate_paper_state(state, increment, convention="post_turn")
    assert not np.allclose(pre[:2], post[:2])
    assert np.isclose(pre[2], post[2])


def test_random_annulus_initialization_respects_radial_bound():
    rng = np.random.default_rng(4)
    auxiliary = np.array([2.0, -1.0])
    z0 = 5.0
    delta_d = 0.4
    particles = initialize_random_annulus_particles(z0, auxiliary, 2000, rng, delta_d)
    radii = np.linalg.norm(particles[:, :2] - auxiliary, axis=1)
    assert particles.shape == (2000, 3)
    assert np.all(radii >= z0 - delta_d)
    assert np.all(radii <= z0 + delta_d)
    assert np.all(particles[:, 2] >= -np.pi)
    assert np.all(particles[:, 2] < np.pi)


def test_paper_pf_runs_and_preserves_initial_yaw_lineage():
    trajectory = generate_paper_trajectory(dt=0.1, duration=2.0)
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    ranges = generate_uwb_ranges(trajectory.state[:, :2], auxiliary)
    rng_init = np.random.default_rng(5)
    particles = initialize_random_annulus_particles(ranges[0], auxiliary[0], 1000, rng_init, 0.1)
    rng_pf = np.random.default_rng(6)
    result = run_paper_bootstrap_pf(
        trajectory.increments,
        ranges,
        auxiliary,
        particles,
        rng_pf,
        sigma_uwb_m=0.05,
        truth_state=trajectory.state,
    )
    assert result.estimate.shape == trajectory.state.shape
    assert result.map_state.shape == trajectory.state.shape
    assert result.map_initial_yaw.shape == (len(trajectory.t),)
    assert result.lineage_initial_yaw_final.shape == (1000,)
    assert np.all(np.isfinite(result.estimate))
    assert np.all((result.correct_mode_mass >= 0.0) & (result.correct_mode_mass <= 1.0))
