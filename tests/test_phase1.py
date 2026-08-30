import numpy as np
from mt_ag.geometry import wrap_angle
from mt_ag.simulation import generate_curved_trajectory, integrate_dr
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.particle_filter import systematic_resample


def test_angle_wrap():
    assert np.isclose(wrap_angle(np.pi), -np.pi)
    assert np.isclose(wrap_angle(3*np.pi), -np.pi)


def test_ideal_dr_reconstructs_truth():
    tr = generate_curved_trajectory(dt=0.1, duration=10.0)
    rec = integrate_dr(tr.state[0], tr.increments)
    assert np.max(np.abs(rec - tr.state)) < 1e-12


def test_ideal_uwb_is_euclidean_range():
    p = np.array([[0., 0.], [3., 4.]])
    a = np.array([[0., 0.], [0., 0.]])
    z = generate_uwb_ranges(p, a)
    assert np.allclose(z, [0., 5.])


def test_systematic_resampling_preserves_particle_count():
    rng = np.random.default_rng(1)
    w = np.array([0.05, 0.15, 0.8])
    idx = systematic_resample(w, rng)
    assert len(idx) == len(w)
    assert np.all((idx >= 0) & (idx < len(w)))
