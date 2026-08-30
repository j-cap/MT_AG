import numpy as np


def generate_uwb_ranges(target_positions, auxiliary_positions, rng=None, sigma=0.0):
    target_positions = np.asarray(target_positions, dtype=float)
    auxiliary_positions = np.asarray(auxiliary_positions, dtype=float)
    ideal = np.linalg.norm(target_positions - auxiliary_positions, axis=1)
    if sigma <= 0.0:
        return ideal.copy()
    if rng is None:
        raise ValueError("rng is required for noisy measurements")
    return ideal + rng.normal(0.0, sigma, len(ideal))
