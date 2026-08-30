import numpy as np


def wrap_angle(angle):
    """Wrap angle(s) to [-pi, pi)."""
    return (np.asarray(angle) + np.pi) % (2.0 * np.pi) - np.pi


def ranges(points, anchor):
    points = np.asarray(points, dtype=float)
    anchor = np.asarray(anchor, dtype=float)
    return np.linalg.norm(points - anchor, axis=-1)
