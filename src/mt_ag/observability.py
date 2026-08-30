import numpy as np


def local_observability_matrix(states, increments, auxiliary_positions, horizon=None):
    """Finite-horizon local linearized observability matrix for x=[px, py, psi]."""
    states = np.asarray(states, dtype=float)
    inc = np.asarray(increments, dtype=float)
    aux = np.asarray(auxiliary_positions, dtype=float)
    K = min(len(states), horizon or len(states))
    Phi = np.eye(3)
    rows = []
    for k in range(K):
        dx = states[k, 0] - aux[k, 0]
        dy = states[k, 1] - aux[k, 1]
        d = np.hypot(dx, dy)
        if d < 1e-12:
            raise ValueError("Target and auxiliary coincide; range Jacobian is singular")
        H = np.array([[dx / d, dy / d, 0.0]])
        rows.append(H @ Phi)
        if k < K - 1:
            dL = inc[k, 0]
            psi = states[k, 2]
            F = np.array([[1.0, 0.0, -dL*np.sin(psi)],
                          [0.0, 1.0, dL*np.cos(psi)],
                          [0.0, 0.0, 1.0]])
            Phi = F @ Phi
    return np.vstack(rows)


def singular_values(states, increments, auxiliary_positions, horizon=None):
    O = local_observability_matrix(states, increments, auxiliary_positions, horizon)
    return np.linalg.svd(O, compute_uv=False)
