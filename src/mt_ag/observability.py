import numpy as np


def observability_matrix_planar_imu(states, ideal_imu, auxiliary_positions, dt):
    """Finite-horizon local observability matrix for [p_x,p_y,v_x,v_y,psi]."""
    states = np.asarray(states, dtype=float)
    ideal_imu = np.asarray(ideal_imu, dtype=float)
    auxiliary_positions = np.asarray(auxiliary_positions, dtype=float)
    phi = np.eye(5)
    rows = []
    for k in range(len(states)):
        displacement = states[k, :2] - auxiliary_positions[k]
        distance = np.linalg.norm(displacement)
        if distance < 1e-12:
            raise ValueError("Target and auxiliary node coincide; range Jacobian undefined")
        h = np.zeros((1, 5), dtype=float)
        h[0, :2] = displacement / distance
        rows.append(h @ phi)
        if k < len(ideal_imu):
            a_x_b, a_y_b, omega_z = ideal_imu[k]
            psi_mid = states[k, 4] + 0.5 * omega_z * dt
            c = np.cos(psi_mid)
            s = np.sin(psi_mid)
            da_dpsi = np.array([-s * a_x_b - c * a_y_b, c * a_x_b - s * a_y_b])
            f = np.eye(5)
            f[0, 2] = dt
            f[1, 3] = dt
            f[0, 4] = 0.5 * dt**2 * da_dpsi[0]
            f[1, 4] = 0.5 * dt**2 * da_dpsi[1]
            f[2, 4] = dt * da_dpsi[0]
            f[3, 4] = dt * da_dpsi[1]
            phi = f @ phi
    return np.vstack(rows)


def singular_values_planar_imu(states, ideal_imu, auxiliary_positions, dt):
    return np.linalg.svd(observability_matrix_planar_imu(states, ideal_imu, auxiliary_positions, dt), compute_uv=False)
