import numpy as np


def _local_transition(states, ideal_imu, k, dt):
    a_x_b, a_y_b, omega_z = ideal_imu[k]
    psi_mid = states[k, 4] + 0.5 * omega_z * dt
    c = np.cos(psi_mid)
    s = np.sin(psi_mid)
    da_dpsi = np.array(
        [-s * a_x_b - c * a_y_b, c * a_x_b - s * a_y_b], dtype=float
    )
    f = np.eye(5)
    f[0, 2] = dt
    f[1, 3] = dt
    f[0, 4] = 0.5 * dt**2 * da_dpsi[0]
    f[1, 4] = 0.5 * dt**2 * da_dpsi[1]
    f[2, 4] = dt * da_dpsi[0]
    f[3, 4] = dt * da_dpsi[1]
    return f


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
            phi = _local_transition(states, ideal_imu, k, dt) @ phi
    return np.vstack(rows)


def singular_values_planar_imu(states, ideal_imu, auxiliary_positions, dt):
    return np.linalg.svd(
        observability_matrix_planar_imu(states, ideal_imu, auxiliary_positions, dt),
        compute_uv=False,
    )


def observability_history_planar_imu(
    states,
    ideal_imu,
    auxiliary_positions,
    dt,
    *,
    rank_rtol=1e-9,
):
    """Return cumulative finite-horizon observability diagnostics at every sample.

    Each cumulative observability matrix is evaluated by direct SVD. This avoids
    the loss of numerical rank information that can occur when forming O^T O,
    whose condition number is squared. ``sigma_ratio`` is sigma_min/sigma_max
    and is used only as a local finite-horizon conditioning diagnostic; it is
    neither an application threshold nor a complete nonlinear observability proof.
    """
    states = np.asarray(states, dtype=float)
    ideal_imu = np.asarray(ideal_imu, dtype=float)
    auxiliary_positions = np.asarray(auxiliary_positions, dtype=float)
    phi = np.eye(5)
    rows = []
    n = len(states)
    singular_values = np.zeros((n, 5), dtype=float)
    rank = np.zeros(n, dtype=int)

    for k in range(n):
        displacement = states[k, :2] - auxiliary_positions[k]
        distance = np.linalg.norm(displacement)
        if distance < 1e-12:
            raise ValueError("Target and auxiliary node coincide; range Jacobian undefined")
        h = np.zeros((1, 5), dtype=float)
        h[0, :2] = displacement / distance
        rows.append((h @ phi).ravel())
        matrix = np.asarray(rows, dtype=float)
        sigma_available = np.linalg.svd(matrix, compute_uv=False)
        sigma = np.zeros(5, dtype=float)
        sigma[: len(sigma_available)] = sigma_available
        singular_values[k] = sigma
        tol = rank_rtol * sigma[0] if sigma[0] > 0.0 else 0.0
        rank[k] = int(np.sum(sigma > tol))
        if k < len(ideal_imu):
            phi = _local_transition(states, ideal_imu, k, dt) @ phi

    sigma_max = singular_values[:, 0]
    sigma_min = singular_values[:, -1]
    ratio = np.divide(
        sigma_min,
        sigma_max,
        out=np.zeros_like(sigma_min),
        where=sigma_max > 0.0,
    )
    return {
        "singular_values": singular_values,
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "sigma_ratio": ratio,
        "rank": rank,
    }


def _paper_local_transition(state, increment, convention):
    """Jacobian of the three-state DR propagation with respect to prior state."""
    delta_l, delta_phi = np.asarray(increment, dtype=float)
    phi = float(state[2])
    if convention == "pre_turn":
        heading = phi
    elif convention == "post_turn":
        heading = phi + delta_phi
    else:
        raise ValueError("convention must be 'pre_turn' or 'post_turn'")
    f = np.eye(3)
    f[0, 2] = -delta_l * np.sin(heading)
    f[1, 2] = delta_l * np.cos(heading)
    return f


def observability_matrix_paper_state(
    states,
    increments,
    auxiliary_positions,
    *,
    convention="pre_turn",
):
    """Finite-horizon local range observability matrix for [x,y,phi].

    The matrix linearizes the exact paper-level dead-reckoning recursion around
    the supplied trajectory. It is a local finite-horizon diagnostic, not a
    nonlinear global-observability proof.
    """
    states = np.asarray(states, dtype=float)
    increments = np.asarray(increments, dtype=float)
    auxiliary_positions = np.asarray(auxiliary_positions, dtype=float)
    if states.ndim != 2 or states.shape[1] != 3:
        raise ValueError("states must have shape (N, 3)")
    if increments.shape != (len(states) - 1, 2):
        raise ValueError("increments must have shape (N-1, 2)")
    if auxiliary_positions.shape != (len(states), 2):
        raise ValueError("auxiliary_positions must have shape (N, 2)")

    phi_total = np.eye(3)
    rows = []
    for k in range(len(states)):
        displacement = states[k, :2] - auxiliary_positions[k]
        distance = np.linalg.norm(displacement)
        if distance < 1e-12:
            raise ValueError("Target and auxiliary node coincide; range Jacobian undefined")
        h = np.zeros((1, 3), dtype=float)
        h[0, :2] = displacement / distance
        rows.append(h @ phi_total)
        if k < len(increments):
            phi_total = (
                _paper_local_transition(states[k], increments[k], convention) @ phi_total
            )
    return np.vstack(rows)


def observability_history_paper_state(
    states,
    increments,
    auxiliary_positions,
    *,
    convention="pre_turn",
    rank_rtol=1e-9,
):
    """Return cumulative SVD/rank diagnostics for the three-state paper model."""
    states = np.asarray(states, dtype=float)
    increments = np.asarray(increments, dtype=float)
    auxiliary_positions = np.asarray(auxiliary_positions, dtype=float)
    phi_total = np.eye(3)
    rows = []
    singular_values = np.zeros((len(states), 3), dtype=float)
    rank = np.zeros(len(states), dtype=int)

    for k in range(len(states)):
        displacement = states[k, :2] - auxiliary_positions[k]
        distance = np.linalg.norm(displacement)
        if distance < 1e-12:
            raise ValueError("Target and auxiliary node coincide; range Jacobian undefined")
        h = np.zeros((1, 3), dtype=float)
        h[0, :2] = displacement / distance
        rows.append((h @ phi_total).ravel())
        sigma_available = np.linalg.svd(np.asarray(rows), compute_uv=False)
        sigma = np.zeros(3, dtype=float)
        sigma[: len(sigma_available)] = sigma_available
        singular_values[k] = sigma
        tol = rank_rtol * sigma[0] if sigma[0] > 0.0 else 0.0
        rank[k] = int(np.sum(sigma > tol))
        if k < len(increments):
            phi_total = (
                _paper_local_transition(states[k], increments[k], convention) @ phi_total
            )

    sigma_max = singular_values[:, 0]
    sigma_min = singular_values[:, -1]
    ratio = np.divide(
        sigma_min,
        sigma_max,
        out=np.zeros_like(sigma_min),
        where=sigma_max > 0.0,
    )
    return {
        "singular_values": singular_values,
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "sigma_ratio": ratio,
        "rank": rank,
    }


def relative_geometry_metrics(target_positions, auxiliary_positions):
    """Summarize range and line-of-sight variation for one geometry."""
    target = np.asarray(target_positions, dtype=float)
    auxiliary = np.asarray(auxiliary_positions, dtype=float)
    relative = auxiliary - target
    ranges = np.linalg.norm(relative, axis=1)
    if np.any(ranges < 1e-12):
        raise ValueError("Target and auxiliary node coincide")
    bearing = np.unwrap(np.arctan2(relative[:, 1], relative[:, 0]))
    bearing_change = np.diff(bearing)
    return {
        "range_min_m": float(np.min(ranges)),
        "range_max_m": float(np.max(ranges)),
        "range_span_m": float(np.ptp(ranges)),
        "range_std_m": float(np.std(ranges)),
        "bearing_span_deg": float(np.rad2deg(np.ptp(bearing))),
        "bearing_total_variation_deg": float(
            np.rad2deg(np.sum(np.abs(bearing_change)))
        ),
    }
