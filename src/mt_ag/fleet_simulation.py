from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from .geometry import rotation_2d, wrap_angle


@dataclass(frozen=True)
class FleetTrajectory:
    t: np.ndarray
    state: np.ndarray  # shape [time, node, 5] = [x,y,vx,vy,psi]
    ideal_imu: np.ndarray  # shape [time-1, node, 3] = [ax_b,ay_b,omega]
    dt: float

    @property
    def n_nodes(self) -> int:
        return int(self.state.shape[1])


def all_pairs(n_nodes: int) -> list[tuple[int, int]]:
    return list(combinations(range(n_nodes), 2))


def _kinematics(node: int, t: np.ndarray) -> tuple[np.ndarray, ...]:
    """Analytic smooth benchmark paths and their first two derivatives."""
    t = np.asarray(t, dtype=float)

    if node == 0:
        w = 2.0 * np.pi / 52.0
        theta = w * t + 0.15
        x = -2.0 + 6.0 * np.cos(theta)
        y = 4.0 * np.sin(theta)
        vx = -6.0 * w * np.sin(theta)
        vy = 4.0 * w * np.cos(theta)
        ax = -6.0 * w**2 * np.cos(theta)
        ay = -4.0 * w**2 * np.sin(theta)
    elif node == 1:
        w = 2.0 * np.pi / 57.0
        theta = -w * t + 2.15
        x = 3.0 + 5.2 * np.cos(theta)
        y = 2.0 + 3.6 * np.sin(theta)
        vx = 5.2 * w * np.sin(theta)
        vy = -3.6 * w * np.cos(theta)
        ax = -5.2 * w**2 * np.cos(theta)
        ay = -3.6 * w**2 * np.sin(theta)
    elif node == 2:
        w = 2.0 * np.pi / 60.0
        theta = w * t - 0.4
        x = -1.0 + 5.0 * np.sin(theta)
        y = -2.0 + 3.0 * np.sin(2.0 * theta)
        vx = 5.0 * w * np.cos(theta)
        vy = 6.0 * w * np.cos(2.0 * theta)
        ax = -5.0 * w**2 * np.sin(theta)
        ay = -12.0 * w**2 * np.sin(2.0 * theta)
    elif node == 3:
        w = 2.0 * np.pi / 64.0
        theta = w * t + 0.85
        x = 5.0 + 4.0 * np.sin(theta)
        y = -4.0 + 3.2 * np.sin(2.0 * theta + 0.55)
        vx = 4.0 * w * np.cos(theta)
        vy = 6.4 * w * np.cos(2.0 * theta + 0.55)
        ax = -4.0 * w**2 * np.sin(theta)
        ay = -12.8 * w**2 * np.sin(2.0 * theta + 0.55)
    else:
        raise ValueError("P2A currently defines exactly four benchmark paths")

    return x, y, vx, vy, ax, ay


def generate_four_node_fleet(dt: float = 0.01, duration: float = 60.0) -> FleetTrajectory:
    """Generate smooth global truth first, then synthesize ideal planar IMU signals."""
    t = np.arange(0.0, duration + 0.5 * dt, dt)
    n_nodes = 4
    state = np.zeros((len(t), n_nodes, 5), dtype=float)
    ideal_imu = np.zeros((len(t) - 1, n_nodes, 3), dtype=float)

    for node in range(n_nodes):
        x, y, vx, vy, _, _ = _kinematics(node, t)
        state[:, node, 0] = x
        state[:, node, 1] = y
        state[:, node, 2] = vx
        state[:, node, 3] = vy
        state[:, node, 4] = wrap_angle(np.arctan2(vy, vx))

        t_mid = 0.5 * (t[:-1] + t[1:])
        _, _, vx_mid, vy_mid, ax_mid, ay_mid = _kinematics(node, t_mid)
        psi_mid = np.arctan2(vy_mid, vx_mid)
        speed_sq = vx_mid**2 + vy_mid**2
        omega = (vx_mid * ay_mid - vy_mid * ax_mid) / speed_sq

        for k in range(len(t_mid)):
            a_nav = np.array([ax_mid[k], ay_mid[k]], dtype=float)
            ideal_imu[k, node, :2] = rotation_2d(-psi_mid[k]) @ a_nav
            ideal_imu[k, node, 2] = omega[k]

    return FleetTrajectory(t=t, state=state, ideal_imu=ideal_imu, dt=dt)


def pairwise_distances(positions: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 3 or positions.shape[2] != 2:
        raise ValueError("positions must have shape [time, node, 2]")
    return {
        pair: np.linalg.norm(positions[:, pair[0]] - positions[:, pair[1]], axis=1)
        for pair in all_pairs(positions.shape[1])
    }
