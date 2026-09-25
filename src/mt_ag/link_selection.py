"""Causal link-selection policies for the four-node P2C ranging campaign."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .fleet_simulation import all_pairs


def exchange_capacity(n_time: int, stride: int, budget: int, n_links: int) -> np.ndarray:
    """Spread a fixed number of exchanges over common maximum-rate opportunities."""
    if stride < 1 or (n_time - 1) % stride or budget < 0:
        raise ValueError("Invalid opportunity grid or budget")
    n_slots = (n_time - 1) // stride
    if budget > n_slots * n_links:
        raise ValueError("Budget exceeds available pairwise opportunities")
    capacity = np.zeros(n_time, dtype=int)
    for slot in range(1, n_slots + 1):
        capacity[slot * stride] = slot * budget // n_slots - (slot - 1) * budget // n_slots
    return capacity


def eligible_mask(capacity: np.ndarray, n_nodes: int) -> np.ndarray:
    mask = np.zeros((len(capacity), n_nodes, n_nodes), dtype=bool)
    for i, j in all_pairs(n_nodes):
        mask[capacity > 0, i, j] = True
        mask[capacity > 0, j, i] = True
    return mask


def baseline_mask(
    capacity: np.ndarray,
    n_nodes: int,
    policy: str,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Fixed cyclic and seeded random controls with the same per-tick capacity."""
    if policy not in ("cyclic", "random"):
        raise ValueError("Unrecognized baseline policy")
    if policy == "random" and rng is None:
        raise ValueError("Random baseline requires an independent RNG")
    pairs = all_pairs(n_nodes)
    mask = np.zeros((len(capacity), n_nodes, n_nodes), dtype=bool)
    next_link = 0
    for t, count in enumerate(capacity):
        if count > len(pairs):
            raise ValueError("Per-tick capacity exceeds number of links")
        if policy == "cyclic":
            chosen = [pairs[(next_link + offset) % len(pairs)] for offset in range(count)]
            next_link += count
        else:
            chosen = [pairs[index] for index in rng.choice(len(pairs), count, replace=False)]
        for i, j in chosen:
            mask[t, i, j] = mask[t, j, i] = True
    return mask


def range_jacobian(state_flat: np.ndarray, pair: tuple[int, int]) -> np.ndarray:
    """Range Jacobian at the predicted state, with the joint-EKF state ordering."""
    n_nodes = len(state_flat) // 5
    i, j = pair
    difference = state_flat[5 * i : 5 * i + 2] - state_flat[5 * j : 5 * j + 2]
    distance = float(np.linalg.norm(difference))
    if distance < 1e-9:
        raise ValueError("Predicted pair distance is too small to choose a link")
    unit = difference / distance
    h = np.zeros(5 * n_nodes)
    h[5 * i : 5 * i + 2] = unit
    h[5 * j : 5 * j + 2] = -unit
    return h


@dataclass
class PredictiveLinkSelector:
    """Greedy geometry or EKF relative-position-variance reduction selector.

    Both policies use only predicted state/covariance and prior decisions. The
    current range value and ground truth are never available to the selector.
    """

    policy: str
    capacity: np.ndarray
    sigma_range_m: float
    dt: float
    n_nodes: int
    geometry_memory_s: float = 2.0
    geometry_prior: float = 1.0
    records: list[tuple[int, int, int, float]] = field(default_factory=list)
    _geometry_gram: np.ndarray = field(init=False, repr=False)
    _last_time_index: int = field(default=0, init=False, repr=False)

    def __post_init__(self) -> None:
        if self.policy not in ("geometry", "information"):
            raise ValueError("Unknown predictive policy")
        if self.geometry_memory_s <= 0 or self.geometry_prior <= 0:
            raise ValueError("Geometry memory and prior must be positive")
        self._geometry_gram = np.zeros((2 * self.n_nodes, 2 * self.n_nodes))

    def __call__(
        self,
        time_index: int,
        state_flat: np.ndarray,
        covariance: np.ndarray,
        eligible: tuple[tuple[int, int], ...],
    ) -> list[tuple[int, int]]:
        count = int(self.capacity[time_index])
        if count > len(eligible):
            raise ValueError("Selector capacity exceeds available links")
        if self.policy == "geometry":
            decay = np.exp(-(time_index - self._last_time_index) * self.dt /
                           self.geometry_memory_s)
            self._geometry_gram *= decay
            self._last_time_index = time_index

        pos_indices = np.array([5 * node + axis
                                for node in range(self.n_nodes) for axis in range(2)])
        candidates = {pair: range_jacobian(state_flat, pair) for pair in eligible}
        chosen = []
        predicted_cov = covariance.copy()
        for _ in range(count):
            scores = {}
            for pair, h in candidates.items():
                if self.policy == "geometry":
                    direction = h[pos_indices] / self.sigma_range_m
                    gram_regularized = self._geometry_gram + self.geometry_prior * np.eye(
                        2 * self.n_nodes
                    )
                    scores[pair] = float(direction @ np.linalg.solve(
                        gram_regularized, direction
                    ))
                else:
                    ph = predicted_cov @ h
                    innovation_var = float(h @ ph + self.sigma_range_m**2)
                    centered = ph[pos_indices].reshape(self.n_nodes, 2)
                    centered = centered - np.mean(centered, axis=0)
                    scores[pair] = float(np.sum(centered**2) /
                                         (self.n_nodes * innovation_var))
            # Candidate insertion order is canonical and resolves exact ties.
            best = max(scores, key=scores.get)
            h = candidates.pop(best)
            chosen.append(best)
            self.records.append((time_index, *best, scores[best]))
            if self.policy == "geometry":
                direction = h[pos_indices] / self.sigma_range_m
                self._geometry_gram += np.outer(direction, direction)
            else:
                ph = predicted_cov @ h
                innovation_var = float(h @ ph + self.sigma_range_m**2)
                predicted_cov -= np.outer(ph, ph) / innovation_var
        return chosen
