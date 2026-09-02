# Phase 1 baseline — first run summary

Configuration: `configs/phase1.yaml`, seed 42.

| Method | Position RMSE [m] | p95 [m] | Final error [m] | Heading RMSE [deg] |
|---|---:|---:|---:|---:|
| Ideal DR | 0.000 | 0.000 | 0.000 | 0.000 |
| Noisy DR | 1.058 | 2.307 | 2.820 | 7.000 |
| PF + UWB | 0.338 | 0.622 | 0.682 | 2.536 |

The PF reduces position RMSE by about 68% relative to noisy DR in this verification run. This is not yet a statistical claim because it is a single seeded realization.

## Observability diagnostic

- Stationary auxiliary singular values: `[176.166, 10.074, 4.94e-15]`; ratio `sigma_min/sigma_max = 2.80e-17`.
- Moving auxiliary singular values: `[110.394, 11.925, 6.150]`; ratio `sigma_min/sigma_max = 5.57e-2`.

The stationary case is numerically rank-deficient for the chosen finite-horizon linearization, whereas the moving case is full rank and much better conditioned. This supports using auxiliary motion/geometry as an explicit experimental factor in the unknown-pose reproduction.

## Decision

The baseline pipeline is accepted as Milestone P1A. Do **not** interpret it as an AACOPF reproduction. Next: multi-seed check, then unknown-pose annular initialization and stationary-vs-moving convergence.
