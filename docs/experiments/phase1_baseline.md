# Phase 1 — Baseline implementation

## Question
Can a minimal, auditable 2-D simulator reproduce the expected progression from ideal dead reckoning, to drift under noisy increments, to drift correction using UWB range measurements in a conventional particle filter?

## Hypothesis
1. Noise-free dead reckoning reconstructs the deterministic ground truth to numerical precision.
2. Perturbed distance/heading increments accumulate position error.
3. A bootstrap PF using the same DR increments and a moving, globally localized UWB auxiliary node reduces position error.
4. Relative geometry affects local observability; observability singular values are recorded before attempting global unknown-pose reproduction.

## Scope
This experiment is deliberately a **debugging baseline**, not yet a numerical reproduction of Han et al.'s AACOPF. The target starts from a known approximate pose. The auxiliary trajectory is known. UWB noise is Gaussian. These restrictions isolate the basic fusion mechanism before global annular initialization, NLOS rejection, and AACOPF are introduced.

## Mathematical model
State: `x = [p_x, p_y, psi]^T`.

Motion model:

`p_x(k+1) = p_x(k) + ΔL(k) cos psi(k)`  
`p_y(k+1) = p_y(k) + ΔL(k) sin psi(k)`  
`psi(k+1) = wrap(psi(k) + Δpsi(k))`.

Range measurement to auxiliary node `a_k`:

`z_k = ||p_k-a_k||_2 + v_k`, with `v_k ~ N(0, sigma_UWB^2)`.

The bootstrap PF samples the DR process model, weights particles with the Gaussian range likelihood, and uses systematic resampling when the effective sample size falls below a threshold.

## Reproduction rule
All stochastic runs use an explicit seed. Configuration is stored in `configs/phase1.yaml`; result metrics and trajectories are written to `results/phase1/`.

## Interpretation policy
Do not interpret one seed as a final performance claim. The purpose of this first run is verification of implementation behavior. Multi-seed studies, update-rate sweeps, NLOS models, and energy claims belong to later experiments.
