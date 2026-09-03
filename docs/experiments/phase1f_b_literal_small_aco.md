# P1F-B — Literal-small all-pairs ACO mechanism validation

## Purpose

P1F-B implements and validates the audited ant-colony particle transition before any tuning or headline PF-vs-AACOPF comparison. The goal is deliberately narrower than P1F-C/P1F-D: establish that the repository implementation follows the P1E contract, expose what the transition actually does to a particle cloud, and measure the computational scaling of the literal all-pairs interpretation.

This is **not** an exact reproduction of Han et al. The paper does not report `alpha`, `beta`, `lambda`, the exact candidate rule, update ordering, or several numerical safeguards. P1F-B therefore follows the explicit interpretation frozen in `configs/phase1e_aacopf_interpretation.yaml`.

## Implemented transition

For source particle `i`, only particles with strictly higher normalized weight are candidate destinations:

`C_i = {j : w_j > w_i}`.

For candidate `j`, the score is

`s_ij = (w_j - w_i + eps_w)^alpha * (1 / (||p_i-p_j|| + eps_d))^beta`.

The candidate scores are normalized over `C_i`. The maximum-probability candidate is selected and particle `i` moves only when

`max_j P_ij > c_lambda / K_i`,

where `K_i = |C_i|`. The update is synchronous: all destination choices are computed from the immutable pre-transition cloud. When a particle moves, its complete `[x,y,phi]` state and initial-yaw lineage are copied from the selected destination. After the ACO transition, all particle weights are reset uniformly.

The current posterior estimate and MAP lineage are recorded **before** the ACO transition, matching the P1E ordering decision.

## Deterministic mechanism checks

The unit tests verify:

1. equal weights imply no strictly higher-weight candidates and therefore no movement;
2. the destination follows the documented weight-difference/inverse-distance score;
3. the transition is synchronous rather than sequential/cascading;
4. initial-yaw lineage is copied together with the destination state;
5. the normalized threshold can block otherwise valid moves;
6. equal-weight neighboring particles cannot select one another because the candidate relation is strictly higher weight;
7. the complete literal-small filtering sequence runs and records finite transition diagnostics.

These tests are the main P1F-B acceptance criterion. End-to-end accuracy is intentionally deferred to P1F-D.

## Controlled bimodal cloud

A secondary mechanism experiment uses a 12 s paper-state trajectory with exact DR increments, moving known auxiliary geometry, and Gaussian UWB noise (`sigma_UWB = 0.12 m`). The particle cloud has 400 particles split evenly between:

- a local mode around the correct initial pose;
- a wrong mode placed 90 degrees around the same first-range circle, with a 90 degree yaw offset.

Both modes are therefore plausible under the first range alone. This controlled setup is not the global-annulus problem used for the later PF-vs-AACOPF comparison. Its purpose is to observe whether the ACO transition reacts to sequential likelihood information and how aggressively it duplicates particles.

P1F-B uses `alpha = beta = 1` and `c_lambda = 0.5` only as a **provisional mechanism-check setting**. These values are not claimed to be optimal and are not frozen for P1F-D. P1F-C performs the actual sensitivity study.

Recorded diagnostics include pre-ACO effective sample size and correct-mode mass, post-ACO correct-mode fraction, moved-particle fraction, unique-parent fraction, destination multiplicity, transition-probability concentration, and transition runtime.

## Runtime scaling

The literal implementation forms dense pairwise distances and candidate scores. A single-transition benchmark is therefore run for particle counts

`N = {50, 100, 200, 400, 800, 1200}`.

For unique weights, the number of admissible directed higher-weight candidate scores is approximately `N(N-1)/2`, while the dense pairwise representation itself contains `N^2` entries. Runtime is measured after warmup and summarized by the median of repeated transitions. The empirical log-log runtime slope is reported only as an implementation diagnostic; the structural all-pairs complexity is already quadratic by construction.

## Acceptance / decision rule

P1F-B is complete when:

- all deterministic transition tests pass;
- the controlled filter executes without numerical failures and produces the required mechanism diagnostics;
- the runtime experiment confirms the practical cost growth of the dense all-pairs interpretation;
- no tuning conclusion is drawn from the provisional `c_lambda` value.

The next step is P1F-C: sensitivity of `alpha`, `beta`, and normalized threshold `c_lambda` on development seeds, followed by freezing a setting for the matched P1F-D comparison.
