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

## Results

The controlled bimodal campaign concentrates at least 90% of the particles on the correct mode in 9/10 runs. Among those successful runs, the median concentration time is 0.60 s. The mean final correct-mode particle fraction and mean final pre-ACO posterior correct-mode mass are both 0.90.

These aggregate numbers hide the most important mechanism finding. The transition is **episodic and potentially catastrophic**: although only 4.42% of particles move on average over all updates, individual early updates can move more than 90% of the cloud and reduce the represented ancestry to only a few percent.

For representative seed 0, at `t=0.5 s` 94.25% of particles move and the unique-parent fraction falls to 6%. The pre-ACO posterior has already assigned 79.7% mass to the correct region. After the transition 66.25% of particles are in that region, and at `t=0.6 s` the cloud reaches 100% correct-mode particles.

For representative seed 4, at `t=0.3 s` 98.5% of particles move and the unique-parent fraction falls to 4%. Despite a 76.4% pre-ACO posterior mass on the correct region, the post-transition correct fraction is only 46.25%. At `t=0.5 s`, another 92.75% of particles move, one destination is copied by 302 particles, and the correct-mode fraction becomes exactly zero. This wrong-mode collapse is then irreversible in the tested horizon.

The practical lesson is that the current literal transition is not simply a mild diversity-preserving resampling alternative. With the provisional threshold it is a strong **mode-selection and cloning operator**. This is exactly why P1F-C must tune the threshold jointly with the score exponents before any PF-vs-AACOPF performance comparison.

A second diagnostic lesson is that effective sample size alone is insufficient after the audited uniform weight reset. A cloud may have nearly uniform weights and therefore high `N_eff` while having very few distinct parents. Unique-parent fraction and destination multiplicity must remain first-class diagnostics.

The runtime benchmark gives median transition times of approximately 0.315 ms, 0.633 ms, 2.08 ms, 7.03 ms, 33.05 ms, and 73.32 ms for 50, 100, 200, 400, 800, and 1200 particles respectively. The empirical log-log slope is approximately 1.76 over this finite range, while the structural dense all-pairs implementation remains `O(N^2)`.

At the `N_p=160000` population implied by the paper's `N=M=400`, a single dense float64 pairwise matrix would contain `2.56e10` entries and require roughly 204.8 GB. Since the implementation requires multiple such arrays, the literal dense interpretation is computationally infeasible at the paper's nominal population. A naive quadratic timing extrapolation from 1200 particles would also imply roughly 22 minutes per update; this extrapolation is illustrative only because memory becomes prohibitive first.

## Acceptance / decision

P1F-B is complete:

- deterministic mechanism tests pass;
- the controlled filter executes and provides the required diagnostics;
- the ACO transition can rapidly amplify the correct mode but can also prematurely destroy it;
- explicit diversity diagnostics are necessary in addition to effective sample size;
- the literal dense all-pairs interpretation is unsuitable for the paper's nominal particle count and for the intended embedded platform;
- no tuning conclusion is drawn from `c_lambda=0.5`.

The next step is P1F-C: sensitivity of `alpha`, `beta`, and normalized threshold `c_lambda` on development seeds. The P1F-C selection criterion must penalize both wrong-mode lock and catastrophic ancestry collapse, not just final position error. Only after that study should a setting be frozen for P1F-D.

Provenance: workflow run `33751243293`, artifact `9891624402`, commit `c24f751891a6e57fe6e75d86c63340854141066a`.
