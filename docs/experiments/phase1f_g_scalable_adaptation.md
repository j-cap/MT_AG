# Phase P1F-G — scalable bounded-candidate AACOPF adaptation

## Purpose

P1F-G is the final planned AACOPF experiment of Phase 1. It addresses the two main implementation boundaries established by P1F-D--P1F-F:

1. the literal all-pairs transition has quadratic cost and is not viable at the particle populations needed for sparse global initialization;
2. aggressive copying can destroy rare useful support or amplify a temporarily wrong likelihood mode.

The goal is therefore **not** to claim a new reproduction of Han et al. The P1F-G method is an explicit repository adaptation that asks whether the useful mode-amplification idea can be made computationally bounded and less destructive to particle diversity.

## Adapted transition

For each source particle, P1F-G scores only a global pool of the `K` highest-weight particles. Strictly higher-weight members of that pool remain eligible destinations and retain the audited score

`(w_j - w_i + eps_w)^alpha * (1 / (d_ij + eps_d))^beta`.

The P1F-C tuple remains frozen:

`(alpha, beta, c_lambda) = (0.5, 0, 2)`.

Because `beta=0`, the selected frozen transition is weight-driven; the bounded pool therefore targets the highest-weight destinations directly. Candidate scoring is `O(N K)` instead of `O(N^2)`.

The movement threshold deliberately preserves the literal normalization `lambda_i = c_lambda / K_i`, where `K_i` is the **total** number of strictly higher-weight particles in the full cloud, not the truncated candidate-pool size. `K_i` is obtained from the sorted weights in `O(N log N)`. This detail matters: using only the bounded candidate count would inflate the threshold and can suppress the ACO transition almost completely. The first implementation smoke campaign exposed exactly that failure mode, and the production P1F-G experiment therefore keeps the original threshold scale while truncating only destination scoring.

Two support/diversity guards are added after the score-threshold test:

- **move cap:** at most a fixed fraction of the source cloud may be replaced in one update;
- **destination capacity:** one destination can receive only a bounded fraction of the cloud as incoming copies.

Eligible moves are prioritized by the probability margin above the movement threshold. All accepted moves are synchronous and copy from the immutable pre-transition cloud.

## Development selection

A compact development sweep uses the same controlled two-mode benchmark as P1F-C with 400 particles and new seeds 400--407. The grid is:

- `K in {8,16,32}`;
- max move fraction in `{0.10,0.25,0.50}`;
- destination-capacity fraction in `{0.01,0.025,0.05}`.

The selected configuration must first meet a 90% correct-lock floor and at most 5% wrong-lock rate across balanced and minority-correct clouds. Among eligible settings the selection is safety-first: minimize catastrophic ancestry collapse, then dominant cloning, then candidate count and guard sizes. The workflow fails explicitly if no setting meets these predefined development criteria. Global random-annulus seeds are not used for tuning.

## Held-out controlled ablation

New seeds 500--519 compare:

- conventional PF;
- literal all-pairs AACOPF;
- bounded candidates only;
- bounded candidates + move cap only;
- bounded candidates + destination cap only;
- the fully guarded bounded transition.

This separates computational truncation from the two diversity guards.

## Held-out sparse global validation

The selected fully guarded transition is then evaluated on the complete random-annulus global-localization problem using new seeds 600--619, informative moving geometry, and particle budgets `Np = 1200, 5000, 10000`. Conventional PF receives the same initial cloud and UWB realization at every budget. A candidate-only bounded variant is retained at `Np=1200` as a direct support-preservation ablation.

No P1F-G parameter is retuned on these global seeds. The main metrics are strict full-pose convergence, late position/yaw RMSE, initial correct-region support, runtime, minimum unique-parent fraction, and maximum destination multiplicity.

## Runtime scaling

Single-transition runtime is measured for the selected bounded method at `Np = 400, 1200, 5000, 10000, 40000`. A literal reference is measured only up to `Np=1200`. Candidate-score count is also reported relative to the hypothetical dense `N^2` pair count.

## Claim boundary

P1F-G can support three distinct conclusions:

- **scalability:** bounded candidate scoring materially changes runtime/memory growth;
- **diversity safety:** the guards reduce ancestry collapse and destination monopolization;
- **global-localization effect:** only held-out random-annulus results can show whether preserving support translates into more robust full-pose recovery.

These claims must remain separate. A faster transition is not automatically a better estimator, and improved diversity is not automatically improved localization accuracy.
