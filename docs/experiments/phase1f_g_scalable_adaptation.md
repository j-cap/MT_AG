# Phase P1F-G — scalable bounded-candidate AACOPF adaptation

## Purpose

P1F-G is the final planned AACOPF experiment of Phase 1. It addresses the two structural boundaries established by P1F-D--P1F-F:

1. the literal all-pairs transition has quadratic cost and is not viable at the particle populations needed for sparse global initialization;
2. aggressive copying can destroy rare useful support or amplify a temporarily wrong likelihood mode.

P1F-G is an explicit repository adaptation rather than a literal reproduction of Han et al. The goal is to determine whether the useful mode-amplification idea can be made computationally bounded and genealogically safer.

## Adapted transition

For each source particle, P1F-G scores only a global pool of the `K` highest-weight particles. Strictly higher-weight members of that pool remain eligible destinations and use the audited score

`(w_j - w_i + eps_w)^alpha * (1 / (d_ij + eps_d))^beta`.

The P1F-C score parameters remain frozen:

`(alpha, beta, c_lambda) = (0.5, 0, 2)`.

Because `beta=0`, the selected transition is weight-driven. Candidate scoring is therefore reduced from dense `O(N^2)` scoring to `O(N K)`, plus sorting.

The movement threshold preserves the literal normalization `lambda_i = c_lambda / K_i`, where `K_i` is the **total** number of strictly higher-weight particles in the full cloud, not the truncated candidate-pool size. `K_i` is obtained from the sorted weights without constructing a dense pair matrix. An initial smoke implementation used the truncated candidate count in the denominator; this inflated the threshold and almost disabled movement, so that version was rejected before the final campaign.

Two diversity guards are applied after the score-threshold test:

- **move cap:** at most a fixed fraction of sources may be replaced in one update;
- **destination capacity:** one destination may receive only a bounded fraction of incoming copies.

All moves remain synchronous from the immutable pre-transition cloud.

## Development selection

The development sweep uses the controlled 400-particle two-mode benchmark and new seeds 400--407. The grid contains 27 settings:

- `K in {8,16,32}`;
- max move fraction in `{0.10,0.25,0.50}`;
- destination-capacity fraction in `{0.01,0.025,0.05}`.

A setting must first achieve at least 90% correct-mode lock and at most 5% wrong-mode lock across balanced and minority-correct clouds. Among eligible settings, the predefined rule minimizes catastrophic collapse, then dominant cloning, then candidate count and guard sizes.

The selected setting is

```
candidate_count = 8
max_move_fraction = 0.10
destination_capacity_fraction = 0.025
```

with the frozen score tuple `(0.5,0,2)`. On the 16 development runs it achieves 16/16 correct locks, 0/16 wrong locks, no catastrophic ancestry collapse, and no dominant-clone event.

## Held-out controlled ablation

Held-out seeds 500--519 compare conventional PF, literal AACOPF, candidate-only bounded AACOPF, each individual diversity guard, and the fully guarded method.

| method | correct lock | wrong lock | collapse | dominant clone | mean runtime [s] |
|---|---:|---:|---:|---:|---:|
| PF | 1.000 | 0.000 | 0.000 | 0.000 | 0.0235 |
| literal AACOPF | 0.950 | 0.000 | 0.100 | 0.750 | 0.9872 |
| candidate only | 0.450 | 0.525 | 1.000 | 1.000 | 0.0404 |
| move guard only | 1.000 | 0.000 | 0.000 | 0.000 | 0.0513 |
| destination guard only | 1.000 | 0.000 | 0.000 | 0.000 | 0.0634 |
| fully guarded | 1.000 | 0.000 | 0.000 | 0.000 | 0.0561 |

Candidate truncation alone is therefore unsafe. It converts the mode-amplification mechanism into an even stronger destination-monopolization operator: all held-out runs undergo catastrophic ancestry collapse and a dominant clone, and 52.5% finish in wrong-mode lock.

The diversity guards are the essential part of the adaptation. The fully guarded method achieves 40/40 correct locks with a minimum unique-parent fraction of 0.975. The destination-capacity guard alone produces the same held-out lock result, showing that incoming-copy capacity is the dominant safety mechanism in this benchmark; the move cap remains as an independent upper bound on transition aggressiveness.

## Held-out random-annulus global validation

The selected fully guarded method is then evaluated on new seeds 600--619 with the full random first-range annulus + random-yaw initialization, informative moving geometry, exact paper-level motion increments, and matched PF/AACOPF initial clouds and UWB realizations.

The final corrected campaign explicitly sets conventional-PF motion-increment noise to zero, matching the P1F-D paper-level reference.

| Np | method | strict pose success | mean late position RMSE [m] | mean late yaw RMSE [deg] | mean runtime [s] |
|---:|---|---:|---:|---:|---:|
| 1200 | PF | 1/20 | 14.812 | 46.687 | 0.255 |
| 1200 | guarded AACOPF | 0/20 | 10.699 | 32.638 | 0.570 |
| 1200 | candidate only | 0/20 | 36.096 | 111.593 | 0.501 |
| 5000 | PF | 2/20 | 10.275 | 28.331 | 0.802 |
| 5000 | guarded AACOPF | 2/20 | 8.000 | 22.672 | 1.797 |
| 10000 | PF | 3/20 | 10.017 | 29.514 | 1.525 |
| 10000 | guarded AACOPF | 1/20 | 7.043 | 19.501 | 3.439 |

Initial diagnostic correct-region support rises with population: 2.75 particles on average at 1200, 9.85 at 5000, and 19.0 at 10000.

The global result is deliberately separated from the controlled result. The guarded adaptation **does not improve strict full-pose acquisition**. At 1200 particles it loses the single PF success; at 5000 the methods tie at 2/20; at 10000 PF succeeds in 3/20 and guarded AACOPF in 1/20.

At the same time, guarded AACOPF reduces mean late position and yaw error at every tested population. This is retained as an error-shaping effect, not a convergence claim. The candidate-only ablation at 1200 particles confirms that computational truncation without diversity protection is actively harmful.

## Runtime scaling

Single-transition runtime is measured separately from complete filtering.

The bounded empirical log-log runtime slope over `Np={400,1200,5000,10000,40000}` is **0.926**. The literal reference over `Np={100,200,400,800,1200}` has slope **1.947**, consistent with dense quadratic growth over the measured range.

At `Np=40000` the bounded transition:

- has median transition runtime about 0.0281 s;
- scores 319,964 candidate pairs;
- would correspond to 1.6e9 dense pairs under an all-pairs implementation;
- therefore evaluates about 0.020% of the hypothetical dense pair count.

This supports the scalability claim independently of the localization result.

## P1F-G decision

P1F-G supports two claims and rejects a third:

- **supported — scalability:** bounded candidate scoring removes the dense all-pairs computational bottleneck;
- **supported — diversity safety:** explicit move/destination guards prevent the catastrophic cloning seen in literal and candidate-only AACOPF on the controlled held-out benchmark;
- **not supported — improved strict global acquisition:** the guarded transition does not outperform conventional PF in held-out random-annulus full-pose convergence.

The Phase-1 conclusion is therefore not that AACOPF solves global initialization. Instead, P1F-G identifies a safe and scalable form of weight-directed redistribution that can reduce late error when useful hypotheses survive, while leaving the fundamental finite-support/global-mode problem unresolved.

## Reproducibility

Canonical files:

- `configs/phase1f_g.yaml`
- `src/mt_ag/scalable_aacopf.py`
- `tests/test_scalable_aacopf.py`
- `experiments/run_phase1f_g_exact.py`
- `results/phase1/p1f_g/summary.json`
- `results/phase1/p1f_g/summary.md`

Final corrected GitHub Actions campaign:

- workflow run: `35068506066`
- artifact: `p1f-g-results`
- artifact ID: `10435397349`
- branch head used by the campaign: `63e2e1edf930222adad064f291e63a01cd4ab706`

Full per-run evidence remains in the Actions artifact; compact final evidence is version controlled.
