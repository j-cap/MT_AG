# Phase P1F-D — matched conventional PF versus frozen literal AACOPF

## Decision

P1F-D does **not** show a convergence-rate advantage for the frozen literal AACOPF under full random-annulus global initialization at equal small particle budgets. At the primary 1,200-particle budget, conventional PF converges in 2/20 runs and AACOPF in 1/20; the successful runs are disjoint. However, AACOPF substantially reduces mean late position/yaw RMSE at 1,200 particles, indicating that it often selects a closer global mode without reliably entering the strict terminal convergence region.

This is a crucial transfer result: the P1F-C tuple was tuned on deliberately well-populated bimodal clouds, whereas random annulus initialization provides only a handful of correct-region particles. The literal ACO transition can amplify represented hypotheses but cannot create missing or insufficiently represented state support.

## Matched design

- state `[x,y,phi]`, exact `[Delta L,Delta phi]`, Equation-(3) pre-turn propagation;
- informative moving globally known auxiliary;
- Gaussian UWB noise `sigma=0.12 m`;
- identical random first-range annulus + random-yaw cloud for PF and AACOPF;
- frozen AACOPF `(alpha,beta,c_lambda)=(0.5,0,2)`;
- 20 unseen seeds at each of 400, 800, and 1,200 particles.

## Aggregate results

| Np | Initial support >0 | PF success | AACOPF success | PF late pos RMSE | AACOPF late pos RMSE | PF late yaw RMSE | AACOPF late yaw RMSE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 400 | 50% | 0% | 0% | 17.678 m | 16.894 m | 56.45° | 50.32° |
| 800 | 95% | 5% | 0% | 21.371 m | 19.080 m | 64.89° | 57.33° |
| 1200 | 95% | 10% | 5% | 19.945 m | 12.782 m | 66.95° | 33.99° |

At 1,200 particles, the paired outcome counts are: **0 both**, **1 AACOPF-only**, **2 PF-only**, and **17 neither**. Therefore the strict convergence result does not support an AACOPF success-rate improvement.

Nevertheless, AACOPF has lower late position RMSE in **65%** of matched 1,200-particle runs and lower late yaw RMSE in **70%**. The mean paired differences are -7.162 m in position RMSE and -32.96° in yaw RMSE (AACOPF minus PF).

## Initial-support diagnosis

- **400 particles:** mean correct-region support 0.80; positive support in 50% of runs.
- **800 particles:** mean correct-region support 2.05; positive support in 95% of runs.
- **1200 particles:** mean correct-region support 2.35; positive support in 95% of runs.

The primary 1,200-particle cloud contains only **2.35 correct-region particles on average**, with a range of 0–5. This is fundamentally different from P1F-C, where the correct mode initially contained 80 particles in the 20/80 case and 200 particles in the 50/50 case. The P1F-C tuning result therefore does not transfer directly to sparse global support.

At 1,200 particles, success conditioned on positive initial support is 2/19 for PF and 1/19 for AACOPF. For runs with three or more initial correct-region particles, each method succeeds in 1/9 runs. Thus merely having a few correct particles is necessary but far from sufficient.

## Diversity and runtime

| Np | AACOPF catastrophic collapse | AACOPF dominant clone | PF runtime/run | AACOPF runtime/run |
|---:|---:|---:|---:|---:|
| 400 | 5% | 85% | 0.130 s | 6.19 s |
| 800 | 30% | 80% | 0.190 s | 21.83 s |
| 1200 | 10% | 85% | 0.248 s | 47.08 s |

The frozen setting reduces catastrophic-collapse incidence compared with the provisional P1F-B point, but **dominant cloning remains present in 80–85% of runs**. At 1,200 particles, the literal AACOPF is about **190× slower** than the conventional PF in the current Python implementation. This reinforces that P1F-G is required before any deployment-oriented interpretation.

## Interpretation

1. **P1F-C solved the controlled mode-selection tuning problem, not the complete global-initialization problem.** The tuned operator works when the correct mode has substantial initial population, but random annulus initialization leaves it severely underrepresented.
2. **Literal AACOPF shows a useful error-shaping effect at 1,200 particles.** It often ends closer to the true mode, approximately halving mean late yaw RMSE and reducing mean late position RMSE by roughly one third, but this does not translate into a higher strict convergence fraction.
3. **The result does not justify retuning on P1F-D seeds.** The frozen tuple must remain frozen to preserve the experiment's validity. Any algorithmic change motivated by this outcome belongs to a new experiment.
4. **P1F-E remains necessary.** Rank-deficient geometry negative controls must verify that AACOPF does not appear to overcome missing observability.
5. **P1F-G becomes more important, not less.** A scalable adaptation should address both computational cost and sparse-support preservation/rejuvenation rather than simply approximating the same aggressive all-pairs cloning rule.

## Reproducibility

- workflow run: `33881931650`
- artifact: `p1f-d-results`, ID `9941206807`
- configuration: `configs/phase1f_d.yaml`
- runner: `experiments/run_phase1f_d.py`
- full raw results retained in the workflow artifact; aggregate evidence is version controlled here.
