# Phase P1F-C — AACOPF parameter-tuning results

## Decision

The predeclared safety-aware tuning and held-out validation select and freeze the repository setting

`alpha = 0.5`, `beta = 0.0`, `c_lambda = 2.0`

for P1F-D. These are repository-selected values; Han et al. do not report numerical values for these parameters.

The selected tuple passes the predeclared validation freeze criterion. It is therefore suitable for the next matched conventional-PF versus literal-small-AACOPF comparison. This does **not** imply that the diversity problem is solved or that the literal mechanism is deployment-ready.

## Development result

The full sweep contains 60 tuples:

- `alpha in {0.5, 1, 2}`;
- `beta in {0, 0.5, 1, 2}`;
- `c_lambda in {0.5, 1, 2, 4, 8}`.

Each tuple was evaluated on 16 matched development runs: eight balanced 50/50 correct/wrong clouds and eight harder minority-correct 20/80 clouds, each with 400 particles.

For the selected tuple:

| Metric | Development result |
|---|---:|
| Correct-mode terminal lock | 15/16 = 93.75% |
| Wrong-mode terminal lock | 0/16 = 0% |
| Catastrophic ancestry collapse | 1/16 = 6.25% |
| Dominant-clone event | 10/16 = 62.5% |
| Mean final correct-mode fraction | 0.9931 |
| Mean minimum unique-parent fraction | 0.4023 |
| Mean moved-particle fraction | 0.0399 |

The next two safety-ranked development settings also had `beta=0`: `(alpha,beta,c_lambda)=(1,0,4)` and `(2,0,8)`. This is evidence that the inverse-distance term is not required for good controlled mode selection in this experiment and can make the literal transition more aggressive. It is not a general statement about all AACOPF applications.

## Held-out validation

The selected tuple was frozen before validation and then evaluated on 40 unseen runs: 20 balanced and 20 minority-correct.

| Metric | Validation result |
|---|---:|
| Correct-mode terminal lock | **40/40 = 100%** |
| Wrong-mode terminal lock | **0/40 = 0%** |
| Catastrophic ancestry collapse | **9/40 = 22.5%** |
| Dominant-clone event | **34/40 = 85%** |
| Mean final correct-mode fraction | 0.999625 |
| Mean minimum unique-parent fraction | 0.3414 |
| Mean moved-particle fraction | 0.0392 |
| Mean full-horizon position RMSE | 2.846 m |
| Mean full-horizon yaw RMSE | 21.812 deg |

Both scenario subsets reach 20/20 correct lock. Catastrophic-collapse incidence is 5/20 in the balanced case and 4/20 in the minority-correct case.

The full-horizon RMSE values include the deliberately ambiguous initial interval and should not be used as final localization-performance claims. P1F-D owns the matched accuracy comparison.

## Comparison with the P1F-B mechanism-check point

The provisional P1F-B tuple `(alpha,beta,c_lambda)=(1,1,0.5)` is clearly too aggressive on the P1F-C development set:

- correct lock: 75%;
- wrong lock: 25%;
- catastrophic ancestry collapse: 100%;
- dominant clone: 100%;
- mean minimum unique-parent fraction: 0.0316.

The selected `c_lambda=2` requires substantially more concentrated destination evidence than the P1F-B `c_lambda=0.5` threshold. Together with `beta=0`, this reduces premature local cloning while retaining reliable eventual mode selection in the controlled experiment.

## Important residual warning

The freeze criterion is passed because validation obtains 100% correct lock, 0% wrong lock, and 22.5% catastrophic-collapse incidence, below the predeclared 25% maximum. However, dominant-clone events remain frequent at 85% of validation runs. The selected literal transition is therefore still genealogically aggressive.

This distinction must remain explicit:

- **accepted for P1F-D comparison** means the parameter tuple is frozen without using the P1F-D evaluation data;
- it does **not** mean the literal ACO operator has satisfactory particle diversity for deployment;
- P1F-G remains necessary for a scalable and diversity-aware adaptation.

## Reproducibility

- configuration: `configs/phase1f_c.yaml`
- runner: `experiments/run_phase1f_c.py`
- method note: `docs/experiments/phase1f_c_aacopf_tuning.md`
- workflow run: `33803059007`
- artifact: `p1f-c-results`, ID `9912216810`

## Next step

Proceed to **P1F-D** using the frozen literal-small setting `alpha=0.5`, `beta=0`, `c_lambda=2`. The comparison must use matched trajectories, UWB noise, initial clouds, and particle budgets, and must retain ancestry/diversity diagnostics alongside accuracy and convergence metrics.
