# P1F-G summary — scalable bounded-candidate AACOPF

Final corrected campaign: exact paper-level PF increment model, workflow run `35068506066`, artifact `10435397349`.

Selected development setting:

- candidate_count: 8
- max_move_fraction: 0.1
- destination_capacity_fraction: 0.025
- frozen score parameters: `(alpha,beta,c_lambda)=(0.5,0,2)`

The selected setting met the predefined development criterion on 16 controlled runs: 16/16 correct locks, 0/16 wrong locks, no catastrophic ancestry collapse, and no dominant-clone events.

## Controlled held-out validation

| method | correct lock | wrong lock | collapse | dominant clone | runtime [s] |
|---|---:|---:|---:|---:|---:|
| pf | 1.000 | 0.000 | 0.000 | 0.000 | 0.0235 |
| literal | 0.950 | 0.000 | 0.100 | 0.750 | 0.9872 |
| candidate_only | 0.450 | 0.525 | 1.000 | 1.000 | 0.0404 |
| move_guard_only | 1.000 | 0.000 | 0.000 | 0.000 | 0.0513 |
| destination_guard_only | 1.000 | 0.000 | 0.000 | 0.000 | 0.0634 |
| guarded | 1.000 | 0.000 | 0.000 | 0.000 | 0.0561 |

Candidate truncation alone is unsafe: it gives 45% correct lock, 52.5% wrong lock, and catastrophic collapse / dominant cloning in every held-out run. Either diversity guard restores 40/40 correct locks in this controlled benchmark; the destination-capacity guard limits the minimum unique-parent fraction to 0.975.

## Random-annulus global validation

| Np | method | pose success | late pos RMSE [m] | late yaw RMSE [deg] | runtime [s] |
|---:|---|---:|---:|---:|---:|
| 1200 | pf | 0.050 | 14.812 | 46.687 | 0.255 |
| 1200 | guarded | 0.000 | 10.699 | 32.638 | 0.570 |
| 1200 | candidate_only | 0.000 | 36.096 | 111.593 | 0.501 |
| 5000 | pf | 0.100 | 10.275 | 28.331 | 0.802 |
| 5000 | guarded | 0.100 | 8.000 | 22.672 | 1.797 |
| 10000 | pf | 0.150 | 10.017 | 29.514 | 1.525 |
| 10000 | guarded | 0.050 | 7.043 | 19.501 | 3.439 |

Initial correct-region support rises from 2.75 particles on average at Np=1200, to 9.85 at Np=5000, and 19.0 at Np=10000.

The guarded adaptation therefore does **not** improve strict global-pose acquisition. It does, however, reduce mean late position/yaw error at every evaluated budget while preventing the catastrophic candidate-only collapse mechanism.

## Runtime scaling

Bounded empirical log-log runtime slope: **0.926**.

Literal empirical log-log runtime slope: **1.947**.

At Np=40000 the bounded transition has median single-transition runtime 0.0281 s and scores 319,964 candidate pairs, compared with 1.6e9 hypothetical dense pairs (about 0.020% of the dense count).

## Decision

P1F-G supports a **scalability** and **diversity-safety** claim, but not a strict global-acquisition improvement claim. The bounded/guarded adaptation removes the dense all-pairs bottleneck and prevents catastrophic cloning in the controlled validation. The remaining global-localization limitation is finite support and mode survival under random-annulus initialization.

Full per-run evidence remains in GitHub Actions artifact `10435397349`.