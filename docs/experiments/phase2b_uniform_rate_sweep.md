# Phase P2B — uniform all-pairs UWB rate sweep

## Question and controlled comparison

How much UWB communication is needed to retain the P2A localization benefit when *every* pair ranges at the same reduced rate? P2B measures the absolute fleet error together with its centroid and relative-shape components. It is a uniform-rate reference for the later fixed-budget link-selection study (P2C), not yet a claim that any scheduling policy is optimal.

The 60 s, four-node trajectory-first truth, six available undirected links, 100 Hz IMU, 20 seeds, IMU noise and bias, 0.12 m UWB noise, known initial position and yaw, 0.05 m/s initial-velocity perturbation, and joint EKF settings all come directly from `configs/phase2a.yaml`. Per seed, every rate shares the **same** IMU observations, initial perturbation, and pre-generated noisy range values. Only the activation mask for those ranges changes. The 0 Hz and 10 Hz runs are checked against the committed P2A per-seed metrics.

## Range schedule and cost

`configs/phase2b.yaml` specifies rates of 0, 0.5, 1, 2, 5, and 10 Hz. At each positive-rate tick, all six pairs range. The first tick is one period after initialization; the last is at 60 s. All ticks coincide with the 10 Hz opportunities on the 0.01 s simulation grid. The cost counts **individual undirected pairwise exchanges**, not six-link rounds, packets, bytes, or electrical energy.

| Uniform rate | Six-link rounds | Pairwise exchanges | Ratio to P2A full rate |
|---:|---:|---:|---:|
| 0 Hz | 0 | 0 | 0% |
| 0.5 Hz | 30 | 180 | 5% |
| 1 Hz | 60 | 360 | 10% |
| 2 Hz | 120 | 720 | 20% |
| 5 Hz | 300 | 1800 | 50% |
| 10 Hz | 600 | 3600 | 100% |

Every rate runs for all 20 matched seeds. The fleet and worst-node position RMSE, centroid position RMSE, and centered relative-shape RMSE are calculated over the same entire 60 s trajectory, including time zero. `time_profile.csv` stores the pointwise across-seed mean, sample standard deviation, and 95th percentile of the instantaneous fleet RMS error, centroid error, and centered relative-shape RMS error at 10 Hz output resolution. These pointwise bands show between-seed variability, not uncertainty in the estimated mean.

The squared fleet RMSE decomposes **within each run** into squared centroid RMSE plus squared shape RMSE; means of the RMSE values across seeds do not generally obey that identity.

## Reproduction

From the repository root, with the project dependencies installed:

```bash
PYTHONPATH=src python experiments/run_phase2b.py
PYTHONPATH=src python experiments/generate_phase2b_figures.py
```

The first command writes `results/phase2/p2b/seed_metrics.csv`, `time_profile.csv`, `summary.json`, and `summary.md` and fails if the P2A endpoints no longer match. The second creates `report/figures/phase2/p2b_rate_sweep.pdf` and `p2b_time_profiles.pdf` from the compact results.

## Results and interpretation

The 0 Hz and 10 Hz seed-level endpoints reproduce P2A (the campaign verifies every stored fleet, worst-node, centroid, shape, yaw, and individual-node RMSE within numerical precision). All five positive rates improve fleet RMSE against the matched IMU-only seed in all 20 runs.

| Rate | Exchanges | Fleet RMSE mean ± sample SD | Worst-node mean | Centroid mean | Shape mean |
|---:|---:|---:|---:|---:|---:|
| 0 Hz | 0 | 2.322 ± 0.605 m | 3.238 m | 1.057 m | 2.027 m |
| 0.5 Hz | 180 | 1.064 ± 0.447 m | 1.097 m | 1.056 m | 0.114 m |
| 1 Hz | 360 | 1.061 ± 0.450 m | 1.085 m | 1.056 m | 0.086 m |
| 2 Hz | 720 | 1.059 ± 0.449 m | 1.078 m | 1.057 m | 0.067 m |
| 5 Hz | 1800 | 1.058 ± 0.451 m | 1.073 m | 1.056 m | 0.054 m |
| 10 Hz | 3600 | 1.058 ± 0.451 m | 1.071 m | 1.056 m | 0.046 m |

At 0.5 Hz, 5% of full communication yields 54.2% lower mean fleet RMSE than IMU only, compared with 54.5% at 10 Hz. The mean *paired* absolute-error difference between 0.5 and 10 Hz is 0.0066 ± 0.0088 m (sample SD across seeds); full-rate fleet RMSE is lower for 17 of the 20 matched seeds. The relative-shape difference is more visible: 0.114 m at 0.5 Hz versus 0.046 m at 10 Hz. Centroid RMSE remains about 1.056 m throughout, and worst-node RMSE improves gradually from 1.097 to 1.071 m. The small whole-fleet RMSE difference therefore reflects a common-translation error floor, rather than identical quality of relative localization.

The results give a **sampled uniform-rate curve**, not a globally optimal Pareto frontier. Rates below 0.5 Hz, adaptive time/link scheduling, packet loss, outliers, other trajectories and actual radio energy have not been evaluated here. In particular, a 0.5 Hz all-link round is six exchanges every 2 s, not one exchange every 2 s. P2C should compare link/time selection with this baseline at exactly matched individual-exchange counts.

Detailed analysis, the 20-seed error bars, and mean ± SD time histories are in `report/sections/03p_phase2b.tex`.
