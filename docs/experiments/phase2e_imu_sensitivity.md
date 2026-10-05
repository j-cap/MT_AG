# P2E: IMU sensitivity of the accuracy-communication curve

## Frozen question and protocol

How do sensor white noise and slowly varying bias change the UWB effort needed to maintain relative localization accuracy? The protocol is `configs/phase2e.yaml` and inherits the frozen P2A configuration. No stress amplitude, EKF adjustment, or accuracy target is selected using campaign outcomes.

The same four trajectories run for 60 s at 100 Hz. Initial position and heading are known; the same 0.05 m/s standard deviation initial velocity errors are retained even for improved sensors. UWB ranges are clean Gaussian with sigma 0.12 m. There is no gate, loss, outage, or external position reference. Cyclic scheduling uses the P2C allocation at 600 common 10 Hz opportunities.

Budgets are 0, 180, 360, 720, and 3600 attempted exchanges (0%, 5%, 10%, 20%, and 100% of the full reference). All attempts are received and accepted. The 3600 case activates all six links at every 10 Hz tick and reproduces P2A; the intermediate schedules have P2C timing, not the all-link rounds of P2B.

| Condition | White-noise amplitude multiplier | Bias random-walk amplitude multiplier |
|---|---:|---:|
| Reference | 1 | 1 |
| Noise half | 0.5 | 1 |
| Noise x4 | 4 | 1 |
| Bias half | 1 | 0.5 |
| Bias x10 | 1 | 10 |
| Combined poor | 4 | 10 |

White multipliers apply jointly to accelerometer and gyroscope per-sample standard deviations. Reference values are 0.02 m/s^2 and 0.001 rad/s. Bias multipliers apply jointly to their random-walk diffusion amplitudes, 5e-5 m/s^2/sqrt(s) and 5e-6 rad/s/sqrt(s). Bias starts at zero; the experiment does not introduce initial constant offsets. Multipliers change standard deviations/diffusion amplitudes, not variances. They are controlled sensitivity levels, not calibrated grades of commercial IMUs.

All 20 seeds reuse the same latent sensor realizations, range noise, and initial velocity errors. A reference IMU realization is decomposed into truth + random-walk bias + white noise, then the two error components are scaled separately. The reference condition retains the original measured array byte-for-byte. Acceleration and gyro quality are varied together; this study does not isolate the contribution of each sensor axis or type.

## Estimator assumptions

The primary sweep keeps P2A process standard deviations at 0.03 m/s^2 and 0.0015 rad/s. This measures sensitivity of the existing estimator under changed sensor quality, including any process-noise mismatch that creates.

A separately labeled companion uses the identical measurements and schedules for noise-half, noise-x4, and combined-poor, with both process standard deviations multiplied by the corresponding white-noise multiplier. This retains P2A's 1.5 ratio of process to sensor white-noise amplitudes. The adjustment scales the existing Q, including its dt powers; it does not change the mechanization or estimate bias states. Correlated random-walk bias remains absent from the five-state EKF model. Consequently `white_rescaled_q` is **not a fully calibrated or bias-matched estimator**, particularly in combined-poor.

Six conditions x five budgets x 20 seeds = 600 primary runs. Three companion conditions x five budgets x 20 seeds = 300 additional runs. Total: 900. Zero-UWB state estimates must be identical across the two Q treatments for each sensor condition, though covariances differ.

## Metrics and predefined target

Report whole-run fleet, worst-node, centroid, centered relative-shape, per-node, and yaw RMSE from native 100 Hz histories. Shape subtracts instantaneous mean fleet position error without rotation alignment. Fleet squared RMSE equals centroid squared RMSE plus shape squared RMSE within each run.

The descriptive target is **0.10 m mean whole-run shape RMSE**, fixed before the campaign. Report the smallest tested budget reaching it, with no interpolation or minimum-rate claim outside the tested grid. Also report the empirical across-seed 95th percentile and number of seeds at or below 0.10 m. The percentile from 20 fixed-path realizations is not a deployment guarantee or a confidence bound. If even 3600 does not meet the target, record it as not attained.

Time profiles sample at 2 Hz and show pointwise mean +/- sample SD (ddof=1), clipped at zero. Mean curves are averages of instantaneous RMS errors, not the curve whose time RMS equals the mean of whole-run RMSEs. Paired differences at equal budget and against each condition's full-budget result use 30,000 matched-seed bootstrap resamples (seed 20261005) and descriptive percentile 95% intervals, without multiple-comparison correction.

## Reproduction and verification

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python experiments/run_phase2e.py --workers 4
MPLCONFIGDIR=/tmp/mt_ag_mpl python experiments/generate_phase2e_figures.py
```

Checkpoints are atomic, ignored local intermediates, fingerprinted against both configuration and sensor/filter/metric implementation. Committed seed metrics, time profiles, and summary reproduce figures. All five reference-condition endpoints in all 20 seeds must reproduce P2C cyclic at 180/360/720 and P2B at zero/full to 1e-10 tolerance. Every run verifies exchange counts and the error decomposition. Protocol tests verify separate scaling of matched latent components and that changing Q changes covariance without changing zero-UWB state propagation.

## Results and decision

All 900 runs completed. All 100 reference endpoints matched inherited P2B/P2C metrics. All 60 zero-UWB companion trajectories exactly matched their fixed-Q counterparts. Ten Phase 2 test functions passed when invoked directly (pytest and Ruff were unavailable in the execution environment), including separate noise/bias scaling and covariance-only changes during dead reckoning. The checkpoint-resume path rebuilt identical summaries without rerunning seeds.

| Condition | Shape at 360, mean [m] | Shape at full 3600, mean [m] | Centroid at full, mean [m] |
|---|---:|---:|---:|
| Reference | 0.0824 | 0.0456 | 1.056 |
| White noise x0.5 | 0.0770 | 0.0404 | 1.054 |
| White noise x4 | 0.1569 | 0.1096 | 1.188 |
| Bias drift x0.5 | 0.0816 | 0.0446 | 1.057 |
| Bias drift x10 | 0.1562 | 0.1119 | 1.089 |
| Combined poor | 0.2070 | 0.1527 | 1.222 |

Under frozen Q, all three degraded cases were worse than reference in all 20 seeds at 360 exchanges. Fourfold noise and tenfold bias showed similar relative-shape degradation despite different disturbance structure. More ranges helped but did not recover the 0.10 m mean target even at full budget. This is a limit of the tested settings/estimator, not an information-theoretic limit or a demonstration that a bias-state filter cannot recover it.

Lower white noise improved shape by a mean 0.00544 m at 360, paired bootstrap CI [-0.00716,-0.00378], 18/20 wins. Halving the already small baseline bias changed it by only -0.00088 m, CI [-0.00156,-0.00020]. Improved sensor conditions retained the same initialization uncertainty.

Companion Q rescaling reduced noise-half shape at 360 from 0.0770 to 0.0729 m (19/20 wins) and combined-poor from 0.2070 to 0.1817 m (16/20 wins). For noise-x4, shape differences at 360 and full had intervals spanning zero. Combined-poor's full-budget difference also spanned zero. Q rescaling is not a universal cure and does not model correlated bias.

The target is reached by reference and improved conditions at 360 under frozen Q. Noise-half plus Q rescaling reaches the mean target at 180 (0.0989 +/- 0.0160 m) but only 11/20 seeds pass and the empirical p95 is 0.1234 m. Its p95 criterion needs 360, versus reference's 720. The half-count mean crossing is descriptive and close to the threshold; it is not a reliable 50% communication saving.

Shape and yaw can rank process settings differently. Noise-x4's full-budget mean yaw error decreases from about 0.163 to 0.119 degrees with Q rescaling while shape stays nearly unchanged. Noise-half's yaw error instead increases from about 0.045 to 0.065 degrees while shape improves. No consistency claim or yaw-tuned optimal setting is made.

Proceed to P2F consolidation. Keep sensor condition, process-covariance treatment, schedule timing, and P2D reliability explicit. Pairwise ranges do not introduce an absolute reference; centroid drift dominates fleet error. The fixed paths, 60 s duration, zero initial bias, jointly varied acceleration/gyro quality, and stylized multipliers limit generalization. No hardware energy or sensor-grade comparison is claimed.
