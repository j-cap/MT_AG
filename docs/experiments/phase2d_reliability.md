# P2D: reliability at equal attempted communication

The frozen protocol is `configs/phase2d.yaml`, inheriting all P2A trajectory, IMU, Gaussian ranging, initialization, process model, and 20-seed settings. Configuration and gate were fixed before the campaign; none were tuned to outcomes.

## Protocol

Each run attempts 360 undirected pairwise exchanges, 10% of full P2A communication. Cyclic and information use the identical P2C time allocation; the information score is unchanged. P2B's 1 Hz all-six-link schedule is a contextual equal-count reference with different timing. Policies do not receive current reception outcomes or current ranges until after selecting the pair. No retries or free replacements occur.

Cases: clean; 30% and 50% independent packet loss; network outage on [20,30) s; node 1 incident-link outage on [20,30) s; independent 5% positive outliers uniform +1 to +3 m; persistent +0.75 m bias on link 1-2 on [20,40) s. Cases are separate, not combined. Gaussian range sigma remains 0.12 m.

The ordinary EKF runs in every case. A fixed gate accepts NIS <= 9 in clean, network outage, outliers, and persistent bias: 20 seeds x 3 schedules x (7 ordinary + 4 gated) = 660 runs. Rejection preserves state and covariance; the selector receives no rejection penalty or reliability estimate. In simultaneous rounds, gating uses the current posterior after earlier accepted pairs in sorted order.

Potential noise/failure fields exist for every pair and time regardless of selection. IMU, range, and initial-velocity streams match P2A/P2B/P2C. Loss RNG offset is 9501, outlier RNG offset 9601, added to base_seed + 10000*seed. A shared uniform field gives nested potential 30/50% failures. Outcomes are matched across policies at equal link/time; realized successful counts can differ because selections differ.

Attempted, received, accepted, and rejected counts are separate. Received = accepted + rejected. Lost attempts consume budget. These are exchange counts, not measurements of airtime, retry cost, or hardware energy.

## Metrics and interpretation

Whole-run RMSE uses 100 Hz histories. Time profiles sample 2 Hz and aggregate instantaneous errors across seeds with sample SD (ddof=1). Shape is position error after subtracting instantaneous fleet mean error; no rotational alignment is applied. Pointwise mean +/- SD bands are not confidence intervals.

Window metrics use [10,20), [20,30), [30,40), and [40,60] s. The CSV `recovery_shape_rmse_m` names the fixed 30-40 s post-outage window; **it is still inside the bias interval in the NLOS case**. Post-bias recovery uses `late_shape_rmse_m` (40-60 s). No threshold-based recovery time is claimed.

Paired differences use matched seeds, 30,000 bootstrap resamples, seed 20260419, and percentile 95% intervals. Intervals are descriptive and unadjusted for multiple comparisons. They concern fixed paths under stochastic realizations, not diverse real deployments.

## Reproduction

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python experiments/run_phase2d.py --workers 4
MPLCONFIGDIR=/tmp/mt_ag_mpl python experiments/generate_phase2d_figures.py
```

Seed checkpoints are atomic and fingerprint configuration and core code; they are local, ignored intermediates. Committed metrics, link counts, profiles, and summary regenerate figures. All 60 clean ordinary-EKF endpoints (20 seeds x 3 schedules) reproduce the committed P2B/P2C metrics within 1e-10 tolerance. Accounting invariants are checked for every run. Tests cover failure-blind eligibility, lost-attempt accounting, rejection preserving prediction, clean acceptance, and nonfinite observations.

## Main findings

- Cyclic shape RMSE: clean 0.0824 m; 30% loss 0.0919 m; 50% loss 0.1100 m; network outage 0.0917 m. Late post-outage error returns close to clean.
- Information versus cyclic has small, uncertain differences under missing data. Under node outage it loses 54.9 of 360 attempts on average, versus cyclic's 30: roughly 91.5% of 60 outage-period choices involve the unavailable node.
- Positive outliers raise ordinary cyclic shape RMSE to 0.2955 m. Gating lowers it to 0.0838 m (71.6% lower mean), with 342.8 accepted ranges and 20/20 paired improvements. Gated information gives 0.0826 m. Rejection has a much larger effect than specialized selection here.
- Ordinary information improves on cyclic with persistent bias: 0.1478 vs 0.1690 m; paired difference -0.02117 m, CI [-0.02632,-0.01577], 18/20 wins. It happens to select fewer biased ranges (12.9 vs 20 during the bias interval); it does not classify NLOS.
- With the gate, information is worse than cyclic: 0.1305 vs 0.0956 m; paired +0.03484 m, CI [0.01997,0.04804], 1/20 wins. Gated information allocates 59.7/120 bias-period attempts to the biased link (49.8%), versus cyclic's 20 (16.7%). Rejection leaves covariance uncertainty high, consistent with repeated unproductive selection. This is a mechanism supported by attempt logs, not a full causal ablation.
- Post-bias gated shape RMSE over 40-60 s is 0.1489 m for information vs 0.0810 m for cyclic. The fixed gate has uneven protection under sustained bias; it is not a general robust estimator guarantee.

Retain cyclic as practical reference and proceed to P2E IMU sensitivity. A reliability-aware score requires a separately frozen equal-attempt comparison against gated cyclic. One budget, fixed paths, stylized stress levels, and one gate do not establish hardware performance, energy savings, estimator consistency, or an NLOS detector.
