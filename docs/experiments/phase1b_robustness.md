# Experiment P1B — Robustness of the known-pose baseline

## Question

Is the known-pose IMU/UWB particle-filter result from P1A robust, or is the improvement caused by a favorable random seed or a narrowly tuned parameter configuration?

## Hypotheses

1. Across independent sensor/PF random seeds, PF+UWB should consistently reduce position RMSE relative to open-loop IMU dead reckoning.
2. The local PF should tolerate moderate initial position/velocity/yaw error, but sufficiently large initial uncertainty should expose the transition toward the global-localization problem addressed in P1C.
3. UWB should partially compensate increasing inertial noise, but performance should degrade when the IMU propagation becomes too inaccurate.
4. Persistent accelerometer/gyro bias should be more difficult than zero-mean noise because the current PF does not estimate bias states.
5. Increasing particle count should show diminishing returns, with runtime increasing roughly with particle count.
6. The P1A PF process-noise choice should lie inside a useful region rather than at an isolated optimum; process noise that is too small is expected to cause particle under-dispersion.

## Controlled variables

To isolate robustness of the P1A estimator, the following remain unchanged:

- the 60 s deterministic curved target trajectory;
- the moving auxiliary-node trajectory;
- one UWB range every 0.1 s (10 Hz);
- Gaussian UWB range noise with sigma = 0.12 m;
- synchronous measurement timing;
- planar, gravity-compensated IMU model.

These assumptions are intentionally not relaxed in P1B. Update-rate, packet-loss, NLOS, and asynchronous-communication questions belong to Phase 2 or later.

## Randomness and matched comparisons

Sensor randomness and PF randomness use separate generators. The same integer seed produces matched normalized sensor-noise realizations across parameter conditions where this is meaningful. Initial-state-error draws use a third generator. This reduces comparison variance and avoids changing multiple random effects when one parameter is varied.

## Experiment groups

### P1B.1 — Multi-seed baseline

Run the unchanged P1A configuration for 20 seeds. Report the distribution of DR and PF position RMSE, p95/final error, heading RMSE, and the fraction of seeds in which PF+UWB improves position RMSE.

### P1B.2 — Initial-state uncertainty

Draw an actual initial-estimate error from four finite uncertainty levels. Use the same scale as the PF prior around that erroneous estimate. The largest two levels are deliberately boundary/stress cases, not an attempted solution to global localization.

### P1B.3 — IMU quality

Scale accelerometer/gyro white-noise standard deviations and bias-random-walk intensities together by 0.5, 1, 2, and 4. The underlying normalized random draws are matched across scales.

### P1B.4 — Persistent-bias stress test

Add synthetic constant accelerometer and gyroscope offsets while leaving the PF state unchanged. The purpose is to determine whether unestimated bias is a potentially important failure mechanism. The chosen bias levels are not hardware-calibrated and must not be presented as representative of the final device.

### P1B.5 — Particle-count/runtime trade-off

Evaluate 500, 1000, 2500, 5000, and 10000 particles over 20 matched sensor realizations. Record both localization performance and wall-clock runtime on the current execution platform. Runtime is only a relative software-complexity indicator; embedded timing must be measured separately in Phase 4.

### P1B.6 — PF process-noise sensitivity

Jointly scale the acceleration and yaw-rate process perturbations by 0.25, 0.5, 1, 2, and 4. This tests whether P1A relies on narrow process-noise tuning.

## Result

The aggregate numbers and decisions are recorded in `results/phase1/p1b/summary.md` and `summary.json`. The main findings are:

- over 20 baseline seeds, mean position RMSE is 1.491 m for noisy DR and 0.375 m for PF+UWB; PF is better in 20/20 seeds;
- moderate initial errors remain recoverable, while the 1 m / 5 deg and 2 m / 10 deg uncertainty regimes show catastrophic failures in some seeds;
- PF performance degrades progressively with IMU quality, especially at the 4x noise/bias-random-walk scale;
- persistent unestimated biases create substantial residual error despite UWB correction;
- 5000 particles are a conservative operating point; 10000 gives almost no mean-RMSE improvement for roughly twice the runtime;
- process-noise scales 1x and 2x perform similarly, whereas too-small process noise causes occasional failures.

## Interpretation

P1B supports the P1A position-correction result as a robust property of the controlled known-pose problem rather than a single-seed artifact. It also identifies two boundaries that should not be hidden by additional tuning: large initial-state uncertainty changes the problem qualitatively and should be handled in P1C, while persistent IMU bias may require either calibration or state augmentation depending on the real hardware characteristics.

## Decision

P1B is sufficient to move to P1C after review. The selected Phase-1 defaults remain 5000 particles and the original P1A process-noise values. No UWB communication assumptions are changed yet.
