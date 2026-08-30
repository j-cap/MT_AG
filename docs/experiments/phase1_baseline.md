# Phase 1A — IMU dead-reckoning and UWB/PF baseline

## Question
Can we build a minimal, auditable planar inertial-navigation baseline in which raw IMU-like acceleration and yaw-rate measurements generate dead-reckoning drift, and UWB range updates reduce the resulting position error?

## Important change from the first bootstrap
The original bootstrap injected noise directly into displacement and heading increments `(ΔL, Δψ)`. That matches the high-level interface used by Han et al., but it does not match the intended thesis hardware, where dead reckoning is to be built from IMU measurements. Phase 1A therefore now uses body-frame acceleration plus gyroscope yaw rate as the propagation input.

## State and sensor assumptions
The planar navigation state is `x = [p_x, p_y, v_x, v_y, ψ]^T`. The IMU signal is `u_m = [a_x^b, a_y^b, ω_z]^T`. We assume a level platform and treat the horizontal accelerometer channels as gravity-compensated. This is a deliberate Phase-1 simplification; a real 3-D IMU requires roll/pitch estimation and explicit gravity removal.

Accelerometer-only dead reckoning is not sufficient for general 2-D motion because body-frame acceleration must be rotated into the navigation frame. We therefore use both accelerometer and gyroscope channels.

## Discrete mechanization
With midpoint yaw `ψ_{k+1/2} = ψ_k + 0.5 ω_{z,k} Δt`, body acceleration is rotated as `a_k^n = R(ψ_{k+1/2}) a_k^b`. We then propagate

`p_{k+1} = p_k + v_k Δt + 0.5 a_k^n Δt²`,

`v_{k+1} = v_k + a_k^n Δt`,

`ψ_{k+1} = wrap(ψ_k + ω_{z,k} Δt)`.

The simulated IMU contains white accelerometer/gyro noise and random-walk bias. Exact values are stored in `configs/phase1.yaml`.

## UWB and particle filter
A moving auxiliary node has known position `a_k`. Its range is `z_k = ||p_k-a_k|| + v_k`, with Gaussian range noise. The bootstrap PF propagates each particle with the measured IMU plus process perturbations and weights it with the UWB range likelihood. The known-pose initialization remains a debugging baseline; global unknown-pose initialization is a later Phase-1 task.

## Verification run (seed 42)
Local validation after the IMU reformulation: 5 tests passed. Ideal IMU mechanization reconstructs the generated trajectory exactly by construction. Noisy IMU DR yields position RMSE 1.378 m and final error 2.680 m. PF+UWB yields position RMSE 0.231 m and final error 0.324 m. This corresponds to an 83.3% reduction in position RMSE for this single verification realization.

Heading RMSE is 0.110 deg for noisy DR and 0.354 deg for PF+UWB. This is intentionally not hidden: a scalar range does not directly measure heading, and the current PF is not tuned or augmented with explicit bias states. Phase 1A therefore supports a position-correction claim only.

## Observability diagnostic
For the full local state `[p_x,p_y,v_x,v_y,ψ]`, a finite-horizon linearized observability matrix gives singular values approximately `(931.97,232.84,49.46,7.05,4.12e-14)` for the selected stationary auxiliary and `(916.48,320.81,45.58,8.05,0.602)` for the selected moving auxiliary. The stationary case is numerically rank deficient; the moving case is full rank but not especially well conditioned. This is a local diagnostic, not a complete nonlinear proof.

## Interpretation
The reformulated baseline now matches the intended sensing chain: IMU -> inertial mechanization -> drifting DR, with UWB used as an external range correction. The next tasks are multi-seed verification, explicit bias/initial-state sensitivity, unknown-pose initialization, and only then the audited AACOPF reproduction.
