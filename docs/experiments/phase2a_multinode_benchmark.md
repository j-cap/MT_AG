# Phase P2A — four-node trajectory-first cooperative-localization benchmark

## Purpose

P2A establishes the canonical multi-node simulation benchmark used by the communication-efficiency studies. The benchmark follows the revised fixed assumptions:

- every node starts from a known position and yaw;
- every node has an IMU and a UWB transceiver;
- all nodes are mobile and uncertain;
- UWB measurements are pairwise between mobile nodes rather than between one target and a perfectly known auxiliary.

The P2A question is deliberately simple:

> Does full-rate all-pairs UWB cooperation materially improve the localization of a four-node IMU fleet relative to IMU-only propagation, and does the resulting benchmark expose useful geometry for later communication scheduling?

P2A is not yet a scheduling experiment. It defines the lower-communication and upper-communication reference cases that P2B/P2C will interpolate between.

## Trajectory-first truth model

The benchmark defines four smooth global paths first. Position, velocity, acceleration, yaw and yaw rate are analytic truth quantities. Ideal body-frame IMU signals are then synthesized from the global acceleration and heading before sensor noise and bias are added.

This reverses the earlier Phase-1 debugging construction in which IMU inputs generated the trajectory. The new causal chain is:

global path -> truth kinematics -> ideal IMU -> noisy IMU -> estimator.

The four paths are visually distinct and intentionally generate time-varying pairwise geometry:

1. broad clockwise ellipse;
2. offset counter-clockwise ellipse;
3. smooth figure-eight;
4. offset smooth two-frequency loop.

The benchmark uses 60 s duration at 100 Hz IMU rate. The current path set stays below 1 m/s and avoids pair distances below about 1.5 m.

## Pairwise UWB graph

For four nodes the available undirected links are:

12, 13, 14, 23, 24, 34.

The maximum UWB opportunity rate is 10 Hz. The full-communication reference activates all six links at every opportunity. Over 60 s this gives 3600 pairwise range exchanges.

The normalized communication ratio is therefore 1.0 for P2A full communication and 0 for the IMU-only reference.

## Estimator

P2A uses a centralized joint EKF with state

X = [x_1^T, ..., x_4^T]^T

where every node has [p_x,p_y,v_x,v_y,psi].

Known initial position and yaw are initialized exactly. Initial velocity receives a small matched random error (std 0.05 m/s per component) so P2A does not silently assume the complete kinematic state is exact.

Each range measurement

z_ij = ||p_i-p_j|| + nu_ij

updates both endpoint states and creates cross-node covariance.

The same measured IMUs, initialization and range realization are used for the IMU-only and full-communication comparisons.

## Important infrastructure-free limitation

Pairwise ranges are invariant to a common translation of the entire fleet. With no fixed anchor, full UWB cooperation can strongly reduce differential/relative fleet error but cannot directly observe arbitrary common-mode global translation drift.

P2A therefore records three complementary error measures:

- absolute fleet position RMSE;
- fleet-centroid/common-mode position RMSE;
- centered relative-shape RMSE.

This decomposition is important for interpreting how much of the localization error communication can actually correct.

## Campaign

- 20 matched stochastic seeds;
- independent IMU noise/bias process per node;
- Gaussian UWB range noise sigma = 0.12 m;
- IMU-only reference;
- all-pairs 10 Hz UWB joint-EKF reference.

Primary outputs:

- fleet RMSE;
- worst-node RMSE;
- centroid RMSE;
- relative-shape RMSE;
- yaw RMSE;
- full communication count.

## Detailed report figures

P2A deliberately produces a richer figure set than a single result table:

1. four global truth trajectories + all six true pairwise distances;
2. speed, body-acceleration magnitude and yaw-rate histories;
3. representative per-node position-error histories for IMU-only and all-pair UWB;
4. per-node pointwise mean position error and sample-standard-deviation bands across 20 seeds;
5. aggregate fleet error history + paired seed-level RMSE comparison;
6. centroid/common-mode versus relative-shape error decomposition;
7. representative EKF position error versus covariance-based uncertainty scale.

These figures are the main evidence in the detailed report; the summary table is retained only for exact numerical reference.

## Handoff

If P2A confirms a meaningful gap between IMU-only and full communication, P2B will uniformly reduce the all-pairs UWB rate to trace the first communication--accuracy frontier. P2C will then hold communication count fixed and change which links are selected.


## Final P2A result

Across 20 matched seeds:

- IMU-only fleet RMSE: 2.322 +/- 0.605 m;
- all-pairs 10 Hz UWB fleet RMSE: 1.058 +/- 0.451 m;
- full communication is better in 20/20 seeds;
- mean worst-node RMSE: 3.238 m -> 1.071 m;
- centroid/common-translation RMSE: 1.0566 m -> 1.0564 m;
- relative-shape RMSE: 2.027 m -> 0.0456 m.

Thus full pairwise UWB reduces absolute fleet RMSE by about 54.5% and relative-shape RMSE by about 97.8%, while leaving the common translation mode essentially unchanged.

This is the central P2A structural result: pairwise infrastructure-free communication is highly effective at maintaining the relative fleet geometry, but it does not provide an absolute position reference for common-mode drift.

Canonical campaign:

- GitHub Actions run: 36020999903
- full communication count: 3600 pairwise exchanges per 60 s run
- compact evidence: results/phase2/p2a/
- detailed figures: report/figures/phase2/p2a_*.pdf

## P2A decision

Freeze the P2A trajectories, sensor model, seeds, joint EKF, and error decomposition for P2B. P2B should change only the uniform all-pairs UWB rate so that the first communication--accuracy frontier is not confounded by a new benchmark or estimator.
