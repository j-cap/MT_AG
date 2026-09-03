# Phase P1F-A — paper-state conventional PF baseline

## Purpose

P1F-A implements the conventional particle-filter baseline in the **paper-level state and input representation** audited in P1E. It deliberately does **not** implement the ant-colony transition. The purpose is to establish a clean reference against which P1F-B and later AACOPF variants can be compared.

The estimator state is

\[
\mathbf x_k=[x_k,y_k,\phi_k]^\top,
\]

with dead-reckoning increments

\[
\mathbf u_k=[\Delta L_k,\Delta\phi_k]^\top.
\]

The primary propagation convention follows Han et al. Equation (3): translate using the previous azimuth and then update azimuth. The post-turn convention described in Algorithm 1 is retained as a sensitivity case.

## Questions

1. Can the audited three-state model reconstruct its own deterministic truth exactly when initialized at the true state and supplied exact increments?
2. With unknown initial position and yaw, does a conventional PF recover the pose under the informative moving-auxiliary geometry when the DR increments are exact and only UWB range noise is present?
3. How sensitive is the result to the paper's propagation inconsistency and to random-versus-structured first-range initialization?
4. What conventional-PF diagnostics should be frozen now so that P1F-B changes only the resampling/particle-transition mechanism?

## Controlled model

P1F-A uses a 60 s, 0.1 s paper-state trajectory generated directly from deterministic displacement and yaw increments. This avoids contaminating the paper-level baseline with the five-state IMU mechanization. The increment profiles are chosen to remain close in scale to the Phase-1 curved trajectory while satisfying the paper kinematics exactly.

The globally known moving auxiliary uses the same navigation-frame trajectory already used in P1C/P1D. One UWB range is available every 0.1 s and is corrupted by zero-mean Gaussian noise with sigma 0.12 m.

## Initialization

The primary initialization is the P1E audited interpretation:

- target-position bearing sampled uniformly on [-pi, pi);
- radial offset sampled uniformly within +/- Delta d around the first measured range;
- yaw sampled uniformly on [-pi, pi);
- Delta d = 3 sigma_UWB = 0.36 m;
- uniform initial weights.

Because z_0 is used to construct the annulus, the first likelihood update is z_1 rather than reusing z_0.

A deterministic structured ring/yaw grid is retained as a control variant. It is not labelled as the Han et al. initialization because the paper states that particles are randomly generated.

## Conventional PF

The baseline uses a bootstrap proposal with deterministic DR propagation in the primary experiment. The Gaussian range likelihood is evaluated in log space and normalized over all particles. The posterior state estimate is computed **before** any resampling:

- position: weighted mean;
- yaw: weighted circular mean;
- initial-yaw estimate: initial-yaw lineage of the current MAP particle.

Systematic resampling is applied when N_eff < gamma N_p, with gamma = 0.5. After ordinary resampling, weights are reset uniformly and particle lineage is copied with the selected ancestor indices.

No undefined particle rejection and no auxiliary-node rejection are enabled.

## Primary campaign

The intended full campaign is:

1. exact-model reconstruction test;
2. 20 matched random-initialization seeds with informative moving auxiliary;
3. the same seeds using the Algorithm-1/post-turn propagation convention;
4. structured initialization control on matched seeds;
5. particle-count sensitivity sufficient to select a practical conventional-PF reference for P1F-B.

The principal convergence criterion remains the Phase-1 convention: position error < 1 m and absolute yaw error < 10 deg continuously to the end with at least 5 s remaining. Late-window metrics use the final 20 s.

## Claim boundary

P1F-A is not an AACOPF result and not an exact numerical reproduction of Han et al. It is a transparent implementation of the audited paper-state PF baseline under repository-defined, documented trajectories and noise. Its role is to make the subsequent ACO mechanism comparison attributable to particle management rather than to a change in state representation, propagation, initialization, or likelihood.
