# Phase P1E — Han et al. AACOPF reproduction audit

## Purpose

P1E does **not** implement AACOPF. It converts Han et al. (Sensors 2020, 20(2):467) into an explicit implementation contract for P1F and records where the paper is reproducible, ambiguous, internally inconsistent, or under-specified.

The audit deliberately separates three categories:

- **Source-specified:** directly stated by Han et al. and can be implemented without interpretation.
- **Source-ambiguous:** the paper gives incomplete or internally inconsistent instructions.
- **P1F interpretation:** the explicit choice that will be used in this repository. These choices are never presented as values or formulas reported by Han et al.

This distinction is necessary because P1D already showed that geometry-induced non-observability and particle-approximation failure are different mechanisms. AACOPF may help only in the latter case; no particle-management rule can recover information absent from the measurement geometry.

## 1. Reproduction target

### Paper-level state and input

Han et al. use the three-state dead-reckoning representation

\[
\mathbf x_t=[x_t,y_t,\phi_t]^\top
\]

with inputs

\[
\mathbf u_t=[\Delta L_t,\Delta\phi_t]^\top.
\]

The AACOPF core therefore operates on **displacement and azimuth increments**, not directly on raw accelerometer and gyroscope samples. The paper's vehicle experiment obtains these increments from odometer/INS processing, while the pedestrian experiment obtains them from a separate strapdown/ZUPT pipeline.

**P1F decision:** implement a separate paper-level estimator with `[x,y,phi]` and `Delta L, Delta phi`. Do not embed AACOPF directly into the thesis five-state IMU mechanization until the paper-level mechanism has been tested. The current `[p_x,p_y,v_x,v_y,psi]` implementation remains the hardware-oriented thesis baseline.

## 2. Preconditions that must remain explicit

The paper states four operating assumptions:

1. an auxiliary node has more accurate localization than the target;
2. the target may lack initial global position and azimuth but must provide DR increments;
3. the target can range to the auxiliary node(s);
4. communication is available for ranging and localization data exchange.

Therefore the paper is **not** a demonstration that several completely unreferenced mobile nodes obtain absolute global coordinates from mutual ranges alone. The auxiliary trajectory is treated as known in the navigation frame.

**P1F decision:** the paper-level reproduction will use a globally known auxiliary trajectory. Rank-deficient single-auxiliary geometries from P1D remain negative controls.

## 3. Kinematic propagation audit

### Source-specified equation

Equation (3) gives

\[
x_t=x_{t-1}+\Delta L_{t-1}\cos\phi_{t-1},
\]
\[
y_t=y_{t-1}+\Delta L_{t-1}\sin\phi_{t-1},
\]
\[
\phi_t=\phi_{t-1}+\Delta\phi_{t-1}.
\]

### Internal inconsistency

The explanatory text immediately before Equation (3), and Algorithm 1 lines 7--8, instead describe motion using the updated direction `phi(t-1) + Delta phi(t)`. Thus the paper contains two propagation conventions:

- **pre-turn convention:** translate using `phi_{t-1}` and then update azimuth, consistent with Equation (3) and the later observability Jacobian;
- **post-turn convention:** translate using `phi_{t-1}+Delta phi_t`, consistent with Algorithm 1.

### P1F interpretation

The **primary** reproduction uses the pre-turn convention because Equation (3) is the formal motion model and the observability derivation is based on it. P1F must include a propagation-convention sensitivity check using the post-turn variant so that this ambiguity is not hidden.

## 4. Initial-particle construction audit

### Source-specified information

The paper states that the first UWB measurement defines an annulus centered on the known auxiliary position with radius approximately

\[
d_{12}(1)\pm\Delta d,
\]

and that `N x M` particles represent possible initial positions and azimuths. It also states later that the particles are randomly generated. Initial weights are uniform.

### Missing information

The paper does not specify:

- the numerical value of `Delta d`;
- whether radius is uniform, Gaussian, or sampled by another rule;
- whether position bearing is random or deterministically spaced;
- whether initial azimuths are random or evenly spaced;
- the random seed;
- whether radial uncertainty corresponds to the quoted UWB accuracy or another tuning parameter.

The generic Step 1 statement that particles are sampled from a known prior `P(x_0)` is also awkward relative to the paper's claim that initial position and azimuth are unavailable.

### P1F interpretation

The paper-faithful default will use random first-range initialization:

\[
\theta_i\sim\mathcal U[-\pi,\pi),
\]
\[
r_i=z_0+\delta r_i,\qquad \delta r_i\sim\mathcal U[-\Delta d,\Delta d],
\]
\[
\phi_{0,i}\sim\mathcal U[-\pi,\pi).
\]

The explicit repository default will be

\[
\Delta d=3\sigma_{\mathrm{UWB}},
\]

which is a **repository choice**, not a value reported by Han et al. A structured angular-grid initialization will be retained only as a sensitivity/control variant because P1C used that representation while Han et al. explicitly describe random particle generation.

## 5. Importance sampling and process model

### Source-specified information

Step 2 writes a generic proposal density

\[
x_t^i\sim q(x_t^i\mid x_{t-1}^i,h(t))
\]

and says the propagated particle set is obtained using Equation (3). Equation (12) is the generic sequential-importance weight formula.

### Missing information

The paper does not define the proposal density `q`, a process-noise covariance, or noise distributions for `Delta L` and `Delta phi` inside the AACOPF core.

### P1F interpretation

Use a bootstrap proposal:

\[
q(x_t\mid x_{t-1},z_t)=p(x_t\mid x_{t-1}),
\]

so Equation (12) reduces to a likelihood update. The primary mechanism-validation experiment will first use deterministic DR increments so the ACO transition can be isolated. Controlled `Delta L` and `Delta phi` perturbations will then be added as explicitly configured repository noise models.

## 6. UWB likelihood and weight update audit

### Source-specified information

The UWB model is Euclidean range plus ranging error. Section 3 states that ranging observations are modelled with a Gaussian distribution for weight updating. The particle with predicted range closest to the measured range is intended to receive higher weight.

### Internal/typographical problems

Algorithm 1 line 12 is not a usable probability-density expression as printed. Its parentheses/signs are ambiguous, the measurement variance is absent, and a `sort(...)` operation appears inside the weight computation. Equation (13) also normalizes over `N` in print even though the stated particle population is `N x M`.

### P1F interpretation

Use the standard Gaussian range likelihood

\[
\ell_t^i=-\frac12\sum_j\left(\frac{\hat d_{ij,t}-z_{j,t}}{\sigma_j}\right)^2,
\]

followed by numerically stable log-weight normalization over **all** `N_p=N x M` particles. `sort` is not part of the probability update; particle identity and lineage must remain intact.

This is the most conservative probabilistic interpretation of the paper's stated Gaussian observation assumption, but it is not claimed to reproduce the malformed Algorithm 1 line 12 literally.

## 7. Undefined particle rejection rule

Step 4 states that a particle weight is set to zero when an expression of the form `x_t^i · x < 0` is satisfied. The symbol on the right-hand side is not defined well enough to determine whether this means a coordinate sign constraint, inner product, task-specific admissible region, or a typesetting error.

**P1F decision:** do not implement this rejection rule in the core reproduction. Any scenario-dependent state constraint must be introduced separately and named explicitly. This prevents an undefined condition from silently changing the posterior.

## 8. Adaptive auxiliary-node rejection audit

### Source-specified rule

For auxiliary node `j`, Han et al. define a ranging residual `rho_j` and state that the node is eliminated if both

\[
E[\rho_j]>\Omega_1
\]

and

\[
D[\rho_j]>\Omega_2
\]

hold.

### Missing information

The paper does not specify:

- `Omega_1` or `Omega_2`;
- the sample window over which mean and variance are computed;
- whether the mean should be signed or absolute;
- which position estimate is used to compute the predicted range;
- minimum sample count before gating;
- hysteresis/re-entry logic after a node is rejected.

### P1F interpretation

The **single-auxiliary AACOPF mechanism comparison will disable this gate**. With only one auxiliary it is not needed to test the ACO resampling mechanism, and rejecting the sole measurement would confound the comparison.

A later robustness sub-experiment may add a separately configured rolling residual gate. Its window length and thresholds must be treated as new tuning parameters and may not be attributed to Han et al.

## 9. Ant-colony transition audit

### Source-specified intent

The paper's motivation is to move lower-weight particles toward more promising particles while allowing stronger particles to remain, rather than using ordinary copy/delete resampling. Equation (14) defines a transition probability proportional to a weight-information term raised to `alpha` and a position-information term raised to `beta`. Algorithm 1 defines the position factor using inverse particle distance. If the maximum transition probability exceeds `lambda`, particle `i` is moved to particle `j`; otherwise it is retained. The paper states that `lambda` is selected empirically or by experience.

### Critical missing definitions

The following details are not specified:

- numerical values of `alpha`, `beta`, and `lambda`;
- whether the weight difference is `w_j-w_i`, `|w_j-w_i|`, a ratio, or another transformed difference;
- whether destinations with lower weight are allowed;
- treatment of negative weight differences for non-integer `alpha`;
- the exact summation/candidate set in Equation (14);
- whether self-transitions are included;
- treatment of zero particle distance in `1 / ||p_i-p_j||`;
- whether all particle moves are computed from the pre-transition cloud or sequentially modify the candidate cloud;
- tie handling;
- whether the transition is deterministic to `argmax P_ij` or sampled from `P_ij` (the text uses the maximum, so deterministic movement is the closer reading).

### P1F primary interpretation

For a source particle `i`, define the candidate set

\[
\mathcal C_i=\{j:\;w_j>w_i\}.
\]

For `j in C_i`, use

\[
\Delta w_{ij}=w_j-w_i,
\]

and

\[
\eta_{ij}=\frac{1}{\|p_i-p_j\|_2+\epsilon_d}.
\]

Then

\[
s_{ij}=(\Delta w_{ij}+\epsilon_w)^\alpha\eta_{ij}^\beta,
\qquad
P_{ij}=\frac{s_{ij}}{\sum_{l\in\mathcal C_i}s_{il}}.
\]

Let

\[
j^*=\arg\max_{j\in\mathcal C_i}P_{ij}.
\]

If `C_i` is empty or `P_{ij*} <= lambda`, particle `i` is retained. Otherwise it is replaced by the **pre-transition** state and lineage of `j*`. All moves are computed synchronously from an immutable copy of the pre-transition particle cloud.

This interpretation is chosen because it directly implements the verbal statement that lower-weight particles should move toward higher-weight particles and that shorter moves are preferred. It is not uniquely implied by Equation (14).

## 10. Weight reset and output-state ordering

### Source-specified but inconsistent statements

After the ant-colony movement, the paper states that all weights are reset to `1/(N x M)`. Step 6 then gives a weighted state estimate, while Equation (16) derives initial azimuth from the index of the largest-weight particle. If all weights have already been reset uniformly, there is no unique largest-weight particle. Algorithm 1 also contains an unclear state-estimation line with an unexplained summation limit of three.

### P1F interpretation

Use standard filtering order:

1. propagate particles;
2. compute and normalize measurement weights;
3. compute the current posterior estimate and MAP lineage **before** resampling/ACO;
4. apply the ACO particle transition to construct the next generation;
5. reset post-transition weights uniformly.

The current position estimate is the weighted posterior mean before ACO. The current azimuth estimate uses a weighted circular mean. The recovered **initial azimuth** is the initial-azimuth lineage attached to the pre-ACO MAP particle. This makes Equation (16) meaningful while respecting the paper's stated uniform reset after resampling.

A literal `post-ACO uniform mean` implementation may be retained as a diagnostic variant, but it is not the primary interpretation because it conflicts with the MAP-index logic.

## 11. Parameter values and computational implications

The paper evaluates equal values of `N` and `M` from 50 up to 1500 and recommends `N=M=400` as a practical accuracy choice for its experiment. This corresponds to

\[
N_p=N M=160000
\]

particles. The paper also explicitly notes non-monotonic results caused by random particle generation and UWB error.

However, Equation (14) is written as if every source particle may consider every destination particle. A literal all-pairs calculation with 160000 particles requires on the order of

\[
N_p^2\approx2.56\times10^{10}
\]

pair scores per update, before other filtering work. The paper does not explain how this was implemented or accelerated and does not report AACOPF runtime or hardware for the filter computation.

**P1F decision:** maintain two clearly labelled implementations:

- **literal-small AACOPF:** full candidate evaluation for small particle populations, used to verify the audited transition rule;
- **scalable AACOPF adaptation:** bounded candidate pool for P1C-scale experiments, reported as an adaptation rather than a literal reproduction.

The scalable candidate policy must be defined and benchmarked in P1F before any embedded-feasibility claim.

## 12. Tuning policy for alpha, beta, and lambda

No numerical values for `alpha`, `beta`, or `lambda` are reported by Han et al. Therefore P1F must not select one arbitrary tuple and call it the paper configuration.

The initial neutral reference will be

- `alpha = 1`;
- `beta = 1`.

`lambda` will be parameterized relative to the uniform candidate probability. For a candidate count `K_i`, a useful dimensionless threshold factor is

\[
c_\lambda=\lambda K_i.
\]

P1F will sweep a compact set of `alpha`, `beta`, and `c_lambda` values on development seeds, then freeze the selected setting before evaluating held-out matched seeds. The final report must show the sensitivity rather than hide it.

This tuning protocol is a repository design decision motivated by the paper's under-specification.

## 13. Paper experiment facts that are reproducible

The paper reports the following concrete experimental facts that can be used as qualitative reference points:

- car-to-trolley experiment with the car as target and GPS-equipped trolley as known auxiliary;
- UWB module reported around 10 cm chip-level measurement accuracy and 20--30 cm bidirectional ranging accuracy;
- GPS reference around 5 cm in the vehicle experiment;
- a stationary trolley fails to recover the car trajectory/initial azimuth, while an appropriately moving trolley converges;
- true vehicle initial azimuth reported as 71.0548 deg and moving-case estimate 71.0009 deg;
- moving-case tracking accuracy described as up to about 0.3 m and initial-azimuth error 0.0539 deg;
- Table 2 reports particle-count sensitivity and selects `N=M=400`;
- pedestrian experiment uses two GPS-equipped auxiliary nodes, one with peak-shift corruption and the other with non-Gaussian/excess error, and reports AACOPF localization/azimuth accuracy of 0.66372 m and 0.22658 deg for the shown run.

These values **cannot be exact regression targets** because the paper does not provide the raw trajectories, ranging series, random seeds, AACOPF tuning parameters, measurement update frequency, or complete noise distributions.

## 14. What P1F is allowed to claim

P1F may claim a successful reproduction only at the level supported by the source:

1. the three-state DR + UWB problem is implemented as described;
2. a moving known auxiliary can resolve unknown global position/yaw while the selected rank-deficient controls cannot;
3. the audited ACO transition can be compared against conventional PF on identical simulated data;
4. under specified synthetic outlier/non-Gaussian conditions, AACOPF either does or does not improve convergence, mode preservation, outlier behavior, or particle diversity.

P1F must **not** claim exact reproduction of the paper's numerical errors unless the missing data and tuning information become available.

## 15. Required P1F diagnostics

Because Han et al. motivate AACOPF through particle starvation/diversity, P1F should measure more than position RMSE. At minimum record:

- position and azimuth RMSE;
- terminal convergence fraction and convergence time;
- probability/mass retained near the correct global mode;
- effective sample size before ACO;
- unique-particle fraction after ACO;
- spatial/yaw particle spread;
- fraction of particles moved by ACO;
- destination reuse / multiplicity;
- transition-score concentration or entropy;
- runtime per update and total runtime;
- sensitivity to `alpha`, `beta`, `lambda`, particle count, and initialization seed.

These diagnostics test the mechanism the paper actually claims: improved particle management, not only final tracking error.

## 16. P1F experiment sequence fixed by this audit

P1F should proceed in the following order:

1. **P1F-A — paper-state PF baseline:** implement `[x,y,phi]`, `Delta L`, `Delta phi`, random annular initialization, Gaussian UWB likelihood, no ACO.
2. **P1F-B — literal-small ACO mechanism:** implement the audited all-pairs transition on a small population and verify deterministic unit cases.
3. **P1F-C — parameter sensitivity:** sweep `alpha`, `beta`, and normalized `lambda` on development seeds.
4. **P1F-D — matched PF vs AACOPF:** compare on informative moving geometry using identical truth, DR, UWB, and initialization seeds.
5. **P1F-E — negative geometry controls:** confirm neither PF nor AACOPF can overcome stationary/constant-bearing non-observability.
6. **P1F-F — non-Gaussian/outlier stress:** evaluate the condition most directly connected to Han et al.'s claimed robustness.
7. **P1F-G — scalable transition:** only after the literal-small mechanism is understood, introduce and benchmark a bounded-candidate adaptation for larger particle populations and the later IMU-driven thesis model.

## Decision

P1E concludes that Han et al. provide enough information to reproduce the **problem structure and AACOPF concept**, but not enough information for a unique source-exact implementation. The largest ambiguities concern propagation timing, first-range sampling, the Gaussian weight formula, `alpha/beta/lambda`, the exact weight-difference definition, candidate selection, the undefined particle rejection condition, adaptive-node thresholds, and output ordering after uniform weight reset.

P1E is complete when these choices are committed and treated as the implementation contract for P1F. P1F must preserve a strict distinction between source-specified behavior and repository interpretations.