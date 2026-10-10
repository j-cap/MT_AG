# P2F conditional accuracy-communication consolidation

88 unique designs, 1760 unique seed records, 200 matched endpoint checks.
No new simulation. Empirical Pareto labels minimize attempted cost and mean error within fixed sensor/Q/range/update conditions.
Bootstrap intervals are descriptive and unadjusted. Single-budget reliability evidence is not a communication curve.

| Sensor | Q | Range case | Update | Scope | Tested budgets | Mean 10 cm target | p95 10 cm target |
|---|---|---|---|---|---|---|---|
| bias_half | frozen_q | clean | naive | all_available_designs | 0,180,360,720,3600 | 360 cyclic | 360 cyclic |
| bias_x10 | frozen_q | clean | naive | all_available_designs | 0,180,360,720,3600 |   |   |
| combined_poor | frozen_q | clean | naive | all_available_designs | 0,180,360,720,3600 |   |   |
| combined_poor | white_rescaled_q | clean | naive | all_available_designs | 0,180,360,720,3600 |   |   |
| noise_half | frozen_q | clean | naive | all_available_designs | 0,180,360,720,3600 | 360 cyclic | 360 cyclic |
| noise_half | white_rescaled_q | clean | naive | all_available_designs | 0,180,360,720,3600 | 180 cyclic | 360 cyclic |
| noise_x4 | frozen_q | clean | naive | all_available_designs | 0,180,360,720,3600 |   |   |
| noise_x4 | white_rescaled_q | clean | naive | all_available_designs | 0,180,360,720,3600 |   |   |
| reference | frozen_q | clean | gated | all_available_designs | 360 | 360 information |   |
| reference | frozen_q | clean | naive | all_available_designs | 0,180,360,720,1800,3600 | 360 information | 720 random |
| reference | frozen_q | clean | naive | cyclic_reference | 0,180,360,720,3600 | 360 cyclic | 720 cyclic |
| reference | frozen_q | clean | naive | uniform_reference | 0,180,360,720,1800,3600 | 360 uniform | 720 uniform |
| reference | frozen_q | loss30 | naive | all_available_designs | 360 | 360 cyclic |   |
| reference | frozen_q | loss50 | naive | all_available_designs | 360 |   |   |
| reference | frozen_q | network_outage | gated | all_available_designs | 360 | 360 information |   |
| reference | frozen_q | network_outage | naive | all_available_designs | 360 | 360 information |   |
| reference | frozen_q | nlos | gated | all_available_designs | 360 | 360 cyclic |   |
| reference | frozen_q | nlos | naive | all_available_designs | 360 |   |   |
| reference | frozen_q | node_outage | naive | all_available_designs | 360 | 360 information |   |
| reference | frozen_q | outliers | gated | all_available_designs | 360 | 360 information |   |
| reference | frozen_q | outliers | naive | all_available_designs | 360 |   |   |
