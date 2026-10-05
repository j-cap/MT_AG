# P2D reliability at 360 attempted exchanges

20 matched seeds; no retries or replacements. Counts separate reception from acceptance.

| Case | EKF | Schedule | Received | Accepted | Shape RMSE [m], mean +/- SD | Fleet [m] | Worst [m] |
|---|---|---|---:|---:|---:|---:|---:|
| clean | gated | cyclic | 360.0 | 359.2 | 0.0824 +/- 0.0157 | 1.061 | 1.085 |
| clean | gated | information | 360.0 | 359.1 | 0.0810 +/- 0.0164 | 1.060 | 1.081 |
| clean | gated | p2b_uniform | 360.0 | 359.0 | 0.0870 +/- 0.0160 | 1.061 | 1.085 |
| clean | naive | cyclic | 360.0 | 360.0 | 0.0824 +/- 0.0159 | 1.061 | 1.086 |
| clean | naive | information | 360.0 | 360.0 | 0.0789 +/- 0.0141 | 1.060 | 1.083 |
| clean | naive | p2b_uniform | 360.0 | 360.0 | 0.0858 +/- 0.0164 | 1.061 | 1.085 |
| loss30 | naive | cyclic | 251.5 | 251.5 | 0.0919 +/- 0.0179 | 1.062 | 1.087 |
| loss30 | naive | information | 252.4 | 252.4 | 0.0939 +/- 0.0195 | 1.061 | 1.087 |
| loss30 | naive | p2b_uniform | 253.0 | 253.0 | 0.1002 +/- 0.0240 | 1.063 | 1.092 |
| loss50 | naive | cyclic | 175.0 | 175.0 | 0.1100 +/- 0.0203 | 1.064 | 1.100 |
| loss50 | naive | information | 179.4 | 179.4 | 0.1065 +/- 0.0199 | 1.064 | 1.097 |
| loss50 | naive | p2b_uniform | 178.3 | 178.3 | 0.1164 +/- 0.0300 | 1.065 | 1.106 |
| network_outage | gated | cyclic | 300.0 | 299.3 | 0.0918 +/- 0.0155 | 1.061 | 1.089 |
| network_outage | gated | information | 300.0 | 299.3 | 0.0913 +/- 0.0158 | 1.062 | 1.086 |
| network_outage | gated | p2b_uniform | 300.0 | 299.5 | 0.1002 +/- 0.0212 | 1.062 | 1.091 |
| network_outage | naive | cyclic | 300.0 | 300.0 | 0.0917 +/- 0.0158 | 1.061 | 1.089 |
| network_outage | naive | information | 300.0 | 300.0 | 0.0898 +/- 0.0132 | 1.062 | 1.088 |
| network_outage | naive | p2b_uniform | 300.0 | 300.0 | 0.0993 +/- 0.0211 | 1.062 | 1.090 |
| nlos | gated | cyclic | 360.0 | 337.9 | 0.0956 +/- 0.0624 | 1.067 | 1.094 |
| nlos | gated | information | 360.0 | 293.9 | 0.1305 +/- 0.0444 | 1.069 | 1.122 |
| nlos | gated | p2b_uniform | 360.0 | 338.9 | 0.0993 +/- 0.0288 | 1.063 | 1.100 |
| nlos | naive | cyclic | 360.0 | 360.0 | 0.1690 +/- 0.0119 | 1.073 | 1.148 |
| nlos | naive | information | 360.0 | 360.0 | 0.1478 +/- 0.0127 | 1.069 | 1.122 |
| nlos | naive | p2b_uniform | 360.0 | 360.0 | 0.1718 +/- 0.0133 | 1.074 | 1.149 |
| node_outage | naive | cyclic | 330.0 | 330.0 | 0.0883 +/- 0.0156 | 1.062 | 1.090 |
| node_outage | naive | information | 305.1 | 305.1 | 0.0853 +/- 0.0131 | 1.061 | 1.085 |
| node_outage | naive | p2b_uniform | 330.0 | 330.0 | 0.0892 +/- 0.0151 | 1.062 | 1.087 |
| outliers | gated | cyclic | 360.0 | 342.8 | 0.0838 +/- 0.0164 | 1.061 | 1.086 |
| outliers | gated | information | 360.0 | 340.9 | 0.0826 +/- 0.0147 | 1.059 | 1.085 |
| outliers | gated | p2b_uniform | 360.0 | 342.2 | 0.0888 +/- 0.0157 | 1.061 | 1.086 |
| outliers | naive | cyclic | 360.0 | 360.0 | 0.2955 +/- 0.0664 | 1.102 | 1.191 |
| outliers | naive | information | 360.0 | 360.0 | 0.3024 +/- 0.0558 | 1.106 | 1.193 |
| outliers | naive | p2b_uniform | 360.0 | 360.0 | 0.2830 +/- 0.0505 | 1.100 | 1.190 |

Cyclic remains the practical reference. Missing ranges produce moderate shape degradation.
Gating reduces positive-outlier cyclic shape RMSE from 0.2955 m to 0.0838 m.
Under persistent bias, gated information repeatedly attempts the biased link and is worse than gated cyclic (0.1305 versus 0.0956 m).
Reliability-aware scheduling requires a separate experiment against gated cyclic at equal attempts.
Bootstrap intervals are descriptive, paired by seed, and unadjusted for multiple comparisons.
