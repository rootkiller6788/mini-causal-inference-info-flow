# Causal Effect Estimation Methods

## ATE (Average Treatment Effect)
ATE = E[Y(1) - Y(0)]

## ATT (Average Treatment Effect on Treated)
ATT = E[Y(1) - Y(0) | T=1]

## Estimation Methods

### Oracle
Uses known Y(0), Y(1) for all units. Unavailable in practice.

### Difference in Means
Unbiased under randomization. ATE = mean(Y|T=1) - mean(Y|T=0).

### IPW
Weights units by inverse propensity score. Requires correct propensity model.

### Doubly Robust
Consistent if either propensity OR outcome model is correct.

### Matching
Pairs treated units with similar control units.

## Probabilities of Causation

- PN: Probability of Necessity
- PS: Probability of Sufficiency  
- PNS: Probability of Necessity and Sufficiency

Bounds tighten with additional assumptions:
Trivial < Experimental < Monotonicity < Exogeneity+Monotonicity

## References
- Imbens & Rubin (2015)
- Hernan & Robins (2020)
