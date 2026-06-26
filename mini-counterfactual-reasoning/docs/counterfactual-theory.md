# Counterfactual Reasoning Theory

## Structural Causal Models (Pearl, 2009)

An SCM is a triple (U, V, F):
- U: exogenous (background) variables
- V: endogenous (observed) variables
- F: structural equations v_i = f_i(pa_i, u_i)

## The 3-Step Counterfactual Algorithm

1. **Abduction**: Infer P(U|E=e) from evidence
2. **Action**: Intervene do(X=x) by replacing f_X
3. **Prediction**: Compute Y_x in modified model

## Potential Outcomes (Rubin, 1974)

Each unit i has Y_i(0), Y_i(1). Observed Y = T*Y(1) + (1-T)*Y(0).
Key: We never observe both for the same unit.

## ATE Estimation Methods

- Difference in means (RCT)
- Regression adjustment
- IPW (Inverse Probability Weighting)
- Doubly Robust
- Stratification
- Matching

## Mediation Analysis

TE = NDE + NIE where:
- NDE = Y(1,M(0)) - Y(0,M(0))
- NIE = Y(1,M(1)) - Y(1,M(0))

## References
- Pearl (2009). Causality. Cambridge.
- Imbens & Rubin (2015). Causal Inference. Cambridge.
