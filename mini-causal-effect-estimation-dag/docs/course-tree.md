# Course Tree -- Causal Effect Estimation DAG

## Prerequisites
- Probability theory (conditional probability, expectation, variance)
- Statistics (regression, MLE, hypothesis testing)
- Linear algebra (Gaussian elimination, matrix operations)
- C programming (structs, pointers, dynamic memory)

## Dependency Tree
```
Probability & Statistics
  |
  +-- Causal DAGs (Pearl, 2009) [L1-L3]
       |
       +-- d-separation & back-door paths [L3]
       |    |
       |    +-- Back-door / Front-door criteria [L4]
       |    |    |
       |    |    +-- Causal identification [L4]
       |    |    |    |
       |    |    |    +-- Adjustment set selection [L5]
       |    |    |         |
       |    |    |         +-- Covariate balance checks [L7]
       |    |    |
       |    |    +-- Propensity score estimation [L5]
       |    |         |
       |    |         +-- IPW estimation [L5]
       |    |         +-- Stratification [L5]
       |    |         +-- PS Matching [L5]
       |    |
       |    +-- SCM & do-calculus [L3-L4]
       |         |
       |         +-- G-Computation [L5]
       |         +-- Doubly Robust [L5]
       |         +-- Bootstrap SE [L5]
       |
       +-- Mediation analysis [L8]
       +-- Sensitivity analysis [L8]
       +-- Time-varying confounding / MSM [L8]
       |
       +-- Double/Debiased ML [L9 - Research Frontier]
       +-- Causal representation learning [L9 - Research Frontier]
```

## Course Alignment
- Harvard STAT 234: Causal Inference (DAGs, back-door, IPW, G-comp)
- Stanford MS&E 327: Causal Inference (SCM, do-calculus)
- CMU 36-708: Statistical Methods (graphical models, d-separation)
- UC Berkeley STAT 260: Causal Models (matching, DR, sensitivity)
- UCLA STAT 232: Causal Models (mediation, potential outcomes)
- Columbia P8109: Causal Inference (identifiability, adjustment)
- MIT 6.438: Algorithms for Inference (graphical models)
- Cambridge Part III: Statistical Theory (semiparametric efficiency)
- ETH 401-3620-00L: Causal Representation Learning (L9 frontier)

## L9 Research Frontiers Coverage
The following L9 topics have formal definitions in `causal_effect_dag.lean`:
- Double/Debiased Machine Learning (Chernozhukov et al., 2018)
- Causal Representation Learning (Scholkopf et al., 2021)
