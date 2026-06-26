# Course Tree — Do-Calculus & Intervention

## Prerequisite Dependency Tree

```
Statistics & Probability (core)
├── Conditional Probability, Bayes Theorem
├── Expectation, Variance, Covariance
└── Law of Total Probability

    ↓

Causal Inference Foundations (Pearl, 1995–2009)
├── L1: Structural Causal Models (SCM)
│   └── U, V, F, P(u) — formal definition
├── L1: Causal Graphs (DAGs)
│   └── Nodes, directed edges, adjacency
├── L2: d-Separation Criterion (Pearl, 1988)
│   ├── Chain X→M→Y
│   ├── Fork X←M→Y
│   └── Collider X→M←Y
└── L2: Causal vs. Statistical Assoc.

    ↓

Intervention Calculus (this module)
├── L3: Graph Surgery / Mutilation
│   └── Remove incoming edges to intervened vars
├── L4: Truncated Factorization
│   └── P(v|do(x)) = Π_{i∉X} P(v_i | pa_i)
├── L4: do-Calculus Rules 1–3 (Pearl, 1995)
│   ├── Rule 1: Insertion/deletion of observations
│   ├── Rule 2: Action/observation exchange
│   └── Rule 3: Insertion/deletion of actions
├── L4: Back-Door Criterion (Pearl, 1993)
│   ├── Definition: Z blocks back-door paths, Z ∉ Desc(X)
│   ├── Back-Door Adjustment Formula
│   └── Minimal sufficient sets
├── L4: Front-Door Criterion (Pearl, 1993)
│   └── Mediator intercepts all directed paths
└── L4: Identifiability (Shpitser & Pearl, 2006)
    └── do-calculus completeness

    ↓

Causal Effect Estimation (this module)
├── L5: ATE / ATT / ATU / CATE
├── L5: Propensity Score Methods
│   ├── Inverse Probability Weighting (IPW)
│   ├── Standardization (g-formula)
│   ├── Doubly-Robust Estimation
│   └── Stratification
├── L5: Matching & Regression Adjustment
├── L5: Bootstrap Confidence Intervals
├── L5: Sensitivity Analysis (E-value)
├── L5: DiD and RDD
├── L5: G-Computation & Sequential Back-Door
└── L5: Marginal Structural Models

    ↓

Counterfactual Reasoning (this module)
├── L5: Three-Step Process
│   ├── Abduction: Update P(u | E=e)
│   ├── Action: do(X=x)
│   └── Prediction: Compute P(Y=y)
├── L5: Twin Network Method
├── L5: PN / PS / PNS (Probability bounds)
├── L5: Individual-Level Counterfactuals
├── L5: Causal Mediation via Counterfactuals
│   ├── NDE / NIE / CDE decomposition
│   └── Path-Specific Effects
└── L5: Attributable Fraction

    ↓

Advanced Topics (L8, this module)
├── Counterfactual Fairness (Pearl, 2000)
├── Counterfactual Explanation (Shapley-style)
├── Collider Bias & Safe Adjustment
├── Disjunctive Cause Criterion (VanderWeele, 2019)
└── Multiple Back-Door Admissible Sets

    ↓

Research Frontiers (L9)
├── Do-Calculus for Time Series
│   └── g-computation + sequential back-door (partial impl)
├── Causal Reinforcement Learning
├── Transportability (Bareinboim & Pearl)
└── z-Identifiability

## Downstream Dependencies

This module is a prerequisite for:
- mini-causal-discovery/ (PC algorithm, FCI, score-based methods)
- mini-causal-inference-applications/ (applied causal ML)
- mini-counterfactual-ml/ (counterfactual fairness in ML pipelines)
- mini-transportability/ (generalization across populations)

## Key Textbook Progression

1. Pearl (2009) Causality, Ch. 1–4, 7  →  Core SCM & do-calculus
2. Hernán & Robins (2020), Ch. 1–15   →  Estimation methods
3. Morgan & Winship (2015), Ch. 1–9    →  Counterfactual framework
4. Pearl & Mackenzie (2018) The Book of Why → Conceptual grounding
