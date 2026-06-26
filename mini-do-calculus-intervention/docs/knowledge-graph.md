# Knowledge Graph -- Do Calculus Intervention

## L1: Definitions

- do-calculus, do-operator
- Intervention: do(X=x), replaces f_X with constant
- Graph surgery (mutilation): remove incoming edges
- Identifiability: P(y|do(x)) expressible from observational data
- Structural Causal Model (SCM): M = ⟨U, V, F, P(u)⟩
- Counterfactual: Y_x(u) — potential outcome under intervention X=x

## L2: Core Concepts

- do(X=x): setting X by external action, breaking natural mechanism
- do-calculus: 3 rules transforming interventional to observational probability
- Graph surgery: remove all incoming edges to intervened variables
- Truncated factorization: P(v|do(x)) = Π_{i∉X} P(v_i|pa_i)
- d-separation: graphical criterion for conditional independence
- Back-door criterion: Z is admissible if it blocks all back-door paths
- Front-door criterion: M mediates all causal paths from X to Y

## L3: Math Structures

- Causal DAG: directed acyclic graph encoding causal assumptions
- Adjacency matrix: DCI_Graph with double weights
- Path types: directed, back-door, undirected paths
- Collider/v-structure: X → M ← Y pattern
- Moral graph: marry parents, drop directions
- Topological ordering: causal precedence

## L4: Fundamental Laws

- Pearl do-calculus (1995): Three rules for interventional probability
- Truncated factorization (g-formula): P(v|do(x)) product decomposition
- Back-door adjustment: P(y|do(x)) = Σ_z P(y|x,z)P(z)
- Front-door adjustment: P(y|do(x)) = Σ_m P(m|x) Σ_{x'} P(y|x',m)P(x')
- Shpitser-Pearl completeness (2006): do-calculus is complete for identifiability
- Counterfactual consistency (Robins, 1986): X=x ⇒ Y_x = Y
- PNS bounds: max(0,y1-y0) ≤ PNS ≤ min(y1,1-y0)

## L5: Algorithms/Methods

- do-calculus Rule 1-3 application with d-separation check
- Graph mutilation (surgery): remove incoming edges to X
- Back-door set finding (minimal admissible)
- Front-door mediator set finding
- Identifiability check via back/front-door + do-calculus
- Back-door adjustment (Monte Carlo)
- Front-door adjustment (Monte Carlo)
- Propensity score estimation
- IPW (Inverse Probability Weighting)
- Standardization (g-formula)
- Doubly-robust estimation
- Stratification
- Matching estimator (nearest neighbor)
- Bootstrap confidence intervals
- E-value sensitivity analysis
- Counterfactual probability computation (3-step: abduction, action, prediction)
- Twin network construction
- PN/PS/PNS computation
- Causal mediation (NIE via counterfactuals)
- G-computation (sequential back-door)
- Difference-in-Differences
- Regression Discontinuity Design

## L6: Canonical Problems

- Simpson's Paradox: Gender confounds treatment-outcome relation
- IV model: Instrument Z → X, U → X, U → Y, X → Y
- Mediation model: X → M → Y, X → Y (Baron-Kenny)
- M-bias: Collider bias (U1→X, U1→Z, U2→Y, U2→Z)
- Front-door model: U→X, X→M→Y, U→Y (smoking→tar→cancer)
- Counterfactual necessity: "Would the patient have recovered without treatment?"

## L7: Applications

- L7-1: Clinical Trial Policy Evaluation (302 lines)
  - Simpson's Paradox resolution via back-door adjustment
  - Multiple estimators compared (ATE, ATT, ATU, IPW, standardization,
    doubly-robust, matching, DiD, RDD, bootstrap CI)
  - E-value sensitivity analysis
  - Counterfactual necessity/sufficiency
  - Policy recommendation with causal justification
- L7-2: Causal Mediation for Health Policy (350 lines)
  - Education → Health Literacy → Health mediation analysis
  - NDE/NIE/CDE decomposition with proportion mediated
  - Path-specific counterfactual effects
  - Twin network counterfactual queries
  - Individual-level counterfactual
  - Front-door criterion for Smoking→Tar→Cancer
  - Counterfactual fairness assessment
  - Shapley-style feature importance explanation
  - Population attributable fraction

## L8: Advanced Topics

- Path-specific counterfactual effects
- Counterfactual fairness (Pearl, 2000)
- Counterfactual explanation (Shapley-style importance)
- Attributable fraction (population-level)
- Sequential back-door (g-computation for longitudinal data)
- Marginal structural models (Robins, Hernán, Brumback, 2000)
- Doubly-robust estimation
- Collider bias detection (M-bias, butterfly bias)
- Safe adjustment checking
- Disjunctive cause criterion (VanderWeele, 2019)
- Multiple back-door admissible sets

## L9: Research Frontiers

- Do-calculus for time series (documented; g-computation implemented)
- Causal reinforcement learning (documented)
- Transportability (documented)
- z-identifiability (documented)
- Causal discovery integration (partial: v-structure detection, skeleton)

