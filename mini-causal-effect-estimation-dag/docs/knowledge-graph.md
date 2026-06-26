# Knowledge Graph — Causal Effect Estimation Dag

## L1: Definitions (Complete — 11 items)

- Causal effect: P(Y|do(X)) — distribution under intervention
- Identifiability: effect expressible from observational data
- DAG: Directed Acyclic Graph (V,E) no directed cycles
- SCM: Structural Causal Model (U,V,F)
- d-Separation: graphical conditional independence criterion
- Backdoor criterion: Z blocks all backdoor X->Y paths
- Frontdoor criterion: M intercepts all X->Y paths
- do-calculus: three rules for transforming do-expressions
- Propensity score: e(X) = P(T=1|X), balancing score
- ATE/ACE: Average Treatment/Causal Effect
- Potential outcomes: Y(0), Y(1) — Neyman-Rubin framework
- NDE/NIE/CDE: Natural Direct/Indirect/Controlled Direct Effects
- E-value: minimum confounder strength for nullification
- Rosenbaum Gamma: odds ratio of hidden bias

## L2: Core Concepts (Complete — 8 items)

- P(Y|do(X)) differs structurally from P(Y|X)
- Backdoor: adjust for confounders blocking paths
- Frontdoor: use mediator for unmeasured confounders
- Adjustment set validity conditions
- d-separation as conditional independence in DAGs
- Double robustness: consistent if either model correct
- Stratification removes ~(100/K)% confounding bias
- Balancing property of propensity score

## L3: Math Structures (Complete — 9 items)

- ATE = E[Y|do(X=1)] - E[Y|do(X=0)]
- Backdoor formula: sum_z P(y|x,z) P(z)
- Frontdoor formula: sum_m P(m|x) sum_x' P(y|x',m) P(x')
- IPW: (1/n) sum [T_i Y_i/e_i - (1-T_i)Y_i/(1-e_i)]
- DR: combines IPW + outcome model with correction term
- Propensity: e(X) via logistic regression MLE
- Rosenbaum: p+ = Gamma/(1+Gamma), p- = 1/(1+Gamma)
- E-value: E = RR + sqrt(RR*(RR-1))
- NIE = a1 * b2 (product-of-coefficients)
- PM = NIE/TE (proportion mediated)

## L4: Fundamental Laws (Complete — 8 items)

- Pearl backdoor theorem (1993): conditional independence implies identifiability
- Pearl frontdoor theorem (1995): mediator-based identification
- Rosenbaum-Rubin propensity score theorem (1983): balancing property
- Do-calculus completeness (Shpitser & Pearl, 2006)
- Doubly robust consistency (Robins, Rotnitzky, Zhao, 1994)
- Cornfield conditions (1959): confounder strength requirements
- Total effect decomposition: TE = NDE + NIE
- Rosenbaum sensitivity bounds (2002): p-value interval under Gamma

## L5: Algorithms/Methods (Complete — 17 items)

- Kahn's Topological Sort (O(V+E))
- DFS Cycle Detection with 3-coloring (O(V+E))
- BFS Path Finding (O(V+E))
- Back-Door Path Enumeration (O(V * 2^V))
- d-Separation Test via moral graph (O(V * (V+E)))
- Back-Door Adjustment Subset Search (O(2^V))
- Logistic Regression for Propensity Score (gradient descent)
- IPW Estimator (Horvitz-Thompson, stabilized weights)
- G-Computation via OLS (normal equations + Gaussian elimination)
- Doubly Robust Estimator (combined IPW + outcome)
- Stratification on Propensity Score Quantiles (O(n log n))
- 1:1 Nearest Neighbor PS Matching with Caliper
- Bootstrap SE (Efron & Tibshirani, 1993)
- Baron-Kenny Mediation (product-of-coefficients + delta method)
- IKT Simulation-Based Mediation (quasi-Bayesian MC)
- Rosenbaum Bounds for Paired Studies (Wilcoxon signed-rank)
- E-value Computation (VanderWeele & Ding, 2017)

## L6: Canonical Problems (Complete — 5 items)

- Smoking-Lung Cancer: frontdoor via tar deposits (Pearl, 2009)
- Drug Treatment Effect: backdoor for age/severity (Hernan & Robins, 2020)
- Job Training Program: Detroit workforce evaluation (LaLonde, 1986)
- Workplace Wellness Mediation: NDE/NIE decomposition (Baron & Kenny, 1986)
- Teaching Method Sensitivity: Rosenbaum + E-value (Rosenbaum, 2002)

## L7: Applications (Complete — 4 items)

- Medical: Drug A vs Drug B treatment effect estimation
- Economics: Detroit job training program (AER, 1986 reference)
- Public Health: Workplace wellness mediation (mental health pathway)
- Education: Teaching method sensitivity to unmeasured confounding

## L8: Advanced Topics (Complete — 6 items)

- Sensitivity analysis: Rosenbaum bounds, E-value, Cornfield conditions
- Mediation analysis: Baron-Kenny, IKT simulation, NDE/NIE/CDE
- Multiple mediators: joint indirect effect computation
- Longitudinal mediation: time-varying treatment/mediator effects
- Binary mediator: logistic mediation model
- Mediation sensitivity: rho parameter for M-Y confounding

## L9: Research Frontiers (Partial — documented)

- Causal ML / Double ML (Chernozhukov et al., 2018)
- Causal representation learning (Scholkopf et al., 2021)
- Instrumental variable methods
- Causal discovery algorithms (PC, FCI)
- Time-series causal inference
- Meta-complexity of causal structure learning
