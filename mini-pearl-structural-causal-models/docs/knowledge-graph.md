# Knowledge Graph — Pearl Structural Causal Models

## L1: Definitions (Complete)

| # | Term | Notation / Type | Source |
|---|------|-----------------|--------|
| 1 | Structural Causal Model (SCM) | M = ⟨U, V, F, P(U)⟩ | Pearl 2000, Def. 7.1.1 |
| 2 | Variable | SCM_Variable, SCM_VarType | scm_core.h |
| 3 | Edge/Structural Coefficient | SCM_Edge (from, to, coef) | scm_core.h |
| 4 | Exogenous variable (U) | SCM_EXOGENOUS | scm_core.h |
| 5 | Endogenous variable (V) | SCM_ENDOGENOUS | scm_core.h |
| 6 | Directed Acyclic Graph (DAG) | G = (V, E), no cycles | scm_graph.h |
| 7 | Path | SCM_Path | scm_core.h |
| 8 | Variable Set | SCM_VarSet | scm_core.h |
| 9 | Collider | X→Z←Y with X⟂̸Y | scm_graph.h |
| 10 | Moral Graph | Married parents of common children | scm_graph.h |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Pearl's Ladder of Causation (3 rungs) | Association → Intervention → Counterfactual |
| 2 | Causal ordering / topological sort | `scm_topological_sort()` |
| 3 | d-separation | `scm_d_separated()`, `scm_d_separated_set()` |
| 4 | Markov blanket | `scm_markov_blanket()` |
| 5 | Parents / Children / Ancestors / Descendants | `scm_get_parents()`, etc. |
| 6 | Structural equations (linear & general) | `SCM_Equation`, `scm_linear_eq` |
| 7 | do-operator (graph surgery) | `scm_do_intervention()` |
| 8 | Intervention formula | Truncated factorization, `scm_do_multi_intervention()` |
| 9 | Causal effect (ATE, ATT, ATC) | `SCM_CausalEffect` |
| 10 | Conditional independence testing | `scm_partial_correlation_test()` |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | DAG topology (adjacency matrix) | `bool adjacency[n][n]` |
| 2 | Topological ordering (Kahn's algorithm) | `scm_topological_sort()` |
| 3 | Linear SEM: Y = β·PA + U | `scm_linear_eq()` |
| 4 | General SEM: Y = f(PA, U) | `SCM_Equation` function pointer |
| 5 | Path enumeration (DFS) | `scm_all_directed_paths()` |
| 6 | d-separation via moral graph + BFS | `scm_d_separated()` |
| 7 | Covariance / precision matrix | `scm_partial_correlation_test()` |
| 8 | Fisher z-transform | `scm_fisher_z()` |
| 9 | Binomial discretization / mutual information | `scm_mutual_information_test()` |
| 10 | Gaussian elimination for linear systems | `solve_linear_system()` |

## L4: Fundamental Laws (Complete)

| # | Theorem | Statement | Implementation |
|---|---------|-----------|---------------|
| 1 | Causal Markov Condition | X ⟂_d Y | Z ⇒ X ⟂ Y | Z | `causal_markov_condition` (Lean) |
| 2 | Truncated Factorization | P(v|do(x)) = ∏ P(v_i|pa_i) | `truncated_factorization_applies` (Lean) |
| 3 | do-calculus Rule 1 | Insertion/deletion of observations | `scm_do_calculus_rule1()` |
| 4 | do-calculus Rule 2 | Action/observation exchange | `scm_do_calculus_rule2()` |
| 5 | do-calculus Rule 3 | Insertion/deletion of actions | `scm_do_calculus_rule3()` |
| 6 | Back-door Adjustment | P(y|do(x)) = Σ_z P(y|x,z)P(z) | `scm_back_door()` |
| 7 | Front-door Adjustment | Identification via mediator M | `scm_front_door()` |
| 8 | Counterfactual Consistency | Y_x = Y when X = x | `counterfactual_consistency` (Lean) |
| 9 | do-calculus Completeness | 3 rules suffice for identifiability | `do_calculus_completeness_closed` (Lean) |
| 10 | Mediation Additivity | TE = NDE + NIE | `mediation_additivity_nat` (Lean) |
| 11 | Markov Equivalence | Same skeleton + v-structures | `scm_markov_equivalent()` |
| 12 | Pearl's Causal Hierarchy | L3 ⊃ L2 ⊃ L1 | `ladder_is_hierarchical` (Lean) |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Source | Implementation |
|---|-----------|--------|---------------|
| 1 | PC Algorithm (constraint-based) | Spirtes et al. 2000 | `scm_pc_algorithm()` |
| 2 | GES Algorithm (score-based) | Chickering 2002 | `scm_ges_algorithm()` |
| 3 | FCI Algorithm (latent variables) | Spirtes et al. 2000 | `scm_fci_algorithm()` |
| 4 | Hill-Climbing Search | — | `scm_hill_climbing_search()` |
| 5 | Tabu Search | Glover 1989 | `scm_tabu_search_causal()` |
| 6 | Bootstrap Stability | Friedman et al. 1999 | `scm_bootstrap_stability()` |
| 7 | BIC Scoring (Gaussian SEM) | Schwarz 1978 | `scm_bic_score()` |
| 8 | AIC Scoring | Akaike 1974 | `scm_aic_score()` |
| 9 | Meek Orientation Rules (R1-R3) | Meek 1995 | `scm_orient_rule1/2/3()` |
| 10 | Partial Correlation CI Test | Fisher 1924 | `scm_partial_correlation_test()` |
| 11 | Mutual Information CI Test | Cover & Thomas 2006 | `scm_mutual_information_test()` |
| 12 | 2SLS (Two-Stage Least Squares) | Angrist & Pischke 2009 | `scm_iv_two_stage_least_squares()` |
| 13 | Propensity Score Matching (NN) | Rosenbaum & Rubin 1983 | `scm_propensity_score_matching()` |
| 14 | Double/Debiased ML | Chernozhukov et al. 2018 | `scm_double_ml_ate()` |
| 15 | g-formula (g-computation) | Robins 1986 | `scm_g_formula()` |
| 16 | Doubly Robust Estimation | Scharfstein et al. 1999 | `scm_doubly_robust_estimator()` |
| 17 | Heckman Two-Step Correction | Heckman 1979 | `scm_selection_bias_correction()` |

## L6: Canonical Problems (Complete)

| # | Problem | Description | Implementation |
|---|---------|-------------|---------------|
| 1 | Simpson's Paradox | Aggregate reversal in subgroups | `scm_simpsons_paradox_detect()` |
| 2 | M-bias | Collider bias from conditioning | `scm_m_bias_test()` |
| 3 | Berkson's Paradox | Hospital selection collider bias | `scm_berkson_paradox_test()` |
| 4 | Smoking-Lung Cancer SCM | Front-door via Tar mediator | `smoking_cancer_scm` (Lean) |
| 5 | Complier ACE (LATE) | IV Wald with non-compliance | `scm_complier_ace()` |
| 6 | Back-door Identification | Covariate adjustment set | `scm_back_door()` |
| 7 | Front-door with Unobserved Confounding | Pearl 2000, Ch. 3 | `scm_front_door()` |
| 8 | Prob. of Necessity/Sufficiency | PNS bounds for causation | `scm_prob_necessity_sufficiency()` |

## L7: Applications (Complete — 8+)

| # | Application | Domain | Implementation |
|---|-------------|--------|---------------|
| 1 | Clinical Trial Analysis | Medicine | `scm_clinical_trial_analysis()` (ITT/PP/AT) |
| 2 | Instrumental Variables (IV) | Econometrics | `scm_iv_wald_estimator()`, 2SLS |
| 3 | Propensity Score Matching | Epidemiology | `scm_propensity_score_matching()` |
| 4 | Difference-in-Differences | Policy Evaluation | `scm_difference_in_differences()` |
| 5 | Regression Discontinuity Design | Economics | `scm_regression_discontinuity()` |
| 6 | Mediation Decomposition | Psychology | `scm_mediation_decomposition()` |
| 7 | Double/Debiased ML | Econometrics | `scm_double_ml_ate()` |
| 8 | Simpson's Paradox Resolution | Epidemiology | `scm_simpsons_paradox_resolve()` |

## L8: Advanced Topics (Complete — 6+)

| # | Topic | Description | Implementation |
|---|-------|-------------|---------------|
| 1 | Transportability | Effect transport between populations | `scm_transportability_test()` |
| 2 | Selection Bias (Heckman) | Sample selection correction | `scm_selection_bias_correction()` |
| 3 | Counterfactual Fairness | Algorithmic fairness | `scm_counterfactual_fairness()` |
| 4 | Sensitivity Analysis (E-value) | Unmeasured confounding bounds | `scm_sensitivity_analysis_evalue()` |
| 5 | Manski Bounds | Partial identification | `scm_manski_bounds()` |
| 6 | g-formula / MSM | Time-varying treatments | `scm_g_formula()`, `scm_marginal_structural_model()` |
| 7 | Doubly Robust Longitudinal | Combine g-formula + IPW | `scm_doubly_robust_estimator()` |
| 8 | Monotonicity Check | Simplifies counterfactual bounds | `scm_monotonicity_check()` |

## L9: Research Frontiers (Partial)

| # | Topic | Status | Reference |
|---|-------|--------|-----------|
| 1 | Causal Representation Learning | Documented | Schölkopf et al. 2021 |
| 2 | Neural SCMs | Documented | Xia et al. 2021 |
| 3 | Counterfactual Generative Models | Documented | Pawlowski et al. 2020 |
| 4 | Double ML with Neural Nets | Documented | Chernozhukov et al. 2018 |
| 5 | Causal Discovery from Heterogeneous Data | Documented | Peters et al. 2016 |
| 6 | Causal Reinforcement Learning | Documented | Bareinboim 2020 |

