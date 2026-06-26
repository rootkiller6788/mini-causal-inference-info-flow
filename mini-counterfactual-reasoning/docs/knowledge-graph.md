# Knowledge Graph: Counterfactual Reasoning (L1-L9)

## L1 ? Definitions

| # | Concept | C Implementation | Lean Formalization |
|---|---------|-----------------|-------------------|
| 1 | Structural Causal Model (SCM) | `CFRSCM` struct | `SCM` structure |
| 2 | Variable (endogenous/exogenous) | `CFRVariable` struct | `Variable` structure |
| 3 | Structural Equation | `CFRStructuralEquation` struct | `Equation` structure |
| 4 | Intervention do(X=x) | `CFRIntervention` struct | `Intervention` inductive |
| 5 | Counterfactual Y_x(u) | `CFRCounterfactual` struct | `CounterfactualQuery` / `CFResult` |
| 6 | Potential Outcomes Y(0), Y(1) | `CFRPotentialOutcomes` struct | `PotentialOutcome` structure |
| 7 | Evidence E=e | `CFREvidence` struct | Evidence as List (Nat x Nat) |
| 8 | ATE / ATT / ATU | `CFRCausalEffects` struct | `ate` function |
| 9 | D-separation | `cfr_scm_d_separated` | `d_separation_implies_independence` theorem |
| 10 | Collider | `cfr_scm_is_collider` | `collider_bias` theorem |

## L2 ? Core Concepts

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Abduction-Action-Prediction (3-step) | `cfr_cf_abduction`, `cfr_cf_action`, `cfr_cf_prediction` |
| 2 | SUTVA (no interference + consistency) | Documented in `cfr_potential.h`, `sutva_no_interference` theorem |
| 3 | Ignorability (unconfoundedness) | `ignorability_implies_ate_identified` theorem |
| 4 | Overlap (positivity) | `overlap_bounds` theorem, propensity clamping |
| 5 | Backdoor Criterion | `cfr_scm_satisfies_backdoor`, `cfr_scm_backdoor_adjustment` |
| 6 | Frontdoor Criterion | `cfr_scm_frontdoor_criterion`, frontdoor SCM factory |
| 7 | Do-calculus (Rules 1-3) | Lean formalization of all 3 rules |
| 8 | Identifiability | `cfr_cf_is_identifiable`, identifiability theorems |
| 9 | Monotonicity assumption | `cfr_check_monotonicity`, monotonicity bounds |
| 10 | Randomization | `randomized` flag in CFRBoundsData |

## L3 ? Mathematical Structures

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | SCM as triple (U, V, F) | `CFRSCM`: vars[] + equations[], n_endogenous, n_exogenous |
| 2 | DAG structure | Adjacency via parents/children relations |
| 3 | Moral graph | `cfr_scm_moral_graph` |
| 4 | Topological ordering | `cfr_scm_topological_order` (Kahn's algorithm) |
| 5 | Wright's path coefficients | `cfr_scm_path_coefficient` |
| 6 | Directed paths enumeration | `cfr_scm_all_paths` (DFS with backtracking) |
| 7 | PC algorithm skeleton | `cfr_scm_pc_algorithm_adjacency` |
| 8 | Information flow | `cfr_scm_information_flow` |
| 9 | SCM transformation (clone, marginalize) | `cfr_scm_clone`, `cfr_scm_marginalize` |
| 10 | Propensity score e(X) | `CFRPropensity`, logistic regression |

## L4 ? Fundamental Laws

| # | Theorem/Law | C Code | Lean Theorem |
|---|-------------|--------|-------------|
| 1 | Counterfactual 3-step consistency | `cfr_cf_compute_full` | `counterfactual_consistency` |
| 2 | Mediation decomposition TE = NDE + NIE | `cfr_cf_mediation_effects` | `mediation_additive_decomposition` |
| 3 | Baron-Kenny product NIE = a*b | `cfr_med_barron_kenny` | `barron_kenny_product` |
| 4 | Double robustness | `cfr_eff_ate_doubly_robust` | `doubly_robust_bias_zero` |
| 5 | Manski bounds sharpness | `cfr_po_manski_bounds` | `manski_bounds_sharp` |
| 6 | Balke-Pearl IV bounds | `cfr_bounds_balke_pearl_iv` | `balke_pearl_iv_bounded` |
| 7 | PN under monotonicity | `cfr_bounds_with_monotonicity` | `pn_monotonicity` |
| 8 | Excess fraction identification | `cfr_bounds_excess_fraction` | `excess_fraction_formula` |
| 9 | Counterfactual identifiability | `cfr_cf_is_identifiable` | `counterfactual_deterministic` |
| 10 | ATE oracle consistency | `cfr_eff_ate_oracle` | `ate_single` |

## L5 ? Algorithms/Methods

| # | Algorithm | Function | Complexity |
|---|-----------|----------|------------|
| 1 | IPW estimator | `cfr_eff_ate_ipw` | O(n) |
| 2 | Doubly Robust estimator | `cfr_eff_ate_doubly_robust` | O(n) |
| 3 | Regression adjustment | `cfr_eff_ate_regression` | O(n?p) |
| 4 | Nearest-neighbor matching | `cfr_eff_ate_matching`, `cfr_eff_genetic_matching` | O(n?) |
| 5 | Stratification | `cfr_eff_ate_stratified` | O(n?k) |
| 6 | Barron-Kenny mediation | `cfr_med_barron_kenny` | O(1) |
| 7 | Pearl mediation formula | `cfr_med_pearl_formula` | O(1) |
| 8 | IV Wald estimator | `cfr_iv_estimate` | O(n) |
| 9 | SCM topological sort | `cfr_scm_topological_order` | O(V+E) |
| 10 | Monte Carlo counterfactuals | `cfr_cf_monte_carlo` | O(samples?V) |

## L6 ? Canonical Problems

| # | Problem | Example/Demo |
|---|---------|-------------|
| 1 | ATE estimation from RCT | `example2_cfr.c`, `example_app_medical.c` |
| 2 | Mediation decomposition | `example3_cfr.c` (Barron-Kenny + Pearl) |
| 3 | PN/PS/PNS bounds computation | `example3_cfr.c` (experimental + monotonicity) |
| 4 | IV bounds for non-compliance | `cfr_bounds_balke_pearl_iv`, `example3_cfr.c` |
| 5 | Front-door adjustment | `cfr_create_frontdoor_scm`, `cfr_front_door_adjustment` |
| 6 | Collider bias detection | `cfr_scm_is_collider`, `example1_cfr.c` |

## L7 ? Applications

| # | Application | File |
|---|-------------|------|
| 1 | Medical RCT: Blood pressure drug effect | `example_app_medical.c` |
| 2 | Policy Evaluation: DiD for minimum wage | `example_app_policy.c` |
| 3 | Manski bounds for policy evaluation | `cfr_po_manski_bounds` in policy example |

## L8 ? Advanced Topics

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Twin network for counterfactuals | `cfr_cf_twin_network` |
| 2 | Transportability across populations | `cfr_cf_transportability` |
| 3 | Path-specific effects | `cfr_cf_path_specific_effect` |
| 4 | Sequential counterfactuals | `cfr_cf_sequential` |
| 5 | Counterfactual fairness | `cfr_cf_fairness_measure` |
| 6 | Controlled Direct Effect | `cfr_cf_controlled_direct_effect` |

## L9 ? Research Frontiers

| # | Topic | Status |
|---|-------|--------|
| 1 | Neural causal inference (deep IV) | Documented |
| 2 | Causal representation learning | Documented |
| 3 | Continuous-time causal models | Documented |
| 4 | Causal discovery with latent variables | Documented |

