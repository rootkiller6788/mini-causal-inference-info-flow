# Course Alignment: Counterfactual Reasoning

## Nine-School Curriculum Mapping

### MIT
- **6.438 Algorithms for Inference** (Shah)
  - SCM representation: `cfr_core.h/.c`
  - 3-step counterfactual algorithm: `cfr_counterfactual.c`
  - D-separation: `cfr_scm_d_separated`
- **6.7900 Machine Learning** 
  - Doubly robust estimation: `cfr_eff_ate_doubly_robust`
  - IPW: `cfr_eff_ate_ipw`

### Stanford
- **STATS 361 Causal Inference** (Imbens)
  - Potential outcomes: `cfr_potential.c`
  - Matching estimators: `cfr_eff_ate_matching`
  - IV estimation: `cfr_iv_estimate`
- **CS 229 Machine Learning**
  - Causal effect estimation: `cfr_effects.c`

### Berkeley
- **STAT 239A Causal Inference** (van der Laan)
  - SCM and do-calculus: `cfr_core.c`, Lean `do_calculus_rule*`
  - Backdoor adjustment: `cfr_scm_backdoor_adjustment`
  - TMLE precursor (G-computation): `cfr_g_computation`
- **CS 281A Statistical Learning Theory**
  - Mediation analysis: `cfr_mediation.c`

### CMU
- **36-702 Statistical Machine Learning** (Wasserman)
  - Double robustness proof: Lean `doubly_robust_bias_zero`
  - Propensity score methods: `cfr_po_propensity_logistic`
- **10-708 Probabilistic Graphical Models**
  - DAG analysis: `cfr_scm_is_dag`, topological order

### Princeton
- **ORF 526 Causal Inference** (Small)
  - Mediation: `cfr_mediation.c`
  - Sensitivity analysis: `cfr_eff_sensitivity_analysis`, `cfr_med_sensitivity_rho`
  - IV bounds: `cfr_bounds_balke_pearl_iv`

### Caltech
- **CMS 290 Causality** (Eberhardt)
  - Pearl SCM: `cfr_scm_*` API
  - Counterfactuals: `cfr_cf_*` API
  - Probability of causation: `cfr_bounds_*`

### Cambridge
- **Part III Causal Inference** (Dawid)
  - Potential outcomes vs. decision-theoretic: `cfr_potential.c`
  - Bounds and partial identification: `cfr_bounds_*`
  - Manski bounds: `cfr_po_manski_bounds`

### Oxford
- **Advanced Causal Inference** (Silva)
  - Do-calculus formalization: Lean `do_calculus_*`
  - Identifiability conditions: `cfr_cf_is_identifiable`
  - Front-door criterion: `cfr_scm_frontdoor_criterion`

### ETH Z?rich
- **401-3620 Causal Inference** (B?hlmann/Peters)
  - SCM framework: `cfr_core.c`
  - Counterfactual reasoning: `cfr_counterfactual.c`
  - Causal discovery: `cfr_scm_pc_algorithm_adjacency`

## Cross-Cutting Themes

| Theme | Schools | Module Coverage |
|-------|---------|-----------------|
| SCM vs Potential Outcomes | MIT, Caltech, ETH, Oxford | Both frameworks unified |
| Do-calculus | Berkeley, Oxford, Cambridge | Rules 1-3 formalized |
| Mediation | Princeton, CMU, Stanford | Barron-Kenny + Pearl formula |
| IV Methods | Princeton, Stanford, Berkeley | Wald + Balke-Pearl bounds |
| Double Robustness | MIT, CMU, Berkeley | DR estimator + bias proof |
| Sensitivity Analysis | Princeton, Cambridge, Oxford | E-value, Rosenbaum, rho |
| Bounds/Partial ID | Cambridge, Caltech, Princeton | Manski, Tian-Pearl, Balke-Pearl |
