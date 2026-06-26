# Course Alignment — Pearl Structural Causal Models

## Nine-School Cross-Reference

| School | Course | Topic Coverage | Module Implementation |
|--------|--------|---------------|----------------------|
| **MIT** | 6.045/18.400 Automata, Computability | Causal graphs as computational models | SCM graph as automaton: `scm_sample()` executes causal "computation" |
| **MIT** | 6.S978 Causal Reasoning | Pearl's ladder, counterfactuals, mediation | `scm_counterfactual_compute()`, `scm_mediation_analysis()` |
| **Stanford** | CS228 Probabilistic Graphical Models | d-separation, Markov equiv, structure learning | `scm_d_separated()`, `scm_pc_algorithm()`, `scm_ges_algorithm()` |
| **Stanford** | CS361 Causal Inference | do-calculus, back-door, front-door, IV | `scm_back_door()`, `scm_front_door()`, `scm_find_instrument()` |
| **Berkeley** | CS294 Causal Inference | SCM semantics, identifiability, completeness | `scm_is_identifiable()`, `do_calculus_completeness_closed` (Lean) |
| **CMU** | 10-708 Probabilistic Graphical Models | Constraint/score-based structure learning | `scm_fci_algorithm()`, `scm_hill_climbing_search()`, `scm_bic_score()` |
| **CMU** | 36-708 Statistical Methods in AI | Bootstrap stability, model selection | `scm_bootstrap_stability()`, `scm_aic_score()` |
| **Princeton** | COS 597C Causal Inference | IV, mediation, counterfactuals, PNS | `scm_iv_two_stage_least_squares()`, `scm_prob_necessity_sufficiency()` |
| **Princeton** | COS 522 Computational Complexity | Causality + hardness, causal inference limits | Theoretical foundations in Lean formalization |
| **Caltech** | CS/CNS 155 Learning Systems | Transportability, selection bias, fairness | `scm_transportability_test()`, `scm_selection_bias_correction()`, `scm_counterfactual_fairness()` |
| **Cambridge** | Part III Statistical Learning | Double ML, sensitivity analysis, doubly robust | `scm_double_ml_ate()`, `scm_sensitivity_analysis_evalue()`, `scm_doubly_robust_estimator()` |
| **Oxford** | Advanced Causal Inference | Pearl's hierarchy, PNS, counterfactual fairness | `ladder_is_hierarchical` (Lean), `scm_counterfactual_fairness()` |
| **ETH** | 263-5300 Causality | Formal SCM semantics, identifiability, completeness | Full Lean 4 formalization, 20+ structures, 6 theorems |

## Topic-by-Course Matrix

| Topic | UCLA | Stanford | Berkeley | CMU | Princeton | Caltech | Cambridge | Oxford | ETH |
|-------|------|----------|----------|-----|-----------|---------|-----------|--------|-----|
| SCM definition | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| d-separation | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| do-calculus | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | ✓ | ✓ |
| Back-door | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Front-door | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ✓ |
| Counterfactuals | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | ✓ | ✓ |
| PNS bounds | ✓ | - | - | - | ✓ | - | - | ✓ | ✓ |
| Mediation | ✓ | - | ✓ | - | ✓ | - | ✓ | ✓ | ✓ |
| PC algorithm | - | ✓ | ✓ | ✓ | - | - | - | - | - |
| GES algorithm | - | - | ✓ | ✓ | - | - | - | - | - |
| FCI algorithm | - | ✓ | - | ✓ | - | - | - | - | - |
| IV / 2SLS | - | ✓ | ✓ | - | ✓ | - | - | ✓ | - |
| PSM | - | ✓ | ✓ | - | ✓ | ✓ | - | - | - |
| DiD | - | ✓ | ✓ | - | ✓ | - | - | - | - |
| RDD | - | - | ✓ | - | ✓ | - | - | - | - |
| Double ML | - | ✓ | - | - | - | - | ✓ | - | - |
| g-formula | - | - | ✓ | - | - | ✓ | - | - | - |
| Sensitivity (E-val) | - | - | ✓ | - | - | - | ✓ | - | - |
| Transportability | ✓ | - | - | - | - | ✓ | - | - | - |
| Fairness | - | ✓ | - | - | - | ✓ | - | ✓ | - |

## Key References
- Pearl, J. (2000). *Causality: Models, Reasoning, and Inference*. Cambridge.
- Pearl, J. (2009). Causal inference in statistics: An overview. *Statistics Surveys* 3, 96-146.
- Pearl, J., Glymour, M. & Jewell, N.P. (2016). *Causal Inference in Statistics: A Primer*. Wiley.
- Pearl, J. & Mackenzie, D. (2018). *The Book of Why*. Basic Books.
- Spirtes, P., Glymour, C. & Scheines, R. (2000). *Causation, Prediction, and Search*. MIT Press.
- Peters, J., Janzing, D. & Schölkopf, B. (2017). *Elements of Causal Inference*. MIT Press.
- Hernan, M.A. & Robins, J.M. (2020). *Causal Inference: What If*. CRC Press.
- Angrist, J.D. & Pischke, J.S. (2009). *Mostly Harmless Econometrics*. Princeton.
- Imbens, G.W. & Rubin, D.B. (2015). *Causal Inference*. Cambridge.
- Chickering, D.M. (2002). Optimal Structure Identification. *JMLR* 3, 507-554.
- Meek, C. (1995). Causal inference and causal explanation. *UAI-95*.
- Chernozhukov, V. et al. (2018). Double/debiased machine learning. *The Econometrics Journal* 21(1).
- VanderWeele, T.J. & Ding, P. (2017). Sensitivity analysis. *Annals IM* 167(4).
- Bareinboim, E. & Pearl, J. (2016). Causal inference and data-fusion. *PNAS* 113(27).
- Heckman, J.J. (1979). Sample selection bias. *Econometrica* 47(1).
- Schölkopf, B. et al. (2021). Toward causal representation learning. *Proc. IEEE* 109(5).
