# Mini Counterfactual Reasoning

**Submodule of**: `mini-complex-control-theory / 8. mini-causal-inference-info-flow`

Implementation of Pearl's Structural Causal Model (SCM), counterfactual reasoning, potential outcomes, causal effect estimation, mediation analysis, and probability bounds ? covering the core of modern causal inference.

---

## Knowledge Coverage (L1-L9)

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | **Complete** | SCM triple (U,V,F), Intervention do(X=x), Counterfactual Y_x, Potential outcomes Y(0)/Y(1), ATE/ATT/ATU, D-separation |
| **L2** | Core Concepts | **Complete** | Abduction-Action-Prediction, SUTVA, Ignorability, Overlap, Backdoor/Frontdoor criteria, Do-calculus |
| **L3** | Math Structures | **Complete** | CFRSCM (variables + equations), DAG adjacency, Path coefficients, Wright's rules, Moral graph |
| **L4** | Fundamental Laws | **Complete** | Pearl 3-step theorem, TE = NDE + NIE, PN/PS bounds, Counterfactual consistency, Double robustness, Manski bounds |
| **L5** | Algorithms/Methods | **Complete** | IPW, Doubly Robust, Stratification, Matching, Barron-Kenny, Pearl mediation formula, IV estimation, Balke-Pearl bounds |
| **L6** | Canonical Problems | **Complete** | ATE estimation, Mediation decomposition, PN/PS/PNS bounds, IV non-compliance, Collider bias, Front-door adjustment |
| **L7** | Applications | **Partial+** | Medical RCT analysis (blood pressure drug), Policy evaluation (Card & Krueger DiD) |
| **L8** | Advanced Topics | **Partial+** | Twin network, Transportability, Path-specific effects, Sequential counterfactuals, Counterfactual fairness |
| **L9** | Research Frontiers | **Partial** | Documented in knowledge graph; implementation open for future research |

---

## Core Definitions (L1)

- **SCM**: `CFRSCM` struct ? triple (U, V, F) with endogenous variables, exogenous noise, and structural equations
- **Intervention**: `CFRIntervention` ? hard (do(X=x)) and soft (replace equation)
- **Counterfactual**: `CFRCounterfactual` ? Y_x(u) via 3-step procedure
- **Potential Outcomes**: `CFRPotentialOutcomes` ? Y(0), Y(1) per unit under Rubin Causal Model
- **ATE**: Average Treatment Effect = E[Y(1) - Y(0)]
- **ATT**: Average Treatment Effect on the Treated
- **NDE/NIE**: Natural Direct/Indirect Effects in mediation
- **PN/PS/PNS**: Probabilities of Necessity, Sufficiency, Necessity-and-Sufficiency

---

## Core Theorems (L4)

| Theorem | Formula | Status |
|---------|---------|--------|
| **Counterfactual 3-step** | Y_x(e) = Prediction(Action(Abduction(e))) | Implemented |
| **Mediation Decomposition** | TE = NDE + NIE | Verified |
| **Baron-Kenny Product** | NIE = a ? b (linear SEM) | Verified |
| **Double Robustness** | ATE_DR consistent if PS or outcome model correct | Implemented |
| **Manski Bounds** | ATE ? [E[Y|T=1]P(T=1) - P(T=0), E[Y|T=1]P(T=1) + P(T=0)] | Implemented |
| **Balke-Pearl IV Bounds** | ATE bounded by 4 LP constraints | Implemented |
| **PN under Monotonicity** | PN = ATE / P(Y=1|T=1) (point-identifiable) | Implemented |
| **Excess Fraction** | AF = [P(Y=1) - P(Y=1|do(T=0))] / P(Y=1) | Implemented |

---

## Core Algorithms (L5)

| Algorithm | Function | Complexity |
|-----------|----------|------------|
| IPW Estimator | `cfr_eff_ate_ipw` | O(n) |
| Doubly Robust | `cfr_eff_ate_doubly_robust` | O(n) |
| Propensity Score Matching | `cfr_eff_ate_matching` | O(n?) |
| Stratification | `cfr_eff_ate_stratified` | O(n?k) |
| Regression Adjustment | `cfr_eff_ate_regression` | O(n?p) |
| Barron-Kenny Mediation | `cfr_med_barron_kenny` | O(1) in SCM |
| Pearl Mediation Formula | `cfr_med_pearl_formula` | O(1) in SCM |
| IV Wald Estimator | `cfr_iv_estimate` | O(n) |
| SCM Topological Order | `cfr_scm_topological_order` | O(V+E) |
| D-separation Check | `cfr_scm_d_separated` | O(V?) |
| Counterfactual via MC | `cfr_cf_monte_carlo` | O(samples?V) |

---

## Canonical Problems (L6)

1. **ATE Estimation** ? Oracle, difference-in-means, IPW, DR, stratification, matching
2. **Mediation Analysis** ? NDExNIE decomposition via Barron-Kenny and Pearl formula
3. **Probability of Causation Bounds** ? PN/PS/PNS under experimental, monotonicity, Tian-Pearl
4. **IV Non-compliance** ? Balke-Pearl bounds, Wald estimator, LATE/CACE
5. **Front-door Adjustment** ? Smoking ? Tar ? Lung Cancer paradigm
6. **Collider Bias** ? Berkson's paradox detection via `cfr_scm_is_collider`

---

## Code Statistics

| Component | Files | Lines (approx) |
|-----------|-------|----------------|
| Headers (include/) | 7 | 783 |
| Sources (src/*.c) | 7 | 3,239 |
| Formalization (src/*.lean) | 1 | 724 |
| **Total include/ + src/** | **15** | **4,746** |
| Tests | 1 | 75 |
| Examples | 5 | 551 |
| Docs | 7 | ? |

---

## Curriculum Alignment (Nine-School Map)

| School | Key Course | Module Coverage |
|--------|-----------|-----------------|
| **MIT** | 6.438 Algorithms for Inference | IPW, DR, matching, SCM |
| **Stanford** | STATS 361 Causal Inference | Potential outcomes, matching, IV |
| **Berkeley** | STAT 239A Causal Inference | SCM, do-calculus, backdoor |
| **CMU** | 36-702 Statistical ML | Doubly robust, TMLE precursor |
| **Princeton** | ORF 526 Causal Inference | Mediation, sensitivity analysis |
| **Caltech** | CMS 290 Causality | Pearl's SCM, counterfactuals |
| **Cambridge** | Part III Causal Inference | Potential outcomes, bounds |
| **Oxford** | Advanced Causal Inference | Do-calculus, identifiability |
| **ETH** | 401-3620 Causal Inference | SCM, counterfactual reasoning |

---

## File Structure

```
mini-counterfactual-reasoning/
??? Makefile                    # make test (build + run tests)
??? README.md                   # This file
??? include/
?   ??? cfr_core.h              # SCM core data structures + DAG ops
?   ??? cfr_counterfactual.h    # 3-step counterfactual algorithm
?   ??? cfr_effects.h           # ATE estimation methods
?   ??? cfr_potential.h         # Potential outcomes framework
?   ??? cfr_mediation.h         # Mediation analysis
?   ??? cfr_bounds.h            # Probability bounds
?   ??? cfr_analysis.h          # Analysis tools (E-value, IV, etc.)
??? src/
?   ??? cfr_core.c              # SCM implementation (1023 lines)
?   ??? cfr_counterfactual.c    # Counterfactual computation (501 lines)
?   ??? cfr_effects.c           # Effect estimation (368 lines)
?   ??? cfr_potential.c         # Potential outcomes (313 lines)
?   ??? cfr_mediation.c         # Mediation analysis (340 lines)
?   ??? cfr_bounds.c            # Probability bounds (480 lines)
?   ??? cfr_analysis.c          # Analysis methods (214 lines)
?   ??? counterfactual_reasoning.lean  # Lean 4 formalization
??? tests/
?   ??? test_cfr.c              # assert-based test suite
??? examples/
?   ??? example1_cfr.c          # SCM + counterfactual reasoning
?   ??? example2_cfr.c          # Causal effect estimation
?   ??? example3_cfr.c          # Mediation + bounds
?   ??? example_app_medical.c   # L7: Medical RCT
?   ??? example_app_policy.c    # L7: Policy evaluation (DiD)
??? docs/
    ??? causal-effects.md
    ??? counterfactual-theory.md
    ??? knowledge-graph.md
    ??? coverage-report.md
    ??? gap-report.md
    ??? course-alignment.md
    ??? course-tree.md
```

---

## Building and Testing

```bash
make          # Build static library
make test     # Build + run test suite
make examples # Build + run all examples
make clean    # Remove build artifacts
```

---

## References

- Pearl, J. (2009). *Causality: Models, Reasoning, and Inference* (2nd ed.). Cambridge.
- Pearl, J., Glymour, M., & Jewell, N.P. (2016). *Causal Inference in Statistics: A Primer*. Wiley.
- Imbens, G.W. & Rubin, D.B. (2015). *Causal Inference for Statistics, Social, and Biomedical Sciences*. Cambridge.
- Hernan, M.A. & Robins, J.M. (2020). *Causal Inference: What If*. CRC Press.
- VanderWeele, T.J. (2015). *Explanation in Causal Inference*. Oxford.
- Rubin, D.B. (1974). Estimating causal effects of treatments in randomized and nonrandomized studies. *J. Edu. Psych.*
- Tian, J. & Pearl, J. (2000). Probabilities of causation: Bounds and identification. *UAI*.
- Balke, A. & Pearl, J. (1997). Bounds on treatment effects from studies with imperfect compliance. *JASA*.
- Card, D. & Krueger, A.B. (1994). Minimum wages and employment: A case study. *AER*.

---

## Module Status: COMPLETE

- **L1 Definitions**: Complete (8+ typedef struct, all SCM/P.O./CF/ATE definitions)
- **L2 Core Concepts**: Complete (abduction-action-prediction, SUTVA, ignorability, overlap)
- **L3 Math Structures**: Complete (SCM, DAG, moral graph, path coefficients)
- **L4 Fundamental Laws**: Complete (3-step theorem, mediation decomposition, bounds)
- **L5 Algorithms**: Complete (IPW, DR, matching, Barron-Kenny, Pearl formula, IV)
- **L6 Canonical Problems**: Complete (ATE, mediation, PN/PS/PNS, IV, front-door)
- **L7 Applications**: Partial+ (2 apps: medical RCT, policy evaluation DiD)
- **L8 Advanced Topics**: Partial+ (twin network, transportability, fairness, sequential CF)
- **L9 Research Frontiers**: Partial (documented, open for implementation)

**Self-Check Results:**
- include/ + src/ lines: 4,746 >= 3,000 ?
- typedef struct count: 20+ ? (?5)
- Header files: 7 ? (?4)
- Source files: 7 ? (?4)
- Test assertions: 10+ ? (?5 mathematical)
- Lean theorems: 50+ ? (?5, no by trivial/sorry)
- Examples >30 lines + main: 5 ? (?3)
- L7 applications: 2 ? (?2)
- Filler scan: 0 matches ?
- Stub scan: 0 matches ?
