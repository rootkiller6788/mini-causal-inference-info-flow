# Do-Calculus & Intervention — Causal Inference

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (6 struct typedefs + Lean definitions)
- **L2 Core Concepts**: Complete (6 headers + 6 src files)
- **L3 Math Structures**: Complete (full graph algebra: paths, d-sep, v-structures, moral graph, etc.)
- **L4 Fundamental Laws**: Complete (do-calculus rules, back/front-door, truncated factorization, d-separation, PNS bounds — C + Lean proofs)
- **L5 Algorithms/Methods**: Complete (21+ methods: ATE, IPW, DR, matching, stratification, bootstrap, twin network, etc.)
- **L6 Canonical Problems**: Complete (3 canonical examples >100 lines each: Simpson's paradox resolution with 5 estimators + CI, do-calculus derivation + graph surgery on IV model, counterfactual PN/PS/PNS + mediation + path-specific effects)
- **L7 Applications**: Complete (Clinical Trial Policy Evaluation + Causal Mediation for Health Policy — 652 lines combined)
- **L8 Advanced Topics**: Complete (path-specific effects, fairness, explanation, attributable fraction, sequential back-door, MSM, DiD, RDD)
- **L9 Research Frontiers**: Partial (do-calculus for time series & causal RL documented)

**include/ + src/**: 4,037 lines (467 headers + 3,190 C + 380 Lean) | **Tests**: 21/21 passed | **make**: Clean compile (0 warnings)

## Pearl's Calculus of Intervention

Do-calculus (Pearl, 1995) provides three rules for transforming interventional probabilities into observational ones. This module implements the complete do-calculus framework: Structural Causal Models (SCMs), graph algorithms (d-separation, graph surgery), back-door and front-door criteria, causal effect estimation, and counterfactual reasoning.

### Core Concepts

1. **Structural Causal Model (SCM)**: M = ⟨U, V, F, P(u)⟩ — the mathematical foundation
2. **do-Operator**: do(X=x) replaces the equation for X with X=x (graph surgery)
3. **Truncated Factorization**: P(v|do(x)) = Π P(v_i|pa_i) evaluated at X=x
4. **Three Rules**: Rule 1 (observations), Rule 2 (actions↔observations), Rule 3 (actions)
5. **Back-Door Criterion**: Identifies admissible covariate sets for adjustment
6. **Counterfactuals**: Y_x(u) — what would Y be had X been x in unit u

### Implementation

| File | Purpose |
|------|---------|
| `include/dci_core.h` / `src/dci_core.c` | SCM types, structural equations, SCM factories (Simpson, IV, mediation, m-bias, front-door) |
| `include/dci_graph.h` / `src/dci_graph.c` | Graph algorithms: d-separation, paths, topo sort, moral graph, v-structures, transitive closure |
| `include/dci_do.h` / `src/dci_do.c` | Three do-calculus rules, truncated factorization, identifiability |
| `include/dci_backdoor.h` / `src/dci_backdoor.c` | Back-door criterion, front-door criterion, covariate selection, instrumental variables |
| `include/dci_effect.h` / `src/dci_effect.c` | ATE/ATT/ATU, CATE, NDE/NIE, standardization, IPW, DR, matching, DiD, RDD |
| `include/dci_counterfactual.h` / `src/dci_counterfactual.c` | Counterfactual probability, PN/PS/PNS, twin networks, mediation |
| `src/dci_theory.lean` | Lean 4 formalization: causal graphs, d-separation, mutilation, PNS bounds (decidable, no `sorry`) |
| `src/do_calculus_intervention_app1.c` | L7: Clinical Trial Policy Evaluation (302 lines) |
| `src/do_calculus_intervention_app2.c` | L7: Causal Mediation for Health Policy (350 lines) |

### Key Functions

| Function | Description |
|----------|-------------|
| `dci_scm_add_variable` | Add variable with structural equation to SCM |
| `dci_scm_evaluate` | Evaluate SCM given exogenous values |
| `dci_scm_is_dag` | Verify SCM graph is acyclic |
| `dci_is_d_separated` | Test d-separation (X ⟂ Y | Z) |
| `dci_graph_mutilate` | Remove incoming edges — graph surgery for do(X) |
| `dci_find_v_structures` | Detect collider patterns (X→M←Y) |
| `dci_apply_rule1/2/3` | Apply do-calculus rules |
| `dci_identify_effect` | Determine if causal effect is identifiable |
| `dci_backdoor_find` | Find minimal back-door admissible set |
| `dci_frontdoor_find` | Find front-door mediator set |
| `dci_compute_ate` | Compute Average Treatment Effect |
| `dci_inverse_probability_weighting` | IPW estimator |
| `dci_counterfactual_probability` | Compute P(Y_x = y | evidence) |
| `dci_probability_of_necessity` | Probability that X was necessary for Y |
| `dci_ate_bootstrap` | Bootstrap confidence intervals for ATE |
| `dci_difference_in_differences` | DiD estimator |
| `dci_causal_mediation` | Causal mediation via counterfactuals |
| `dci_counterfactual_fairness` | Counterfactual fairness assessment |
| `dci_path_specific_effect` | Path-specific counterfactual effects |
| `dci_counterfactual_explanation` | Shapley-style feature importance for specific cases |

### Standard SCM Factories

| Factory | Description |
|---------|-------------|
| `dci_scm_simpsons_paradox` | Gender → Treatment, Gender → Outcome, Treatment → Outcome |
| `dci_scm_iv` | Instrument Z → Treatment X → Outcome Y, with confounder U |
| `dci_scm_mediation` | X → M → Y, X → Y (Baron-Kenny mediation) |
| `dci_scm_m_bias` | Collider bias: U1→X, U1→Z, U2→Y, U2→Z |
| `dci_scm_frontdoor` | Front-door: U→X, X→M→Y, U→Y |

### Building

```bash
make && make test && make examples
```

### Module Structure

```
mini-do-calculus-intervention/
  Makefile, README.md
  include/ (6 headers)
  src/      (6 C files + 2 application files + 1 Lean file)
  tests/   test_dci.c (21 tests)
  examples/
    example1_simpsons_paradox.c
    example2_do_calculus.c
    example3_counterfactual.c
  docs/
    do-calculus-theory.md, course-alignment.md, course-tree.md
    knowledge-graph.md, coverage-report.md, gap-report.md
```

### d-Separation in Brief

d-separation (directional separation) is the graphical criterion for reading conditional independence from a causal DAG:

| Structure | Pattern | Blocked by Z? | Unblocked by Z? |
|-----------|---------|---------------|-----------------|
| Chain | X→M→Y | M in Z | M not in Z |
| Fork | X←M→Y | M in Z | M not in Z |
| Collider | X→M←Y | M not in Z, no descendant in Z | M or descendant in Z |

### Dependencies

- C11 compiler (gcc/clang)
- Standard C library (stdlib, math, string, assert)
- Lean 4 (for formal verification only; not required for build)

### References

- Pearl, J. (2009). Causality: Models, Reasoning, and Inference (2nd ed.). Cambridge.
- Pearl, J. (1995). Causal diagrams for empirical research. Biometrika, 82(4):669-710.
- Pearl, J. & Mackenzie, D. (2018). The Book of Why. Basic Books.
- Hernán, M.A. & Robins, J.M. (2020). Causal Inference: What If. Chapman & Hall.
- Spirtes, P., Glymour, C., & Scheines, R. (2000). Causation, Prediction, and Search. MIT Press.
- Morgan, S.L. & Winship, C. (2015). Counterfactuals and Causal Inference (2nd ed.). Cambridge.
- Imbens, G.W. & Rubin, D.B. (2015). Causal Inference for Statistics, Social, and Biomedical Sciences. Cambridge.
