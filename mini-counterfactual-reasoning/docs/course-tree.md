# Course Tree — mini-counterfactual-reasoning

## Prerequisite Dependencies

```
mini-counterfactual-reasoning
│
├── [PREREQ] Probability Theory (L1)
│   ├── Random variables, expectation, variance
│   ├── Conditional probability, Bayes rule
│   └── Joint/marginal/conditional distributions
│       └── Used by: ATE estimation, IPW, Bayesian inference
│
├── [PREREQ] Graph Theory (L1)
│   ├── DAGs, topological ordering
│   ├── Paths, ancestors, descendants
│   └── Moral graphs, separation
│       └── Used by: SCM DAG analysis, d-separation
│
├── [PREREQ] Linear Algebra (L1)
│   ├── Covariance, correlation matrices
│   ├── Linear regression (OLS)
│   └── Matrix operations
│       └── Used by: PC algorithm, Wald estimator
│
├── [PREREQ] Statistics (L2)
│   ├── Hypothesis testing, p-values
│   ├── Confidence intervals
│   ├── Multiple testing correction
│   └── Bootstrap resampling
│       └── Used by: effect estimation, sensitivity analysis
│
├── [PREREQ] Optimization (L1)
│   ├── Linear programming basics
│   └── Gradient descent intuition
│       └── Used by: Balke-Pearl LP bounds
│
└── [PREREQ] Formal Logic (L1)
    ├── Boolean algebra
    └── Quantifiers, implications
        └── Used by: Lean formalization, do-calculus
```

## Internal Dependency Graph

```
cfr_core.h/c (SCM + DAG + interventions)
├── cfr_potential.h/c (RCM + dataset generation)
│   ├── cfr_effects.h/c (ATE estimation methods)
│   │   └── cfr_analysis.h/c (analysis tools)
│   └── cfr_mediation.h/c (mediation analysis)
├── cfr_counterfactual.h/c (3-step CF + twin network)
│   ├── uses: cfr_core (SCM cloning, intervention)
│   └── cfr_bounds.h/c (PN/PS/PNS bounds)
│       └── uses: cfr_core (SCM operations)
└── cfr_mediation.h/c
    └── uses: cfr_core (SCM, path analysis)
```

## Knowledge Progression (Learning Path)

### Stage 1: Foundations (L1-L3)
1. CausalVariable, StructuralEquation, SCM definitions → cfr_core.h
2. SCM creation, evaluation, topological computation → cfr_core.c
3. DAG analysis: parents, children, ancestors, descendants → cfr_core.c
4. D-separation and moral graphs → cfr_core.c

### Stage 2: Interventions & Counterfactuals (L2-L4)
5. Hard/soft interventions and do-calculus → cfr_core.c + Lean
6. 3-step counterfactual algorithm (Abduction → Action → Prediction) → cfr_counterfactual.c
7. Potential outcomes framework + SUTVA → cfr_potential.c
8. Twin networks for cross-world reasoning → cfr_counterfactual.c

### Stage 3: Effect Estimation (L4-L5)
9. Oracle ATE (ground truth) → cfr_effects.c
10. Difference in means (RCT) → cfr_effects.c
11. IPW estimation → cfr_effects.c
12. Doubly robust estimation → cfr_effects.c
13. Stratification and matching → cfr_effects.c
14. Propensity score → cfr_potential.c

### Stage 4: Mediation & Decomposition (L4-L6)
15. Baron-Kenny product method → cfr_mediation.c
16. Pearl mediation formula → cfr_mediation.c
17. Controlled direct effects → cfr_mediation.c
18. Mediation sensitivity analysis → cfr_mediation.c

### Stage 5: Bounds & Sensitivity (L4-L6)
19. Trivial bounds (no assumptions) → cfr_bounds.c
20. Experimental bounds (RCT) → cfr_bounds.c
21. Monotonicity bounds → cfr_bounds.c
22. Tian-Pearl combined bounds → cfr_bounds.c
23. Balke-Pearl IV bounds → cfr_bounds.c
24. Manski trimming → cfr_bounds.c

### Stage 6: Advanced Topics (L7-L8)
25. Instrumental variables + compliance types → cfr_analysis.c
26. Front-door adjustment → cfr_analysis.c
27. Placebo refutation tests → cfr_analysis.c
28. Transportability → cfr_counterfactual.c
29. Counterfactual fairness → cfr_counterfactual.c
30. Sequential counterfactuals → cfr_counterfactual.c
31. Moderated mediation → cfr_mediation.c

### Stage 7: Research Frontiers (L9)
32. PC algorithm skeleton → cfr_core.c
33. Causal discovery literature → docs/
34. Meta-causal analysis → docs/
35. Causal representation learning → docs/

---

## Dependency Notes

- **Minimal external dependencies**: Standard C library only (`stdio.h`, `stdlib.h`, `string.h`, `math.h`, `assert.h`)
- **Lean 4**: Pure core (no Mathlib), uses `rfl`, `cases`, `simp`, `linarith`, `omega`, `decide`
- **Self-contained**: All causal concepts defined within the module
- **No circular dependencies**: All includes are acyclic (cfr_core → specialized modules)
