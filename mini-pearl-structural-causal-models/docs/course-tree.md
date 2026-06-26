# Course Tree — Pearl Structural Causal Models

## Prerequisite Dependency Graph

```
                    ┌──────────────────────────┐
                    │  Probability & Statistics │
                    │  (conditional prob, CI)   │
                    └────────────┬─────────────┘
                                 │
                    ┌────────────▼─────────────┐
                    │   Graph Theory            │
                    │   (DAG, paths, cycles)    │
                    └────────────┬─────────────┘
                                 │
         ┌───────────────────────┼───────────────────────┐
         │                       │                       │
┌────────▼────────┐   ┌─────────▼──────────┐   ┌────────▼────────┐
│  L1: SCM        │   │  L2: d-separation  │   │  L3: Linear SEM │
│  Definitions    │   │  Markov condition  │   │  Cov/Precision   │
└────────┬────────┘   └─────────┬──────────┘   └────────┬────────┘
         │                       │                       │
         └───────────────────────┼───────────────────────┘
                                 │
                    ┌────────────▼─────────────┐
                    │  L4: do-Calculus (3 rules)│
                    │  Back-door / Front-door   │
                    │  Counterfactual Consistency│
                    └────────────┬─────────────┘
                                 │
         ┌───────────────────────┼───────────────────────┐
         │                       │                       │
┌────────▼────────┐   ┌─────────▼──────────┐   ┌────────▼────────┐
│  L5: Structure  │   │  L5: Causal Effect │   │  L5: Scoring     │
│  Learning       │   │  Estimation        │   │  BIC/AIC         │
│  PC/GES/FCI/HC  │   │  IV/PSM/DiD/RDD    │   │  Bootstrap       │
└────────┬────────┘   └─────────┬──────────┘   └────────┬────────┘
         │                       │                       │
         └───────────────────────┼───────────────────────┘
                                 │
                    ┌────────────▼─────────────┐
                    │  L6: Canonical Problems   │
                    │  Simpson/M-bias/Berkson   │
                    │  Smoking-Cancer SCM       │
                    └────────────┬─────────────┘
                                 │
         ┌───────────────────────┼───────────────────────┐
         │                       │                       │
┌────────▼────────┐   ┌─────────▼──────────┐   ┌────────▼────────┐
│  L7: Clinical   │   │  L7: Econometric   │   │  L7: Policy      │
│  Trials/PSM     │   │  IV/2SLS/Double ML │   │  DiD/RDD          │
└────────┬────────┘   └─────────┬──────────┘   └────────┬────────┘
         │                       │                       │
         └───────────────────────┼───────────────────────┘
                                 │
                    ┌────────────▼─────────────┐
                    │  L8: Advanced Topics      │
                    │  Transportability         │
                    │  Selection Bias (Heckman) │
                    │  Counterfactual Fairness  │
                    │  Sensitivity (E-value)    │
                    │  g-formula / MSM / DR     │
                    └────────────┬─────────────┘
                                 │
                    ┌────────────▼─────────────┐
                    │  L9: Research Frontiers   │
                    │  Causal Repr. Learning    │
                    │  Neural SCMs              │
                    │  Causal RL                │
                    └──────────────────────────┘
```

## Linear Learning Path

### Phase 1: Foundations (L1-L2)
1. Read Pearl (2009) §1-2 — SCM definition, causal graphs
2. Implement: `scm_create()`, `scm_add_variable()`, `scm_add_edge()`
3. Understand: DAG property, topological ordering
4. Practice: Build simple 3-variable SCMs

### Phase 2: Graphical Criteria (L2-L3)
1. Read Pearl (2009) §1.2-1.3 — d-separation criterion
2. Implement: `scm_d_separated()`, `scm_moral_graph()`, `scm_find_colliders()`
3. Practice: Verify d-separation on small graphs
4. Lean: Formalize d-separation in `pearl_scm.lean`

### Phase 3: do-Calculus (L4)
1. Read Pearl (2009) §3 — do-operator and do-calculus
2. Implement: `scm_do_intervention()`, 3 do-calculus rule checks
3. Implement: `scm_back_door()`, `scm_front_door()`
4. Practice: Compute causal effects with back-door adjustment

### Phase 4: Counterfactuals (L4-L6)
1. Read Pearl (2000) §7 — 3-step counterfactual process
2. Implement: `scm_counterfactual_compute()`, 3-step helpers
3. Implement: `scm_prob_necessity_sufficiency()`
4. Practice: Smoking-cancer example with front-door

### Phase 5: Structure Learning (L5)
1. Read Spirtes et al. (2000) §5-6 — PC and FCI algorithms
2. Read Chickering (2002) — GES algorithm
3. Implement: PC, GES, FCI, hill-climbing, tabu search
4. Practice: Learn causal graph from synthetic data

### Phase 6: Applications (L7)
1. Read Hernan & Robins (2020) — clinical trials, PSM
2. Read Angrist & Pischke (2009) — IV, DiD, RDD
3. Implement: clinical trial analysis, IV/2SLS, PSM, DiD, RDD
4. Practice: Analyze simulated clinical trial data

### Phase 7: Advanced Topics (L8)
1. Read Bareinboim & Pearl (2016) — transportability
2. Read Heckman (1979) — selection bias
3. Read VanderWeele & Ding (2017) — sensitivity analysis
4. Implement: transportability, Heckman correction, E-value, fairness
5. Practice: Sensitivity analysis on effect estimates

### Phase 8: Frontiers (L9)
1. Read Schölkopf et al. (2021) — causal representation learning
2. Explore: Neural SCMs, causal RL
3. This level is documented but not fully implemented (research-stage)

## Nine-School Course Prerequisites

| School | Key Course | Prerequisites for This Module |
|--------|-----------|-------------------------------|
| **MIT** | 6.S978 | 6.041 (Probability), 6.006 (Algorithms) |
| **Stanford** | CS228 | CS229 (ML), STATS 217 (Stochastic Processes) |
| **Berkeley** | CS294 | CS281A (PGM), STAT 210A (Theoretical Stats) |
| **CMU** | 10-708 | 10-601 (ML), 36-705 (Intermediate Stats) |
| **Princeton** | COS 597C | COS 324 (ML), ORF 405 (Regression) |
| **Caltech** | CS/CNS 155 | CMS/CS 156a (Learning), ACM 116 (Probability) |
| **Cambridge** | Part III SL | Part II Probability, Part II Statistical Modelling |
| **Oxford** | Adv. Causal Inf. | B8.1 Probability, B8.3 Statistical ML |
| **ETH** | 263-5300 | 401-2604 (Probability), 401-3632 (Statistical Modelling) |

## Module Dependency Map (within mini-theory-of-computation)

```
mini-probability-theory/
  └── mini-pearl-structural-causal-models/   ← This module
        └── (future) mini-pcp-theorem/
        └── (future) mini-causal-discovery/
```
