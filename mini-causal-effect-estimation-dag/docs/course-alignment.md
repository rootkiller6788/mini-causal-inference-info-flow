# Course Alignment — Causal Effect Estimation Dag

## Nine-School Course Chapter Mapping

| School | Course | Relevant Chapters |
|--------|--------|-------------------|
| **MIT** | 6.045/18.400 Automata, Computability, Complexity | DAGs as computational models; graph reachability as decision problems; topological sort as total order |
| **MIT** | 6.841 Advanced Complexity | Causal identification as proof complexity; do-calculus as rewriting system completeness |
| **Stanford** | CS254 Computational Complexity | Causal graphs as Boolean circuit analogs; backdoor paths as structural restrictions |
| **Stanford** | CS358 Circuit Complexity | DAG structure as circuit DAG; depth/width tradeoffs in adjustment set search |
| **Berkeley** | CS172 Computability & Complexity | d-separation as graph decision procedure; minimal adjustment as subset search |
| **Berkeley** | CS278 Computational Complexity | PCP-like decomposition in doubly robust estimation; bootstrap as interactive proof analog |
| **CMU** | 15-455 Undergraduate Complexity | NP-hardness of exact covariate selection; O(2^V) adjustment subset enumeration |
| **CMU** | 15-855 Graduate Complexity | Algorithm-complexity bridge: polynomial-time backdoor vs exponential frontdoor search |
| **Princeton** | COS 522 Computational Complexity | Causal identifiability as complexity class; do-calculus as derivation system |
| **Princeton** | COS 551 Advanced Complexity | Cryptographic applications: unmeasured confounding as adversarial hidden state |
| **Caltech** | CS 151 Complexity Theory | Physical causal models as computational systems; thermodynamics of information flow in DAGs |
| **Caltech** | CS 154 Limits of Computation | Information-theoretic bounds on causal effect estimation from finite samples |
| **Cambridge** | Part II Complexity Theory | Structural equation models as deterministic computation; ACE as semantic value |
| **Cambridge** | Part III Advanced Complexity | Mediation analysis as computational decomposition; NDE/NIE as parallel reduction |
| **Oxford** | Computational Complexity | Causal effect bounds as approximation complexity; E-value as hardness of approximation |
| **Oxford** | Advanced Complexity Theory | Quantum causal models; super-dense coding as causal channel capacity |
| **ETH** | 263-4650 Advanced Complexity | Logic of do-calculus; formal proofs in type theory; Lean 4 formalization |
| **ETH** | 252-0400 Logic & Computation | Back-door as logical consequence; d-separation as independence logic |

## Concept-to-Course Mapping

| Concept | Primary Course | Secondary Course |
|---------|---------------|------------------|
| DAG and topological sort | MIT 6.045 | Stanford CS358 |
| d-separation | Berkeley CS172 | ETH 252-0400 |
| Back-door criterion | Harvard STAT 234 | Columbia P8109 |
| Front-door criterion | UCLA STAT 232 | Berkeley STAT 260 |
| do-calculus | MIT 6.841 | ETH 263-4650 |
| Propensity score | Stanford MS&E 327 | Harvard STAT 234 |
| IPW estimator | Harvard STAT 234 | Columbia P8109 |
| G-Computation | Berkeley STAT 260 | UCLA STAT 232 |
| Doubly Robust | Harvard STAT 234 | CMU 15-855 |
| Matching | Stanford MS&E 327 | Harvard STAT 234 |
| Rosenbaum bounds | Harvard STAT 234 | Columbia P8109 |
| E-value | Harvard STAT 234 | Cambridge Part III |
| Mediation (BK/IKT) | UCLA STAT 232 | Berkeley STAT 260 |
| Longitudinal mediation | Harvard STAT 234 | CMU 15-855 |

## Implementation Coverage

All core concepts from the nine-school curriculum are implemented:
- **Graph algorithms**: DAG, topological sort, cycle detection, path finding, d-separation
- **Causal identification**: Back-door, front-door, do-calculus Rule 2
- **Adjustment sets**: Minimal, optimal, parent set, efficient selection
- **Estimation methods**: IPW, G-Computation, Doubly Robust, Stratification, Matching
- **Sensitivity analysis**: Rosenbaum bounds, E-value, Cornfield conditions
- **Mediation analysis**: Baron-Kenny, IKT simulation, binary mediator, longitudinal
- **Formal verification**: Lean 4 with 19 theorem sections
