# Course Alignment — Causal Discovery (PC/FCI)

## Nine-School Curriculum Mapping

### MIT
- **6.867 Machine Learning** — Causal inference, graphical models, d-separation
- **6.438 Algorithms for Inference** — Belief propagation, junction trees, CI testing
- **6.S191 Intro to Deep Learning** — Causal representation learning (guest lectures)

### Stanford
- **CS228 Probabilistic Graphical Models** — DAGs, d-separation, Markov equivalence
- **CS229 Machine Learning** — Causal discovery from observational data
- **STATS 361 Causal Inference** — PC algorithm, FCI, constraint-based methods
- **MS&E 326 Causal Inference** — Counterfactuals, graphical criteria

### Berkeley
- **CS 281A Statistical Learning Theory** — Graphical models, CI tests
- **CS 294 Deep Learning** — Causal discovery extensions
- **STAT 260 Causal Inference** — Pearl's framework, do-calculus, graphical models

### CMU
- **36-708 Statistical Methods for Machine Learning** — Graphical models, CI testing
- **10-708 Probabilistic Graphical Models** — Markov properties, d-separation, PC/FCI
- **36-702 Statistical Machine Learning** — Causal structure learning

### Princeton
- **COS 513 Foundations of Probabilistic Modeling** — Bayesian networks, d-separation
- **COS 597C Causal Inference** — Constraint-based discovery, FCI
- **ORF 526 Probability Theory** — Statistical foundations for CI testing

### Caltech
- **CS 156 Machine Learning** — Graphical models, causal inference
- **CS 159 Advanced Topics in ML** — Causal discovery methods
- **CMS/CS/CNS/EE 155** — Probabilistic graphical models

### Cambridge
- **Part II: Machine Learning and Bayesian Inference** — Graphical models, CI
- **Part II: Statistical Methods** — Fisher's z-test, G² test, partial correlation
- **Part III: Advanced Causal Inference** — PC, FCI, score-based methods

### Oxford
- **Advanced Topics in Machine Learning** — Causal graphical models
- **Statistical Machine Learning** — Structure learning, CI tests
- **Causal Inference Reading Group** — Frontiers in causal discovery

### ETH
- **263-5210 Probabilistic Artificial Intelligence** — PGMs, d-separation
- **401-3620 Causal Inference** — Graphical criteria, PC, FCI
- **263-5300 Causal Representation Learning** — Modern causal discovery

## Key Textbooks

| Textbook | Chapters | Topics |
|----------|----------|--------|
| Pearl, "Causality" (2009) | 1-3 | DAGs, d-separation, Markov, causal graphs |
| Spirtes, Glymour, Scheines, "CPS" (2000) | 4-7 | PC algorithm, FCI, completeness proofs |
| Koller & Friedman, "PGMs" (2009) | 3-4 | d-separation, Markov properties, equivalence |
| Peters, Janzing, Scholkopf (2017) | 4-6 | Constraint-based causal discovery |
| Lauritzen, "Graphical Models" (1996) | 2-3 | Moral graph, d-separation algorithms |

## Module Mapping

| Module Component | Course Topic | Course |
|-----------------|-------------|--------|
| cdf_core.h/c | DAG/CPDAG/PAG data structures | CS228, 36-708 |
| cdf_citest.h/c | Fisher z-test, G², partial corr | STATS 361, Part II Stats |
| cdf_graph.h/c | d-separation, skeleton, v-structures | CS228, 10-708 |
| cdf_orientation.h/c | Meek R1-R4, Zhang R1-R10 | COS 597C, STATS 361 |
| cdf_pc.h/c | PC algorithm (3 phases) | CS228, 36-708 |
| cdf_fci.h/c | FCI, RFCI, PDS phase | COS 597C, APML |
| cdf_theory.lean | d-separation, Markov equiv formalization | COS 513, 401-3620 |

## Prerequisites

- Probability theory (conditional probability, Bayes rule)
- Statistics (hypothesis testing, p-values, correlation)
- Graph theory (DAGs, paths, reachability)
- Linear algebra (covariance matrices, matrix inversion)
- C programming (for implementation)