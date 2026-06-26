# Course Tree -- Causal Discovery PC/FCI

## Prerequisites (this module depends on)
- Probability theory: random variables, conditional independence
- Statistics: hypothesis testing, correlation, regression
- Graph theory: DAGs, paths, reachability
- Linear algebra: matrix inversion, Gaussian elimination
- C programming: dynamic memory, linear algebra

## Dependencies within this project
- `cdf_core` → no dependencies (foundation types)
- `cdf_citest` → `cdf_core` (uses Dataset)
- `cdf_graph` → `cdf_core`, `cdf_citest` (uses CI tests)
- `cdf_orientation` → `cdf_core`, `cdf_graph`
- `cdf_pc` → `cdf_core`, `cdf_citest`, `cdf_graph`, `cdf_orientation`
- `cdf_fci` → `cdf_core`, `cdf_citest`, `cdf_graph`, `cdf_orientation`, `cdf_pc`
- `cdf_bootstrap` → `cdf_core`, `cdf_pc`, `cdf_fci`
- `cdf_ida` → `cdf_core`, `cdf_pc`, `cdf_graph`
- `cdf_kernel_ci` → `cdf_core`, `cdf_citest`
- `cdf_scores` → `cdf_core`

## Course Alignment
| School | Course | Topics |
|--------|--------|--------|
| CMU | 36-708 Statistical Methods | PC/FCI, constraint-based causal discovery |
| CMU | 10-708 Probabilistic Graphical Models | d-separation, Markov equivalence |
| Stanford | STATS 361 Causal Inference | PC, FCI, IDA, do-calculus |
| Berkeley | STAT 260 Causal Inference | Graphical models, identifiability |
| ETH | 401-3620-00 Causal Representation Learning | Constraint-based methods |
| UCLA | STATS 200C Causal Inference | PC algorithm, FCI, causal discovery |
