# Course Tree — Granger Causality Time Series

## Prerequisites

### Immediate Prerequisites
- Linear algebra: matrix multiplication, inverse, determinant, transpose
- Probability theory: random variables, distributions, expectation
- Statistics: hypothesis testing, F-distribution, p-values, OLS regression
- C programming: dynamic memory, structs, Makefiles
- Time series basics: autocorrelation, stationarity, AR models

### Upstream Modules (mini-complex-control-theory)
- None (this is a leaf module in causal-inference-info-flow group)

### External Prerequisites
- `libm` (math library) — for sqrt, log, exp, cos, sin, fabs
- Standard C11 compiler (GCC or Clang)

## Dependency Graph

```
Probability → Statistics → Hypothesis Testing → Granger F-test (L4)
    ↓
Stochastic Processes → Time Series → AR/VAR Models (L3)
    ↓
Linear Algebra → OLS → VAR Fitting → Granger Test (L5)
    ↓
Information Theory → Transfer Entropy (L5)
    ↓
Spectral Analysis → Geweke Decomposition (L5)
    ↓
Causal Graphs → Graph Theory → DAG Detection (L5)
    ↓
Bootstrap Methods → Bootstrap Granger (L8)
    ↓
Kernel Methods → Nonlinear Granger (L8)
    ↓
Confounder Analysis → Conditional/Partial Granger (L8)
    ↓
Non-stationarity → Time-Varying Granger (L8)
```

## Module Internal Dependencies

```
ts_core (Vector, Matrix, TimeSeries, VARModel, OLS)
  ├── granger_test (F-test, p-value, lag selection)
  │   ├── conditional_granger (Z-conditioning)
  │   ├── partial_granger (orthogonalization)
  │   ├── causality_graph (pairwise tests → graph)
  │   └── granger_bootstrap (resampling inference)
  ├── transfer_entropy (histogram TE, surrogates)
  ├── spectral_granger (frequency-domain, Geweke)
  │   └── partial_granger (PDC from VAR coefs)
  ├── fevd (Cholesky, IRF, variance decomposition)
  ├── nonlinear_granger (kernel RBF, ridge regression)
  └── time_varying_granger (sliding/expanding/rolling windows)
```

## Learning Path

### Beginner Path (L1→L2→L5→L6)
1. Read `include/ts_core.h` — understand Vector, Matrix, TimeSeries, VARModel
2. Run `examples/example1_granger` — see basic Granger test
3. Read `src/ts_core.c` — understand OLS and VAR fitting
4. Run `examples/example3_var` — multivariate VAR

### Intermediate Path (+L3→L4)
5. Read `src/granger_test.c` — understand F-distribution and p-value computation
6. Read `src/transfer_entropy.c` — information-theoretic causality
7. Run `examples/example2_te` — transfer entropy demo
8. Read `src/spectral_granger.c` — frequency-domain decomposition

### Advanced Path (+L7→L8)
9. Read `src/fevd.c` — impulse response and variance decomposition
10. Run `examples/example4_fevd` — FEVD demo
11. Read `src/nonlinear_granger.c` — kernel methods for nonlinear causality
12. Read `src/causality_graph.c` — causal network analysis
13. Run `examples/example5_cgraph` — causality graph demo
14. Read `src/granger_bootstrap.c` — bootstrap inference
15. Run `examples/example6_tvgranger` — time-varying Granger demo

### Applications (L7)
16. Read `src/granger_app_economics.c` — macroeconomic Granger
17. Read `src/granger_app_neuroscience.c` — neural connectivity

### Research Frontiers (L9)
18. Read `docs/knowledge-graph.md` L9 section for future directions
19. Consult `docs/gap-report.md` for open problems

## Course Mapping Summary

| Institution | Key Course | Module Coverage |
|-------------|-----------|-----------------|
| CMU | 36-708 Statistical Methods | GT, FEVD (L5) |
| Stanford | STATS 361 Causal Inference | GT, CG (L4, L8) |
| Berkeley | STAT 248 Analysis of Time Series | SG, FEVD (L5) |
| Cambridge | Part II Time Series | SG (L5) |
| Oxford | Network Science | CGRAPH (L8) |
| ETH | 401-4623-00 Time Series Analysis | SG (L5) |
| MIT | 6.437 Inference and Information | TE (L5) |
| Princeton | ECO 518 Econometrics | GT (L5) |
| Caltech | CNS 187 Neural Computation | NLG (L8) |
