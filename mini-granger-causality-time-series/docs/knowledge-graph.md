# Knowledge Graph — Granger Causality Time Series

## L1: Definitions (Complete — 11 type definitions)

| # | Definition | File | Type |
|---|-----------|------|------|
| 1 | Time Series (Vector-valued) | ts_core.h | `TimeSeries` |
| 2 | Vector Autoregression (VAR) | ts_core.h | `VARModel` |
| 3 | Vector (linear algebra) | ts_core.h | `Vector` |
| 4 | Matrix (linear algebra) | ts_core.h | `Matrix` |
| 5 | Granger Test Result | granger_test.h | `GrangerTestResult` |
| 6 | Transfer Entropy | transfer_entropy.h | `TransferEntropy` |
| 7 | Spectral Granger | spectral_granger.h | `SpectralGranger` |
| 8 | Conditional Granger | conditional_granger.h | `ConditionalGranger` |
| 9 | Causality Graph | causality_graph.h | `CausalityGraph` |
| 10 | Forecast Error Variance Decomp. | fevd.h | `FEVDResult` |
| 11 | Nonlinear Granger Result | nonlinear_granger.h | `NonlinearGrangerResult` |
| 12 | Time-Varying Granger | time_varying_granger.h | `TVGrangerResult` |
| 13 | Bootstrap Result | granger_bootstrap.h | `BootstrapGrangerResult` |
| 14 | Partial Directed Coherence | partial_granger.h | `PDCResult` |
| 15 | Partial Granger Result | partial_granger.h | `PartialGrangerResult` |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | VAR(p) model: Y_t = ΣA_k Y_{t-k} + ε_t | ts_core.c: var_fit(), var_forecast() |
| 2 | Restricted vs. unrestricted prediction | granger_test.c: gt_test() |
| 3 | Granger causality: X→Y iff past X improves Y prediction | granger_test.c |
| 4 | Transfer entropy: I(X^t→Y_{t+1} | Y^t) | transfer_entropy.c |
| 5 | Frequency-domain causality (Geweke 1982) | spectral_granger.c |
| 6 | Conditional independence | conditional_granger.c |
| 7 | Partial directed coherence | partial_granger.c |
| 8 | Impulse response analysis | fevd.c |
| 9 | Causal network topology | causality_graph.c |
| 10 | Time-varying causal relationships | time_varying_granger.c |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Euclidean vector space (R^n) | ts_core.c: vec_create, vec_norm, vec_dot |
| 2 | Matrix algebra (R^{m×n}) | ts_core.c: mat_mul, mat_inverse, mat_det |
| 3 | OLS normal equations: β̂ = (X'X)⁻¹X'y | ts_core.c: ols_estimate() |
| 4 | Cholesky decomposition: A = LL' | fevd.c: cholesky_decompose() |
| 5 | F-distribution (regularized incomplete beta) | granger_test.c: beta_cf() |
| 6 | Gaussian RBF kernel: K(x,y)=exp(-||x-y||²/2σ²) | nonlinear_granger.c: kernel_rbf() |
| 7 | Complex transfer function H(ω) | spectral_granger.c: transfer_function() |
| 8 | Graph theory: DAG, adjacency, betweenness | causality_graph.c |
| 9 | Bootstrap distribution (Efron 1979) | granger_bootstrap.c |
| 10 | Sliding window decomposition | time_varying_granger.c |

## L4: Fundamental Theorems (Complete)

| # | Theorem | Reference | Implementation |
|---|---------|-----------|---------------|
| 1 | Granger Representation Theorem | Granger (1969) | granger_test.c, granger_causality.lean |
| 2 | Geweke Spectral Decomposition | Geweke (1982) | spectral_granger.c: sg_geweke_decomposition() |
| 3 | Sims FEVD Theorem | Sims (1980) | fevd.c: fevd_compute() |
| 4 | Schreiber Transfer Entropy | Schreiber (2000) | transfer_entropy.c |
| 5 | Bootstrap Consistency | Efron & Tibshirani (1993) | granger_bootstrap.c |
| 6 | Pesaran-Shin Generalized FEVD | Pesaran & Shin (1998) | fevd.c: fevd_generalized() |
| 7 | AIC/BIC Model Selection | Akaike (1974), Schwarz (1978) | ts_core.c: var_order_select() |
| 8 | Wold Decomposition | Wold (1938) | fevd.c: compute_irf() |
| 9 | Granger Causality Asymmetry | Granger (1969) | granger_causality.lean: granger_false_when_p_ge_T |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Complexity | Implementation |
|---|-----------|------------|---------------|
| 1 | VAR fitting (equation-by-equation OLS) | O(T·k²) | ts_core.c: var_fit() |
| 2 | Granger F-test with p-value | O(T·p) | granger_test.c: gt_test() |
| 3 | Geweke spectral decomposition | O(n_freqs·d³) | spectral_granger.c |
| 4 | Transfer entropy histogram estimation | O(T·bins) | transfer_entropy.c: te_compute() |
| 5 | Kernel ridge regression | O(n³) | nonlinear_granger.c |
| 6 | Kernel median heuristic | O(n²) | nonlinear_granger.c |
| 7 | Block bootstrap resampling | O(B·T·p²) | granger_bootstrap.c |
| 8 | Sliding window Granger | O(n_win·T·p) | time_varying_granger.c |
| 9 | Kahn's topological sort (DAG check) | O(V+E) | causality_graph.c |
| 10 | Brandes betweenness centrality | O(VE) | causality_graph.c |
| 11 | Cholesky decomposition | O(d³) | fevd.c |
| 12 | PDC computation | O(n_freqs·d³) | partial_granger.c |
| 13 | Partial Granger orthogonalization | O(T) | partial_granger.c |
| 14 | Surrogate significance testing | O(S·T) | transfer_entropy.c |

## L6: Canonical Problems (Complete)

| # | Problem | Implementation |
|---|---------|---------------|
| 1 | AR(1) simulation: y_t = φ y_{t-1} + ε_t | ts_core.c: ts_simulate_ar1() |
| 2 | Coupled VAR(1) Granger detection | ts_core.c: ts_simulate_var1_coupled() |
| 3 | GDP vs Money Supply (Economics) | granger_app_economics.c |
| 4 | LFP/EEG connectivity (Neuroscience) | granger_app_neuroscience.c |
| 5 | Causality graph DAG detection | causality_graph.c: graph_is_dag() |
| 6 | Spectral band decomposition | example4_fevd.c |

## L7: Applications (Partial+ — 4 items)

| # | Application | Implementation | Keywords |
|---|-------------|---------------|----------|
| 1 | Macroeconomic forecasting | granger_app_economics.c | GDP, money supply |
| 2 | Brain connectivity (EEG/LFP) | granger_app_neuroscience.c | LFP, EEG |
| 3 | End-to-end VAR forecasting demo | example3_var.c | VAR, AIC |
| 4 | FEVD impulse response demo | example4_fevd.c | FEVD, Cholesky |

## L8: Advanced Topics (Complete — 6 items)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Nonlinear Granger (kernel methods) | nonlinear_granger.c |
| 2 | Conditional Granger (confounder control) | conditional_granger.c |
| 3 | Partial Granger (exogenous control) | partial_granger.c |
| 4 | Time-varying Granger (sliding/expanding/rolling) | time_varying_granger.c |
| 5 | Bootstrap inference (parametric + block) | granger_bootstrap.c |
| 6 | Causality graph analysis (hubs, DAG, betweenness) | causality_graph.c |

## L9: Research Frontiers (Partial — documented)

| # | Frontier | Status |
|---|----------|--------|
| 1 | Deep Granger (neural network-based) | Documented |
| 2 | Neural Granger causality | Documented |
| 3 | Multi-scale transfer entropy | Documented |
| 4 | Causal time series ML | Documented |
