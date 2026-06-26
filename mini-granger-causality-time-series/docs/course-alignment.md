# Course Alignment — Granger Causality Time Series

Nine-school course mapping for Granger causality and related time-series
causal inference methods.

## MIT

| Course | Topic | Coverage |
|--------|-------|----------|
| 6.045/18.400 Automata, Computability, Complexity | Information theory foundations | TE (L2) |
| 6.437 Inference and Information | Directed information, transfer entropy | TE (L5) |
| 6.841 Advanced Complexity | Causal inference complexity | — (L9 future) |
| 15.077 Statistical Learning | VAR models, AIC/BIC selection | VAR (L5) |

## Stanford

| Course | Topic | Coverage |
|--------|-------|----------|
| STATS 361 Causal Inference | Granger causality, potential outcomes | GT (L4), CG (L8) |
| EE 278 Probability | Stochastic processes, time series | SIM (L6) |
| CS 229 Machine Learning | Feature-based causality | NLG (L8) |
| STATS 315B Modern Applied Statistics | Bootstrap, resampling | BOOT (L8) |

## Berkeley

| Course | Topic | Coverage |
|--------|-------|----------|
| STAT 248 Analysis of Time Series | VAR, spectral analysis | SG (L5), FEVD (L5) |
| STAT 210A Theoretical Statistics | F-test theory, asymptotics | GT (L4) |
| INFO 290 Causal Inference | DAGs, causal graphs | CGRAPH (L5) |

## CMU

| Course | Topic | Coverage |
|--------|-------|----------|
| 36-708 Statistical Methods | Granger, VAR, FEVD | GT, FEVD (L5) |
| 10-708 Probabilistic Graphical Models | Causal discovery, DAGs | CGRAPH (L5) |
| 36-705 Intermediate Statistics | Bootstrap inference | BOOT (L8) |
| 11-785 Deep Learning | Neural Granger | — (L9 future) |

## Princeton

| Course | Topic | Coverage |
|--------|-------|----------|
| ORF 524 Statistical Theory | Time series asymptotics | GT (L4) |
| ECO 518 Econometrics | Granger causality tests | GT (L5) |
| COS 513 Foundations of Probabilistic Modeling | Graphical models | CGRAPH (L5) |

## Caltech

| Course | Topic | Coverage |
|--------|-------|----------|
| CMS 117 Probability and Random Processes | Stochastic processes | SIM (L6) |
| CNS 187 Neural Computation | Neural time series | NLG (L8) |
| ACM 216 Markov Chains | State-space models | TVG (L8) |

## Cambridge

| Course | Topic | Coverage |
|--------|-------|----------|
| Part II Time Series | VAR, spectral analysis | SG (L5) |
| Part III Statistical Theory | Bootstrap, resampling | BOOT (L8) |
| MPhil Economics | Applied Granger | GT (L5) |

## Oxford

| Course | Topic | Coverage |
|--------|-------|----------|
| Statistical Machine Learning | Causal discovery | CGRAPH (L5) |
| Time Series and Forecasting | VAR, forecasting | FEVD (L5) |
| Network Science | Causal networks | CGRAPH (L8) |

## ETH

| Course | Topic | Coverage |
|--------|-------|----------|
| 401-3628-14 Bayesian Time Series | Dynamic models | TVG (L8) |
| 401-4623-00 Time Series Analysis | Spectral methods | SG (L5) |
| 401-6282-00 Statistical Modeling | Causal inference | CG, PG (L8) |

## Key: Course → Module Mapping

| Abbrev | Module |
|--------|--------|
| GT | granger_test (Granger F-test) |
| TE | transfer_entropy (Schreiber TE) |
| SG | spectral_granger (Geweke) |
| CG | conditional_granger |
| PG | partial_granger (Guo et al.) |
| FEVD | fevd (Sims, Pesaran-Shin) |
| NLG | nonlinear_granger (kernel) |
| CGRAPH | causality_graph (Eichler) |
| TVG | time_varying_granger |
| BOOT | granger_bootstrap |
| VAR | ts_core (VAR fitting, AIC/BIC) |
| SIM | ts_core (AR(1)/VAR(1) simulation) |
