# Mini Causal Inference & Information Flow

A collection of **from-scratch, zero-dependency C implementations** of causal inference theory and information-theoretic causality. Each module translates Pearl's structural causal models, do-calculus, counterfactual reasoning, Granger causality, and Schreiber's transfer entropy equations into runnable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-causal-discovery-pc-fci](mini-causal-discovery-pc-fci/) | PC/FCI algorithms, conditional independence tests, constraint-based causal structure discovery, bootstrap stability | Stanford STATS 366, CMU 10-708 |
| [mini-causal-effect-estimation-dag](mini-causal-effect-estimation-dag/) | SCM, DAG adjustment sets, backdoor criterion, frontdoor criterion, causal effect identification | Harvard EPI 289, MIT 6.438 |
| [mini-counterfactual-reasoning](mini-counterfactual-reasoning/) | Abduction-action-prediction, counterfactual bounds, potential outcomes, mediation analysis | Harvard EPI 289, Stanford MS&E 338 |
| [mini-do-calculus-intervention](mini-do-calculus-intervention/) | Do-calculus 3 rules, graph surgery, intervention, identifiability, causal effect estimation | Stanford CS 229M, Harvard EPI 289 |
| [mini-granger-causality-time-series](mini-granger-causality-time-series/) | Granger F-test, spectral/nonlinear/conditional Granger, FEVD, VAR causality | MIT 14.384, LSE EC484 |
| [mini-info-flow-in-networks](mini-info-flow-in-networks/) | Transfer entropy, Granger causality, mutual information on network topologies, causal network reconstruction | MIT 6.441, MIT 6.207 |
| [mini-pearl-structural-causal-models](mini-pearl-structural-causal-models/) | SCM formalism, transportability, selection bias correction, g-formula, counterfactual fairness | Stanford STATS 362, Harvard EPI 289 |
| [mini-transfer-entropy-schreiber](mini-transfer-entropy-schreiber/) | Schreiber transfer entropy, information-theoretic causality, embedding, statistical significance | MIT 6.441, MPI-PKS |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module maps directly to canonical references (Pearl 2009, Spirtes et al. 2000, Schreiber 2000, Granger 1969)
- **Practical demos** — causal graph discovery from CSV data, do-calculus identifiability checker, Granger causality significance tests, transfer entropy estimator, and more

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-causal-discovery-pc-fci
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-causal-inference-info-flow/
├── mini-causal-discovery-pc-fci/       # PC & FCI constraint-based causal discovery
├── mini-causal-effect-estimation-dag/  # DAG-based causal effect estimation & adjustment
├── mini-counterfactual-reasoning/      # Pearl's 3-step counterfactual reasoning
├── mini-do-calculus-intervention/      # Do-calculus, interventions, graph surgery
├── mini-granger-causality-time-series/ # Granger causality for time series analysis
├── mini-info-flow-in-networks/         # Information flow & causality on networks
├── mini-pearl-structural-causal-models/# Pearl's SCM formalism & advanced topics
└── mini-transfer-entropy-schreiber/    # Schreiber transfer entropy & information transfer
```

## License

MIT
