# Course Alignment — Do-Calculus Module

## UCLA: Pearl's Causality Course

| Topic | Coverage |
|-------|----------|
| SCM & Causal Diagrams | dci_core.h/c, dci_graph.h/c |
| do-Calculus Rules | dci_do.h/c |
| Back-Door Criterion | dci_backdoor.h/c |
| Front-Door Criterion | dci_backdoor.h/c |
| Counterfactuals | dci_counterfactual.h/c |
| Causal Effect Estimation | dci_effect.h/c |

## Stanford AA203 / EE363

| Topic | Coverage |
|-------|----------|
| Causal Inference | ATE, ATT, ATU, CATE |
| Instrumental Variables | IV criterion (dci_backdoor) |
| Sensitivity Analysis | E-value, robustness checks |
| Doubly Robust Estimation | dci_effect (IPW + standardization) |

## Harvard EPI 207 / EPI 289

| Topic | Coverage |
|-------|----------|
| Confounding Control | Back-door, front-door, covariate selection |
| Mediation Analysis | NDE, NIE, CDE, causal mediation |
| G-Computation | G-computation, sequential back-door |
| Propensity Scores | Propensity score estimation |

## Textbook Mapping

| Textbook | Chapter | Module |
|----------|---------|--------|
| Pearl (2009) | Ch.3-4 | dci_do, dci_backdoor, dci_graph |
| Pearl (2009) | Ch.7 | dci_counterfactual |
| Hernán & Robins (2020) | Ch.13-15 | dci_effect, dci_backdoor |
| Morgan & Winship (2015) | Ch.4-7 | dci_effect, estimation methods |
