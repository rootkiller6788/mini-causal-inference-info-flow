# Coverage Report -- Do Calculus Intervention

| Level | Name | Coverage | Items |
|-------|------|----------|-------|
| L1 | Definitions | **Complete** | 6 struct typedefs (DCI_SCM, DCI_Variable, DCI_Graph, DCI_Intervention, DCI_CausalEffect, DCI_Counterfactual) + Lean SCM/CausalGraph definitions |
| L2 | Core Concepts | **Complete** | 6 header files (dci_core, dci_graph, dci_do, dci_backdoor, dci_effect, dci_counterfactual) + 6 corresponding src files |
| L3 | Math Structures | **Complete** | Full adjacency matrix, path finding, d-separation, topological sort, v-structures, moral graph, transitive closure |
| L4 | Fundamental Laws | **Complete** | Pearl do-calculus (1995), truncated factorization, back-door theorem (C + Lean), front-door theorem (C + Lean), d-separation soundness (C + Lean), counterfactual consistency (C + Lean), PNS bounds (C + Lean) |
| L5 | Algorithms/Methods | **Complete** | 3 do-calculus rules, graph surgery, back-door adjustment, front-door adjustment, IPW, standardization, doubly-robust, matching, stratification, bootstrap CI, E-value, twin network, PN/PS/PNS, attributable fraction |
| L6 | Canonical Problems | **Complete** | 3 expanded examples (105/125/126 lines): Simpson's Paradox with 5 estimators + CI (example1), Do-calculus derivation + graph surgery on IV model (example2), Counterfactual PN/PS/PNS + twin network + mediation + path-specific effects (example3) |
| L7 | Applications | **Complete** | Clinical Trial Policy Evaluation (app1: 302 lines, 10 analysis sections), Causal Mediation for Health Policy (app2: 350 lines, 13 analysis sections) |
| L8 | Advanced Topics | **Complete** | Path-specific counterfactual effects, attribution, sequential back-door (g-computation), counterfactual fairness, marginal structural models, difference-in-differences, regression discontinuity, doubly-robust estimation, collider bias detection, safe adjustment checking |
| L9 | Research Frontiers | **Partial** | Documented: do-calculus for time series, causal RL; Counterfactual explanation (Shapley-style) implemented |
