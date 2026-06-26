# Gap Report — Pearl Structural Causal Models

## Current Gaps (L9 Only)

All L1-L8 levels are complete. Only L9 (Research Frontiers) has gaps,
which is expected per SKILL.md §6.1 (L9 requires only Partial).

| # | Gap | Level | Priority | Effort | Why Not Implemented |
|---|-----|-------|----------|--------|---------------------|
| 1 | Causal Representation Learning | L9 | Low | High | Research-stage, no stable algorithm consensus |
| 2 | Neural SCMs (VAE-based) | L9 | Low | High | Requires deep learning framework dependencies |
| 3 | Counterfactual Generative Models | L9 | Low | High | GAN/VAE-based, beyond scope of C library |
| 4 | Double ML with Neural Networks | L9 | Low | Medium | Simplified linear version implemented |
| 5 | Causal Discovery from Heterogeneous Data | L9 | Low | Medium | Pooled data assumption in current impl |
| 6 | Causal Reinforcement Learning | L9 | Low | High | Requires RL environment + policy optimization |

## Resolved Gaps

All previous gaps (L1-L8) have been resolved:

| Gap | Resolution |
|-----|------------|
| L4 Fundamental Laws (was Missing) | Full do-calculus (3 rules), back-door, front-door, counterfactual consistency implemented in C + Lean |
| L5 Algorithms (was Partial) | 17 algorithms: PC, GES, FCI, Hill-Climbing, Tabu, Bootstrap, BIC/AIC, 2SLS, PSM, Double ML, g-formula, DR |
| L6 Canonical Problems (was Missing) | 8 problems: Simpson, M-bias, Berkson, Smoking-Cancer, Complier ACE, PNS, Back-door, Front-door |
| L7 Applications (was Missing) | 8 applications: Clinical trials, IV, PSM, DiD, RDD, Mediation, Double ML, Simpson resolution |
| L8 Advanced Topics (was Missing) | 8 topics: Transportability, Selection bias, Fairness, Sensitivity, Manski, g-formula, MSM, DR |

## No Known TODOs, FIXMEs, stubs, or placeholders

All source files compile cleanly with `-Wall -Wextra` (no warnings).
All Lean proofs are substantive (no `by trivial`, no `:= 0.0`, no `sorry`).