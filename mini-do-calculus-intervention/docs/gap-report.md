# Gap Report -- Do Calculus Intervention

## Completed Gaps (2026-06-18)

- L4: Lean 4 formalization updated — removed `sorry`, added non-trivial theorems with `decide`/`cases`/`rfl` proofs (graph mutilation, d-separation, v-structures, PNS bounds)
- L6: Examples expanded from stub (<30 lines) to substantial end-to-end demonstrations (105/125/126 lines each with multiple estimation methods + full derivations)
- L7: Two real application files (app1: Clinical Trial Policy Evaluation, app2: Causal Mediation for Health Policy) replacing 4-line stubs
- L8: Advanced methods implemented (path-specific effects, counterfactual fairness, Shapley-style explanation, attributable fraction, sequential back-door)

## Remaining Gaps

| Level | Gap | Priority | Notes |
|-------|-----|----------|-------|
| L9 | Do-calculus for time series | Low | Partial — g-computation and sequential back-door implemented; full time-series do-calculus needs longitudinal DAG expansion |
| L9 | Causal reinforcement learning | Low | Documented only; implementation requires RL framework beyond Mini scope |
| L8 | Transportability (Bareinboim & Pearl) | Medium | Selection diagrams not implemented; requires multi-population SCM framework |
| L8 | z-identifiability (Bareinboim & Pearl) | Medium | Surrogate experiments framework not yet implemented |

## Resolution Plan

1. Transportability and z-identifiability are medium priority for L8 completeness
2. Time-series and causal RL are L9 frontiers — documentation sufficient per SKILL.md requirements

