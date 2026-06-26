# Algorithms — PC & FCI Causal Discovery

## 1. Correlation Matrix

**Function:** `cdf_citest_corr_matrix(ds, corr)`

Computes p×p Pearson correlation matrix from N×p dataset.

**Algorithm:** Standard two-pass: means → standard deviations → covariances → correlations.

**Complexity:** O(N·p²)

## 2. Partial Correlation

**Function:** `cdf_citest_partial_corr(ds, x, y, Z, nZ)`

Computes ρ_{XY|Z} via precision matrix inversion.

**Algorithm:**
1. Build (2+|Z|) × (2+|Z|) correlation submatrix for {X,Y}∪Z
2. Invert via Gaussian elimination with partial pivoting
3. ρ_{XY|Z} = -Ω_{01} / √(Ω_{00}·Ω_{11})

**Complexity:** O(m³) where m = 2 + |Z|

## 3. Fisher Z-Test

**Function:** `cdf_citest_fisher_z(ds, x, y, Z, nZ, alpha)`

**Algorithm:**
1. Compute ρ = partial correlation
2. z_fisher = 0.5·ln((1+ρ)/(1-ρ))
3. se = 1/√(N-|Z|-3)
4. z_stat = z_fisher / se ∼ N(0,1) under H₀
5. p = 2·(1-Φ(|z_stat|))

**Complexity:** O((2+|Z|)³)

## 4. d-Separation

**Function:** `cdf_graph_d_separated(g, x, y, Z, nZ)`

Checks if X and Y are d-separated by Z.

**Algorithm (moralized graph BFS):**
1. Build ancestral graph of {X,Y}∪Z (keep ancestors)
2. Moralize: connect parents of common children, make edges undirected
3. Remove Z nodes
4. BFS from X avoiding Z nodes — if Y reachable, NOT d-separated

**Complexity:** O(p³) for moralization, O(p²) for BFS

## 5. PC Skeleton (Adjacency Search)

**Function:** `cdf_graph_skeleton_pc(g, ds, config)`

**Algorithm:**
```
For ℓ = 0, 1, 2, ... (conditioning set size):
  For each adjacent pair (i, j):
    Adj_i = neighbors(i) \ {j}
    For each subset S ⊆ Adj_i, |S| = ℓ:
      Test i ⊥ j | S
      If independent: remove edge i—j, record SepSet(i,j) = S
```

**Optimizations:**
- Stop early when no edges removed at current ℓ
- Limit subset enumeration (max 1000 subsets per pair)
- Only test subsets from shared adjacency

**Complexity:** O(p·d^{ℓ_max}·m³) where d = max degree, m = ℓ+2

## 6. V-Structure Detection

**Function:** `cdf_graph_find_vstructures(g)`

**Algorithm:**
For each unshielded triple (u—v—w) with u,w not adjacent:
  If v ∉ SepSet(u,w): orient u → v ← w

**Complexity:** O(p³)

## 7. Meek Orientation Rules (R1-R4)

**Functions:** `cdf_orient_rule_r1` through `cdf_orient_rule_r4`

**Algorithm:** Fixed-point iteration — apply all rules repeatedly until no more orientations possible.

**R1:** For each (a→b—c) with a,c not adj: orient b→c
**R2:** For each (a→b→c) with a—c: orient a→c
**R3:** For each (a—b→c, a—d→c) with b,d not adj and a—c: orient a→c
**R4:** Discriminating path heuristic

**Complexity:** O(p³) per iteration, typically converges in < 10 iterations

## 8. FCI Orientation Rules (R5-R10)

**Functions:** `cdf_orient_fci_rules_5_7`, `cdf_orient_fci_rules_8_10`

Additional rules for PAG orientation handling partially directed edges (∘→).

**R5:** a∘→b—c, a,c not adj → b→c
**R6:** a∘→b∘→c, a—c → b→c
**R7:** a∘→b→c, a—c → a→c
**R8:** a→b→c, a∘→c → a→c
**R9:** a∘→c, a∘→b→c → a→c
**R10:** a∘→c, a→b←c → a→c

**Complexity:** O(p³) per iteration

## 9. Possible-D-SEP

**Function:** `cdf_graph_possible_d_sep(g, v, max_len, pdsep)`

Finds nodes reachable from v via possibly-directed paths.

**Algorithm:** BFS along ∘→, →, ∘—∘, and — edges, limited to max_len steps.

**Complexity:** O(p·d^{max_len})

## 10. PC Algorithm (Full)

**Function:** `cdf_pc_run(ds, config)`

**Algorithm:**
1. Init complete undirected graph
2. Skeleton phase (adjacency search with CI tests)
3. V-structure detection
4. Meek rules (R1-R4)

**Output:** CPDAG

**Complexity:** O(p·d^{ℓ_max}·m³) for skeleton + O(p³) for v-structures

## 11. FCI Algorithm (Full)

**Function:** `cdf_fci_run(ds, config)`

**Algorithm:**
1. Init complete undirected graph
2. PC-like initial skeleton
3. V-structure detection
4. Convert undirected to ∘—∘ (PAG representation)
5. PDS phase (CI tests with PDS conditioning sets)
6. FCI orientation rules (R1-R10)

**Output:** PAG

## 12. RFCI

**Function:** `cdf_fci_rfci_run(ds, config)`

Simplified FCI: runs PC skeleton + v-structures, converts to PAG, applies FCI orientation rules. Skips expensive PDS phase.

**Complexity:** Same as PC, much faster than full FCI.

## Complexity Summary

| Algorithm | Complexity | Notes |
|-----------|-----------|-------|
| Correlation matrix | O(N·p²) | Two-pass |
| Partial correlation | O(m³) | m = 2 + \|Z\| |
| Fisher z-test | O(m³) | Dominated by partial corr |
| d-separation | O(p³) | Moralized graph |
| PC skeleton | O(p·d^ℓ·m³) | Exponential in ℓ |
| V-structures | O(p³) | Triple enumeration |
| Meek rules (R1-R4) | O(p³·iter) | ~5 iterations typical |
| FCI rules (R1-R10) | O(p³·iter) | ~10 iterations typical |
| PDS computation | O(p·d^L) | BFS limited by L |
| PC algorithm | O(p·d^ℓ·m³ + p³) | Dominated by skeleton |
| FCI algorithm | O(PC + p·d^L + p³) | Plus PDS phase |
| RFCI algorithm | O(PC complexity) | Skips PDS phase |