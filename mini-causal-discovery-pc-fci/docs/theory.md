# Causal Discovery: PC & FCI Algorithms — Theory

## Causal Graphical Models

A **causal DAG** (Directed Acyclic Graph) G = (V, E) represents
causal relationships among variables V. An edge X → Y means X is a
direct cause of Y.

### Markov Property

The joint distribution P factorizes according to G:
```
P(X₁,...,Xₚ) = Πᵢ P(Xᵢ | parents(Xᵢ))
```

### Faithfulness

The population distribution P is faithful to G if all conditional
independencies in P are entailed by d-separation in G (and vice versa).

## d-Separation

A path π between X and Y is **d-separated** (blocked) by set Z if:
1. π contains a chain i → m → j or fork i ← m → j with m ∈ Z, OR
2. π contains a collider i → m ← j with m ∉ Z and no descendant of m in Z.

X and Y are d-separated by Z if Z blocks every path between them.

## Markov Equivalence

Two DAGs are Markov equivalent if they encode the same set of
conditional independencies. The equivalence class is represented by
a **CPDAG** (Completed Partially Directed Acyclic Graph).

CPDAG has:
- Directed edges: common to all DAGs in the class
- Undirected edges: orientation differs across the class

## PC Algorithm (Spirtes & Glymour, 1991)

### Assumptions
1. **Causal Markov Condition**: d-separation ⇒ CI
2. **Faithfulness**: CI ⇒ d-separation
3. **Causal Sufficiency**: No latent confounders

### Algorithm Phases

**Phase 1 — Skeleton**: Start with complete undirected graph. For ℓ = 0,1,2,...:
- For each adjacent pair (X,Y), test X ⊥ Y | Z for all Z ⊆ Adj(X)\{Y} with |Z| = ℓ
- If independent, remove edge X—Y and record Z as SepSet(X,Y)

**Phase 2 — V-Structures**: For each unshielded triple X—Z—Y (X,Y not adj):
- If Z ∉ SepSet(X,Y), orient X → Z ← Y (collider at Z)

**Phase 3 — Meek Rules**: Apply orientation rules R1-R4 exhaustively:
- R1: Avoid new v-structures
- R2: Avoid cycles (transitive closure)
- R3: Complex v-structure avoidance
- R4: Discriminating paths

## FCI Algorithm (Spirtes et al., 1999; Zhang, 2008)

FCI relaxes the causal sufficiency assumption, allowing latent
confounders and selection bias.

### Key Differences from PC
1. Outputs a **PAG** (Partial Ancestral Graph) instead of CPDAG
2. PAG represents the equivalence class of **MAGs** (Maximal Ancestral Graphs)
3. Introduces **Possible-D-SEP** sets for additional CI tests
4. Uses **10 orientation rules** (R1-R10) instead of 4

### Edge Types in PAG
| Symbol | Meaning |
|--------|---------|
| A → B  | A is a cause of B (direct or indirect) |
| A ↔ B  | Latent confounder (no causal relation, but correlated) |
| A ∘→ B | B is NOT a cause of A |
| A ∘—∘ B | Uncertain relationship |

### FCI Phases

1. Initial skeleton (like PC)
2. V-structure detection
3. Possible-D-SEP computation and CI testing
4. Apply R1-R10 orientation rules

### RFCI (Really Fast Causal Inference)
A simplified variant that skips the expensive PDS phase — uses only
adjacency sets. Faster but potentially less accurate.

## Orientation Rules

### PC Rules (R1-R4)

- **R0**: Orient v-structures
- **R1**: a→b—c, a,c not adj → b→c
- **R2**: a→b→c, a—c → a→c
- **R3**: a—b→c, a—d→c, b,d not adj, a—c → a→c
- **R4**: Discriminating path → definite orientation

### FCI Additional Rules (R5-R10)

- **R5**: a∘→b—c, a,c not adj → b→c
- **R6**: a∘→b∘→c, a—c → b→c
- **R7**: a∘→b→c, a—c → a→c
- **R8**: a→b→c, a∘→c → a→c
- **R9**: a∘→c, a∘→b→c → a→c
- **R10**: a∘→c, a→b←c → a→c

## Statistical Tests

### Fisher's Z-Test (Gaussian Data)
Tests ρ_{XY|Z} = 0 via:
```
z = 0.5·ln((1+ρ)/(1-ρ)) · √(N - |Z| - 3)  ∼ N(0,1)
```

### G² Test (Discrete Data)
Tests independence via contingency table chi-square.

## References

- Spirtes, P., Glymour, C., & Scheines, R. (2000). *Causation, Prediction, and Search*. MIT Press.
- Pearl, J. (2009). *Causality*, 2nd ed. Cambridge.
- Zhang, J. (2008). On the completeness of orientation rules for causal discovery.
- Meek, C. (1995). Causal inference and causal explanation with background knowledge.
- Colombo, D. et al. (2012). Learning high-dimensional DAGs with latent and selection variables.