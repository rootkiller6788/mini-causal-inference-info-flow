# Do-Calculus & Intervention — Theory Reference

## 1. Structural Causal Models (SCM)

A Structural Causal Model (Pearl, 2000) is a triple M = ⟨U, V, F, P(u)⟩:
- **U**: Exogenous (unobserved) variables
- **V**: Endogenous (observed) variables
- **F**: Structural equations V_i = f_i(PA_i, U_i)
- **P(u)**: Joint distribution over U

### Causal Diagram
Each SCM induces a directed graph G where V_i → V_j if V_i ∈ PA_j.

## 2. Interventions and the do-Operator

### Intervention do(X=x)
The operation do(X=x) replaces the structural equation for X with the constant x:
```
f_X replaced by: X = x
```

### Truncated Factorization (g-formula)
In the mutilated graph, the joint distribution factorizes as:
```
P(v | do(x)) = Π_{i: V_i ∉ X} P(v_i | pa_i)  evaluated at X = x
```

## 3. Three Rules of do-Calculus

### Rule 1 (Insertion/deletion of observations)
```
P(y | do(x), z, w) = P(y | do(x), w)
if Y ⟂ Z | X, W in G_{X̄}
```

### Rule 2 (Action/observation exchange)
```
P(y | do(x), do(z), w) = P(y | do(x), z, w)
if Y ⟂ Z | X, W in G_{X̄Z̲}
```

### Rule 3 (Insertion/deletion of actions)
```
P(y | do(x), do(z), w) = P(y | do(x), w)
if Y ⟂ Z | X, W in G_{X̄Z(W)̄}
```

## 4. Confounding Control

### Back-Door Criterion
A set Z satisfies the back-door criterion relative to (X, Y) if:
1. No node in Z is a descendant of X
2. Z blocks every back-door path from X to Y

### Back-Door Adjustment
```
P(y | do(x)) = Σ_z P(y | x, z) P(z)
```

### Front-Door Criterion
M satisfies the front-door criterion if:
1. M intercepts all directed paths from X to Y
2. No back-door path from X to M
3. All back-door paths from M to Y blocked by X

### Front-Door Adjustment
```
P(y | do(x)) = Σ_m P(m | x) Σ_{x'} P(y | x', m) P(x')
```

## 5. Counterfactuals

A counterfactual Y_x(u) = y is the value Y would take in unit u, had X been x.

### Three-Step Procedure
1. **Abduction**: Update P(u) given evidence E=e
2. **Action**: Modify model with do(X=x)
3. **Prediction**: Compute P(Y=y) in modified model

### Probability of Necessity and Sufficiency
- **PN**: P(Y_{X=0}=0 | X=1, Y=1) — was X necessary?
- **PS**: P(Y_{X=1}=1 | X=0, Y=0) — was X sufficient?
- **PNS**: P(Y_{X=1}=1, Y_{X=0}=0) — necessity AND sufficiency

## References
- Pearl, J. (2009). Causality: Models, Reasoning, and Inference (2nd ed.). Cambridge.
- Pearl, J. (1995). Causal diagrams for empirical research. Biometrika, 82(4):669-710.
- Pearl, J. & Mackenzie, D. (2018). The Book of Why. Basic Books.
- Spirtes, P., Glymour, C., & Scheines, R. (2000). Causation, Prediction, and Search. MIT Press.
