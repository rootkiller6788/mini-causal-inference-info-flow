# Information Flow in Networks

## From Causality to Information Dynamics

Information flow quantifies how knowledge about one variable reduces uncertainty about another, accounting for direction and dynamics.

## Transfer Entropy (Schreiber, 2000)
TE(X->Y) = I(Y_{t+1}; X_t | Y_t) — the directed information transfer from X to Y.

## Granger Causality (Granger, 1969)
X Granger-causes Y if past X improves prediction of Y beyond past Y alone.

## Causal Graphs
- d-separation: graphical criterion for conditional independence
- Do-calculus: intervention-based causal effects
- PC algorithm: constraint-based causal discovery

## Partial Information Decomposition
Decomposes I(T;S1,S2) into unique, redundant, and synergistic components.

