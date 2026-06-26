# Knowledge Graph -- Causal Discovery PC/FCI

## L1: Definitions
- DAG: Directed Acyclic Graph (true causal structure)
- CPDAG: Completed Partially Directed Acyclic Graph (PC output)
- MAG: Maximal Ancestral Graph (with latent variables)
- PAG: Partial Ancestral Graph (FCI output)
- d-separation: graphical criterion for conditional independence
- SepSet: separating set d-separating two nodes
- Markov equivalence: same d-separation relations
- Causal effect: E[Y | do(X=x)]

## L2: Core Concepts
- Constraint-based causal discovery: learn edges from CI tests
- Causal Markov Condition: d-sep ⇒ CI in the population
- Faithfulness: CI ⇒ d-sep (no accidental cancellations)
- Causal sufficiency: no unobserved common causes (PC assumption)
- Latent confounders: unobserved variables causing observed ones
- Selection bias: non-random sampling affecting dependence patterns

## L3: Math Structures
- Adjacency matrix: edges[u*p+v] = edge type
- d-separation moralization: connect parents, make undirected, remove Z, check connectivity
- Partial correlation via precision matrix: rho_{XY|Z} = -Omega_01 / sqrt(Omega_00 * Omega_11)
- Fisher z-transform: z = 0.5 * ln((1+rho)/(1-rho)) ~ N(0, 1/(N-|Z|-3))
- Kernel Gram matrix: K_ij = k(x_i, x_j) in RKHS
- Regression design matrix: X [N x (k+1)] with intercept

## L4: Fundamental Laws
- d-Separation criterion (Pearl, 1988): X ⊥_d Y | Z iff Z blocks all paths
- PC algorithm correctness (SGS, 2000, Theorem 5.1): Under CMC + faithfulness + sufficiency, PC recovers the true CPDAG as N→∞
- Meek's rules soundness (1995): R1-R4 never introduce orientations contradicting the Markov equivalence class
- Verma-Pearl characterization (1990): Two DAGs Markov equivalent iff same skeleton + same unshielded colliders
- Zhang FCI completeness (2008, Theorem 3): R1-R10 maximally orient all orientable edges in a PAG
- IDA consistency (Maathuis et al., 2009, Theorem 1): The multiset of possible effects contains the true effect with probability → 1

## L5: Algorithms/Methods
- PC algorithm: skeleton → v-structures → Meek rules
- FCI algorithm: skeleton → v-structures → PDS → FCI rules R1-R10
- RFCI: FCI without expensive PDS phase
- Skeleton adjacency search: conditioning set size ℓ = 0, 1, 2, ...
- V-structure detection: unshielded colliders u→v←w when v ∉ SepSet(u,w)
- Meek rules R1-R4: no new v-structures, acyclicity, discriminating paths
- FCI rules R5-R10: latent variable structures in PAGs
- Conservative PC (Ramsey et al.): only orient if all same-size CI tests agree
- Greedy DAG search (score-based): iteratively apply best edge operation
- Tabu search: greedy + tabu list to escape local optima
- Bootstrap edge stability: resample B times, count edge frequency
- IDA causal effect estimation: enumerate parent sets, regression per DAG
- Kernel HSIC: Tr(Kx H Ky H) / n^2 test statistic
- Kernel KCI: Tr(Kx_z Ky_z) / n test statistic with conditional centering

## L6: Canonical Problems
- Chain discovery: X1→X2→X3 → skeleton X1—X2—X3
- Collider identification: X1→X2←X3 → v-structure orientation
- Fork detection: X1←X2→X3 → d-separation X1 ⊥ X3 | X2
- Latent confounder: X1↔X2 in PAG (bidirected edge)
- Markov equivalence class: skeleton + colliders uniquely determine equivalence
- CPDAG recovery from Gaussian data

## L7: Applications
- Climate science: causal links among Temp, CO2, Solar, Ocean, Ice
- Supply chain: causal drivers from Demand, LeadTime, Inventory, Production, Delivery
- Gene regulatory networks: reverse-engineering from expression data
- fMRI brain connectivity: functional connectivity from BOLD signals
- Epidemiology: risk factor identification from observational studies

## L8: Advanced Topics
- IDA: Intervention-calculus when DAG is Absent
- Kernel CI tests: HSIC for non-linear dependence, KCI for conditional
- Bootstrap stability selection with FDR control (Meinshausen-Buhlmann)
- Score-based discovery: BIC, BGe, AIC scores with greedy/tabu search
- Linear-Gaussian SEM: parameter estimation via normal equations
- Back-door criterion: admissible adjustment sets for causal effect estimation
- Front-door criterion: mediation-based causal effect estimation
- Discriminating paths: resolving ambiguous orientations in PAGs

## L9: Research Frontiers
- Differentiable causal discovery (NOTEARS, DAG-GNN)
- Neural-network-based structure learning
- High-dimensional consistency of PC (Kalisch & Buhlmann, 2007)

