# Knowledge Graph — Transfer Entropy Schreiber

## L1: Definitions (Complete)

- Transfer Entropy: TE_{Y->X} = I(X_f; Y_p | X_p)
- Mutual Information: I(X;Y) = H(X) + H(Y) - H(X,Y)
- Conditional Transfer Entropy: TE_{Y->X|Z}
- Effective Transfer Entropy: ETE = TE - bias
- Symbolic Transfer Entropy: permutation-encoded TE
- Active Information Storage: AIS(X) = I(X_{n+1}; X_n^{(k)})
- Predictive Information: PI = I(X_future; X_past)
- Information dynamics: Lizier (2012) framework

## L2: Core Concepts (Complete)

- TE asymmetry: TE(X->Y) != TE(Y->X)
- Directed information: Massey (1990)
- Conditional independence testing
- Surrogate data testing for significance
- Embedding: delay, dimension (k, l)
- Bin width selection (Freedman-Diaconis)
- Information flow in coupled systems

## L3: Mathematical Structures (Complete)

- 3D joint probability tables: p(x_f, x_p, y_p)
- Factorial number system / Lehmer code for permutations
- k-NN distance metrics: Euclidean, Chebyshev, Manhattan, Minkowski, cosine
- Kernel Density Estimation (KDE) for continuous entropy
- Digamma function (psi) asymptotic expansion
- KL divergence, Jensen-Shannon divergence, Renyi entropy, Tsallis entropy
- Gibbs inequality

## L4: Fundamental Laws (Complete)

- Schreiber (2000): TE definition and binning estimator
- Kaiser-Schreiber (2002): Continuous TE
- Marschinski-Kantz (2002): Effective TE with bias correction
- Kraskov-Stogbauer-Grassberger (2004): KSG k-NN MI estimator
- Staniek-Lehnertz (2008): Symbolic TE
- Lizier et al. (2012): AIS and information dynamics
- Fadlallah et al. (2013): Weighted permutation entropy
- Transfer entropy non-negativity (KL divergence form)
- TE zero condition (conditional independence)

## L5: Algorithms/Methods (Complete)

1. Schreiber binning TE estimator
2. KSG k-NN TE estimator (with Chebyshev distance)
3. Symbolic/permutation TE estimator
4. Effective TE (surrogate subtraction)
5. Conditional/multivariate TE
6. Multi-scale TE (coarse-graining)
7. Time-varying TE (sliding window)
8. Bootstrap confidence intervals
9. Jackknife bias correction
10. Permutation test for significance
11. KS test (two-sample)
12. Mann-Whitney U test
13. BDS test for nonlinearity
14. Phase randomization surrogates
15. IAAFT iterative surrogates
16. Block bootstrap / stationary bootstrap
17. Granger F-test (linear comparison)
18. Partial Information Decomposition (Williams-Beer 2010)
19. Holm-Bonferroni / Benjamini-Hochberg correction
20. PageRank on causal networks
21. Label propagation community detection

## L6: Canonical Problems (Complete)

1. Coupled Henon map: TE detects asymmetric coupling
2. Coupled AR(1): TE proportional to coupling strength c
3. Logistic map: information flow during period-doubling
4. Rossler system: phase synchronization detection
5. ENSO teleconnection: lagged climate causality
6. Neural spike trains: synaptic connectivity from recordings

## L7: Applications (Complete)

- L7.1: Neural information transfer (NASA/Ames pipeline)
  - Poisson spike train generation with synaptic coupling
  - ISI statistics, firing rate analysis
  - TE-based functional connectivity
  - Vicente et al. (2011) J. Comput. Neurosci.
- L7.2: Climate teleconnection analysis
  - ENSO-like index with ~2-7 year periodicity
  - Regional temperature response with lagged coupling
  - Bootstrap confidence for link detection
  - Runge et al. (2019) Science Advances
- L7.3: Econophysics / financial TE
  - Stock market index causality

## L8: Advanced Topics (Complete)

- Conditional/multivariate TE with confounder control
- Partial Information Decomposition (redundancy, unique, synergy)
- Theiler window correction for temporal correlations
- Multi-scale coarse-graining for long-range detection
- KSG estimator with adaptive k selection
- Jackknife and bootstrap variance estimation
- Causal network analysis (PageRank, label propagation)
- Granger causality (linear alternative) comparison
- False Discovery Rate (Benjamini-Hochberg) control
- Transfer entropy spectrum over embedding dimensions
- Amplitude-aware permutation entropy (AAPE)

## L9: Research Frontiers (Partial)

- Deep transfer entropy (neural estimation)
- Causal discovery (PCMCI framework)
- Quantum transfer entropy
- Reservoir computing TE
- Meta-causal learning
- Continuous-time TE estimation
