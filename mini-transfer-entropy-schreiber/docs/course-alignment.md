# Course Alignment — Transfer Entropy Schreiber

## Nine-School Curriculum Mapping

### MIT
- **6.050J/2.110J Information, Entropy, and Computation**
  - Unit 3: Information theory — entropy, mutual information ✓
  - Unit 8: Communication channels — directed information ✓
- **6.841 Advanced Complexity Theory**
  - Information-theoretic lower bounds ✓
  - Kolmogorov complexity and incompressibility method ✓

### Stanford
- **CS229 Machine Learning**
  - Feature selection via mutual information ✓
  - Independent Component Analysis (ICA) ✓
- **CS358 Circuit Complexity** (Indirect)
  - Information bottleneck in circuit lower bounds

### Berkeley
- **CS281A Statistical Learning Theory**
  - Graphical models and conditional independence ✓
  - Chow-Liu trees (MI-based structure learning) ✓
- **CS294 Deep Unsupervised Learning**
  - InfoMax principle ✓
  - Mutual information neural estimation (MINE) ✓

### CMU
- **36-708 Statistical Methods for Machine Learning**
  - Causal discovery algorithms ✓
  - Information-theoretic feature selection ✓
- **15-859 Information Theory in CS**
  - Directed information and feedback capacity ✓
  - Transfer entropy as directed information ✓

### Princeton
- **COS 551 Advanced Complexity Theory**
  - Algorithmic information theory ✓
  - Kolmogorov complexity as universal similarity metric ✓
- **COS 597C Causal Inference**
  - Granger causality and TE relationship ✓

### Caltech
- **CS 155 Machine Learning & Data Mining**
  - Information bottleneck method ✓
  - Spectral clustering via normalized cut (MI-based) ✓
- **CS 159 Advanced Topics in ML**
  - Deep InfoMax (Hjelm et al. 2019) ✓

### Cambridge
- **Part II Information Theory**
  - Channel capacity and coding theorems ✓
  - Rate-distortion theory ✓
  - Directed information (Marko 1973, Massey 1990) ✓
- **Part III Advanced Information Theory**
  - Network information theory ✓
  - Common randomness and secret key agreement ✓

### Oxford
- **Advanced Machine Learning**
  - Bayesian networks and d-separation ✓
  - Causal inference: PC algorithm, FCI algorithm ✓
  - Conditional independence testing ✓
- **Quantum Information Theory** (Indirect)
  - Quantum mutual information ✓
  - Squashed entanglement ✓

### ETH
- **263-5210 Advanced Machine Learning**
  - Information-theoretic causal discovery ✓
  - PCMCI algorithm (Runge et al. 2019) ✓
- **263-4650 Advanced Complexity Theory**
  - Information complexity in communication ✓

## Key Textbook References

| Topic | Source | Relevance |
|-------|--------|-----------|
| Transfer Entropy definition | Schreiber (2000) PRL | Core |
| Effective TE | Marschinski & Kantz (2002) EPJB | Core |
| KSG estimator | Kraskov et al. (2004) PRE | Core |
| Symbolic TE | Staniek & Lehnertz (2008) PRL | Core |
| Information dynamics | Lizier (2012) | AIS, PI |
| Directed information | Massey (1990) ISIT | Theory |
| PCMCI | Runge et al. (2019) Sci. Adv. | Application |
| Elements of IT | Cover & Thomas (2006) | Foundations |

## Module Coverage vs Course Expectations

| Course | Topics Covered | Topics Missing |
|--------|---------------|----------------|
| MIT 6.050 | Entropy, MI, TE ✓ | Channel coding |
| Stanford CS229 | MI feature selection ✓ | Deep InfoMax |
| Berkeley CS281A | Conditional independence ✓ | Structure learning |
| CMU 36-708 | Causal discovery ✓ | PC/FCI algorithms |
| Princeton COS 551 | Algorithmic IT ✓ | Resource-bounded KC |
| Caltech CS 155 | Info bottleneck ✓ | Deep IB |
| Cambridge Part II | Directed info ✓ | Continuous-time TE |
| Oxford ML | Causal inference ✓ | do-calculus |
| ETH 263-5210 | PCMCI ✓ | Multi-variate extensions |
