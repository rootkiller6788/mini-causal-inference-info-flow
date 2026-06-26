#ifndef CDF_SCORES_H
#define CDF_SCORES_H

/*
 * cdf_scores.h — Score-Based Causal Discovery
 *
 * Constraint-based methods (PC, FCI) use CI tests to learn graphs.
 * Score-based methods evaluate graph structures by a goodness-of-fit
 * score penalized for complexity, searching over the space of DAGs.
 *
 * Scores implemented:
 *   1. BIC (Bayesian Information Criterion) for linear-Gaussian models
 *   2. BGe (Bayesian Gaussian equivalent) marginal likelihood
 *   3. AIC (Akaike Information Criterion)
 *   4. Graph search via greedy and tabu hill-climbing
 *
 * References:
 *   - Chickering (2002), "Optimal Structure Identification With Greedy Search"
 *   - Geiger & Heckerman (2002), "Parameter priors for linear-Gaussian DAGs"
 *   - Schwarz (1978), "Estimating the dimension of a model" (BIC)
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "cdf_core.h"

/* ── Score Types ───────────────────────────────────────────────────── */

typedef enum {
    CDF_SCORE_BIC = 0,   /* Bayesian Information Criterion */
    CDF_SCORE_BGE = 1,   /* Bayesian Gaussian equivalent score */
    CDF_SCORE_AIC = 2    /* Akaike Information Criterion */
} CdfScoreType;

/* ── BIC Score ──────────────────────────────────────────────────────── */

/**
 * BIC score for a linear-Gaussian DAG.
 *
 * For each node v with parent set Pa(v): BIC(v,Pa(v)) = -N·ln(σ̂²) - (|Pa|+1)·ln(N)
 *
 * σ̂² is the residual variance of regressing X_v on its parents.
 * Total BIC = Σ_v BIC(v, Pa(v)) — decomposable over nodes.
 *
 * @param ds  dataset
 * @param g   DAG (directed edges = parent sets)
 * @return    total BIC score (higher is better)
 */
double cdf_score_bic(const CdfDataset *ds, const CdfGraph *g);

/**
 * Local BIC score for node v given parents.
 * Used for ΔBIC computation in greedy search.
 */
double cdf_score_bic_local(const CdfDataset *ds, int v,
                            const int *parents, int n_parents);

/* ── BGe Score ──────────────────────────────────────────────────────── */

/**
 * BGe score (Bayesian Gaussian equivalent).
 *
 * Uses normal-Wishart prior: integrates over parameters analytically.
 * P(Data|G) = ∏_v P(X_v | Pa(v), Data)
 *
 * @param ds       dataset
 * @param g        DAG
 * @param alpha_w  Wishart df (default: ds->p + 2)
 * @param alpha_mu prior precision on mean (default: 1.0)
 * @return         log marginal likelihood (BGe score)
 */
double cdf_score_bge(const CdfDataset *ds, const CdfGraph *g,
                      double alpha_w, double alpha_mu);
double cdf_score_bge_local(const CdfDataset *ds, int v,
                            const int *parents, int n_parents,
                            double alpha_w, double alpha_mu);

/* ── AIC Score ──────────────────────────────────────────────────────── */

double cdf_score_aic(const CdfDataset *ds, const CdfGraph *g);
double cdf_score_aic_local(const CdfDataset *ds, int v,
                            const int *parents, int n_parents);

/* ── Linear Regression ──────────────────────────────────────────────── */

/**
 * Fit y = Xβ + ε via normal equations XᵀX β̂ = Xᵀy.
 * Uses Gaussian elimination with partial pivoting.
 *
 * @param y      [N] response vector
 * @param X      [N × k] design matrix, column-major
 * @param N      samples
 * @param k      predictors
 * @param beta   [k] output coefficients
 * @param sigma2 output residual variance
 * @return       true on success
 */
bool cdf_score_linear_regression(const double *y, const double *X,
                                  int N, int k,
                                  double *beta, double *sigma2);

/* ── Greedy DAG Search ──────────────────────────────────────────────── */

/**
 * Greedy hill-climbing over DAG space (GES-style).
 *
 * Operations at each iteration:
 *   - Add edge u→v if no cycle and u,v not adjacent
 *   - Delete edge u→v
 *   - Reverse edge u→v to v→u if no cycle
 *
 * Takes the operation with maximum positive Δscore.
 * Stops when no improving operation exists (local optimum).
 *
 * @param ds        dataset
 * @param max_iter  max iterations (0 = unlimited)
 * @param score_type BIC, BGE, or AIC
 * @param verbose   print progress
 * @return          best DAG (caller frees)
 */
CdfGraph* cdf_score_greedy_search(const CdfDataset *ds, int max_iter,
                                   CdfScoreType score_type, bool verbose);

/**
 * Tabu search: maintains a tabu list of recently visited edge changes
 * to escape local optima and explore higher-scoring regions.
 */
CdfGraph* cdf_score_tabu_search(const CdfDataset *ds, int max_iter,
                                 int tabu_len, CdfScoreType score_type,
                                 bool verbose);

/* ── Score Delta ────────────────────────────────────────────────────── */

/**
 * ΔScore = Score(G') - Score(G) for local edge operation.
 *
 * Score decomposability means only nodes with changed parent sets
 * need recomputation. For add u→v: only v changes.
 * For delete u→v: only v changes.
 * For reverse u→v to v→u: both u and v change.
 *
 * @param ds  dataset
 * @param g   current DAG
 * @param u   first node
 * @param v   second node
 * @param op  'a' add, 'd' delete, 'r' reverse
 * @param type score type
 * @return    score difference (positive = improvement)
 */
double cdf_score_delta(const CdfDataset *ds, const CdfGraph *g,
                        int u, int v, char op, CdfScoreType type);

/* ── Bootstrap Stability ────────────────────────────────────────────── */

/**
 * Bootstrap-based edge stability.
 * Resample B times, run greedy search, count edge frequency.
 * stability[u*p+v] ∈ [0,1] = frequency of edge (any direction).
 *
 * @param ds        dataset
 * @param B         bootstrap replicates
 * @param stability [p*p] output matrix (symmetric)
 * @param score_type score to optimize
 * @param verbose   progress output
 */
void cdf_score_bootstrap_stability(const CdfDataset *ds, int B,
                                    double *stability,
                                    CdfScoreType score_type,
                                    bool verbose);

/* ── Utilities ──────────────────────────────────────────────────────── */

double cdf_score_compare(const CdfDataset *ds, const CdfGraph *g1,
                          const CdfGraph *g2, CdfScoreType type);
void cdf_score_print_report(const CdfDataset *ds, const CdfGraph *g,
                             CdfScoreType type);
int cdf_score_num_params(const CdfGraph *g);

#endif /* CDF_SCORES_H */
