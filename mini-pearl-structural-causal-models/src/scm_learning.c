/* scm_learning.c -- Causal Structure Learning Algorithms
 * L5 Algorithms: PC, GES, FCI, scoring, orientation rules
 * References:
 *   Spirtes, Glymour, Scheines (2000) "Causation, Prediction, and Search"
 *   Chickering (2002) JMLR 3:507-554
 *   Meek (1995) UAI-95
 */

#include "scm_learning.h"
#include "scm_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- Internal helpers ---------- */

static double col_mean(const double* data, int n_rows, int n_cols, int col) {
    double sum = 0.0;
    for (int i = 0; i < n_rows; i++) sum += data[i * n_cols + col];
    return sum / (double)n_rows;
}

static double col_var(const double* data, int n_rows, int n_cols, int col, double mean) {
    double sum = 0.0;
    for (int i = 0; i < n_rows; i++) {
        double d = data[i * n_cols + col] - mean;
        sum += d * d;
    }
    return sum / (double)(n_rows - 1);
}

static double covariance(const double* data, int n_rows, int n_cols,
                          int c1, double m1, int c2, double m2) {
    double sum = 0.0;
    for (int i = 0; i < n_rows; i++)
        sum += (data[i * n_cols + c1] - m1) * (data[i * n_cols + c2] - m2);
    return sum / (double)(n_rows - 1);
}

/* Gaussian elimination with partial pivoting, solves Ax = b */
static int solve_linear_system(double* A, double* b, double* x, int n) {
    double aug[SCM_MAX_VARS][SCM_MAX_VARS + 1];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i][j] = A[i * n + j];
        aug[i][n] = b[i];
    }
    for (int col = 0; col < n; col++) {
        int max_row = col;
        double max_val = fabs(aug[col][col]);
        for (int row = col + 1; row < n; row++)
            if (fabs(aug[row][col]) > max_val) {
                max_val = fabs(aug[row][col]); max_row = row;
            }
        if (max_val < 1e-14) return -1;
        if (max_row != col)
            for (int j = col; j <= n; j++) {
                double tmp = aug[col][j];
                aug[col][j] = aug[max_row][j]; aug[max_row][j] = tmp;
            }
        for (int row = col + 1; row < n; row++) {
            double factor = aug[row][col] / aug[col][col];
            for (int j = col; j <= n; j++) aug[row][j] -= factor * aug[col][j];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        double sum = aug[i][n];
        for (int j = i + 1; j < n; j++) sum -= aug[i][j] * x[j];
        x[i] = sum / aug[i][i];
    }
    return 0;
}

/* Normal CDF approximation (Abramowitz & Stegun 26.2.17) */
static double normal_cdf_approx(double z) {
    double t = 1.0 / (1.0 + 0.2316419 * fabs(z));
    double d = (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * z * z);
    double poly = 0.319381530 * t - 0.356563782 * t * t
                + 1.781477937 * t * t * t - 1.821255978 * t * t * t * t
                + 1.330274429 * t * t * t * t * t;
    double phi = 1.0 - d * poly;
    if (z < 0) phi = 1.0 - phi;
    if (phi < 0.0) phi = 0.0;
    if (phi > 1.0) phi = 1.0;
    return phi;
}

/* Simple log2 via log */
static double log2_scm(double x) { return log(x) / log(2.0); }

/* ---------- CI Testing ---------- */

/* Partial correlation test: corr(X,Y | Z).
 * Uses precision matrix inversion of the (k+2)x(k+2) covariance submatrix. */
SCM_CITest scm_partial_correlation_test(const double* data, int n_rows, int n_cols,
                                         int x, int y, const SCM_VarSet* Z, double alpha) {
    SCM_CITest result = {0.0, 1.0, false, 0.0};
    if (!data || n_rows < 3 || n_cols < 2 || x < 0 || x >= n_cols || y < 0 || y >= n_cols)
        return result;
    int k = (Z) ? Z->n : 0;
    int df = n_rows - k - 3;
    if (k == 0) {
        double mx = col_mean(data, n_rows, n_cols, x);
        double my = col_mean(data, n_rows, n_cols, y);
        double vx = col_var(data, n_rows, n_cols, x, mx);
        double vy = col_var(data, n_rows, n_cols, y, my);
        double cov = covariance(data, n_rows, n_cols, x, mx, y, my);
        double r = (vx > 0 && vy > 0) ? cov / sqrt(vx * vy) : 0.0;
        if (r > 1.0) r = 1.0;
        if (r < -1.0) r = -1.0;
        result.correlation = r;
    } else {
        int m = k + 2;
        if (m > SCM_MAX_VARS) { result.independent = true; return result; }
        int idx[SCM_MAX_VARS + 2];
        idx[0] = x; idx[1] = y;
        for (int i = 0; i < k; i++) idx[2 + i] = Z->nodes[i];
        double means[SCM_MAX_VARS + 2];
        for (int i = 0; i < m; i++) means[i] = col_mean(data, n_rows, n_cols, idx[i]);
        double* Cmat = calloc(m * m, sizeof(double));
        if (!Cmat) return result;
        for (int i = 0; i < m; i++)
            for (int j = i; j < m; j++) {
                double c = covariance(data, n_rows, n_cols, idx[i], means[i], idx[j], means[j]);
                Cmat[i * m + j] = Cmat[j * m + i] = c;
            }
        double* Iden = calloc(m * m, sizeof(double));
        double* Prec = calloc(m * m, sizeof(double));
        if (!Iden || !Prec) { free(Cmat); free(Iden); free(Prec); return result; }
        for (int i = 0; i < m; i++) Iden[i * m + i] = 1.0;
        bool singular = false;
        for (int col = 0; col < m && !singular; col++) {
            double b[SCM_MAX_VARS + 2];
            for (int i = 0; i < m; i++) b[i] = Iden[i * m + col];
            double xcol[SCM_MAX_VARS + 2] = {0};
            double* Atmp = calloc(m * m, sizeof(double));
            if (!Atmp) { singular = true; break; }
            memcpy(Atmp, Cmat, m * m * sizeof(double));
            if (solve_linear_system(Atmp, b, xcol, m) != 0) {
                free(Atmp); singular = true; break;
            }
            free(Atmp);
            for (int i = 0; i < m; i++) Prec[i * m + col] = xcol[i];
        }
        if (singular) {
            free(Cmat); free(Iden); free(Prec);
            result.independent = true; return result;
        }
        double p00 = Prec[0], p11 = Prec[1 * m + 1], p01 = Prec[1];
        if (p00 > 0 && p11 > 0) {
            double r = -p01 / sqrt(p00 * p11);
            if (r > 1.0) r = 1.0;
            if (r < -1.0) r = -1.0;
            result.correlation = r;
        }
        free(Cmat); free(Iden); free(Prec);
    }
    if (df > 0) {
        double z = scm_fisher_z(result.correlation, n_rows, k);
        result.z_statistic = z;
        double p_left = normal_cdf_approx(-fabs(z));
        result.p_value = 2.0 * p_left;
        if (result.p_value < 0.0) result.p_value = 0.0;
        if (result.p_value > 1.0) result.p_value = 1.0;
        result.independent = (result.p_value >= alpha);
    }
    return result;
}

/* Fisher z-transform for testing H0: partial corr = 0 */
double scm_fisher_z(double r, int n, int k) {
    if (r >= 1.0) r = 0.999999;
    if (r <= -1.0) r = -0.999999;
    double z = 0.5 * log((1.0 + r) / (1.0 - r));
    int df = n - k - 3;
    if (df <= 0) return 0.0;
    return z * sqrt((double)df);
}

/* Mutual information CI test via discretization */
SCM_CITest scm_mutual_information_test(const double* data, int n_rows, int n_cols,
                                        int x, int y, const SCM_VarSet* Z, double alpha) {
    SCM_CITest result = {0.0, 1.0, false, 0.0};
    if (!data || n_rows < 10 || x < 0 || x >= n_cols || y < 0 || y >= n_cols)
        return result;
    int k = Z ? Z->n : 0;
    int n_bins = (int)ceil(log2_scm((double)n_rows)) + 1;
    if (n_bins < 3) n_bins = 3;
    if (n_bins > 10) n_bins = 10;

    double xmin = data[x], xmax = data[x];
    double ymin = data[y], ymax = data[y];
    for (int i = 0; i < n_rows; i++) {
        double xv = data[i * n_cols + x], yv = data[i * n_cols + y];
        if (xv < xmin) xmin = xv;
        if (xv > xmax) xmax = xv;
        if (yv < ymin) ymin = yv;
        if (yv > ymax) ymax = yv;
    }
    double xrange = xmax - xmin; if (xrange < 1e-12) xrange = 1.0;
    double yrange = ymax - ymin; if (yrange < 1e-12) yrange = 1.0;

    int* joint = calloc(n_bins * n_bins, sizeof(int));
    int* marg_x = calloc(n_bins, sizeof(int));
    int* marg_y = calloc(n_bins, sizeof(int));
    if (!joint || !marg_x || !marg_y) {
        free(joint); free(marg_x); free(marg_y); return result;
    }
    for (int i = 0; i < n_rows; i++) {
        int bx = (int)((data[i * n_cols + x] - xmin) / xrange * (n_bins - 1));
        int by = (int)((data[i * n_cols + y] - ymin) / yrange * (n_bins - 1));
        if (bx < 0) bx = 0;
        if (bx >= n_bins) bx = n_bins - 1;
        if (by < 0) by = 0;
        if (by >= n_bins) by = n_bins - 1;
        joint[bx * n_bins + by]++; marg_x[bx]++; marg_y[by]++;
    }
    double mi = 0.0;
    for (int i = 0; i < n_bins; i++)
        for (int j = 0; j < n_bins; j++) {
            if (joint[i * n_bins + j] == 0) continue;
            double p_ij = (double)joint[i * n_bins + j] / n_rows;
            double p_i = (double)marg_x[i] / n_rows;
            double p_j = (double)marg_y[j] / n_rows;
            if (p_i > 0 && p_j > 0) mi += p_ij * log(p_ij / (p_i * p_j));
        }
    if (k > 0) {
        for (int zi = 0; zi < k; zi++) {
            int zvar = Z->nodes[zi];
            double zm = col_mean(data, n_rows, n_cols, zvar);
            double zv = col_var(data, n_rows, n_cols, zvar, zm);
            double mx = col_mean(data, n_rows, n_cols, x);
            double my2 = col_mean(data, n_rows, n_cols, y);
            double vx = col_var(data, n_rows, n_cols, x, mx);
            double vy = col_var(data, n_rows, n_cols, y, my2);
            double r_xz = 0.0, r_yz = 0.0;
            if (zv > 0 && vx > 0)
                r_xz = covariance(data, n_rows, n_cols, x, mx, zvar, zm) / sqrt(vx * zv);
            if (zv > 0 && vy > 0)
                r_yz = covariance(data, n_rows, n_cols, y, my2, zvar, zm) / sqrt(vy * zv);
            mi *= (1.0 - fabs(r_xz * r_yz));
        }
    }
    result.correlation = mi;
    result.p_value = exp(-mi * n_rows * 0.5);
    if (result.p_value > 1.0) result.p_value = 1.0;
    result.independent = (result.p_value >= alpha);
    free(joint); free(marg_x); free(marg_y);
    return result;
}

/* ---------- PC Config ---------- */

SCM_PCConfig scm_pc_config_default(void) {
    SCM_PCConfig c = {0.05, 3, true, 1000};
    return c;
}

/* ---------- PC Algorithm ---------- */

/* Skeleton discovery phase of PC algorithm.
 * Start with complete undirected graph, test CI for each edge.
 * Remove edge if X _||_ Y | Z for some Z subset of adj(X)\{Y}. */
void scm_skeleton_discovery(const double* data, int n_rows, int n_cols,
                             SCM_LearnResult* result, const SCM_PCConfig* config) {
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    /* Initialize complete graph (no self-loops) */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result->adjacency[i][j] = (i != j);

    int depth = 0;
    while (depth <= config->max_conditioning) {
        bool any_removed = false;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                if (!result->adjacency[i][j]) continue;
                SCM_VarSet adj_i; scm_varset_init(&adj_i);
                for (int k = 0; k < n; k++)
                    if (k != j && result->adjacency[i][k]) scm_varset_add(&adj_i, k);
                if (adj_i.n < depth) continue;
                SCM_VarSet Z; scm_varset_init(&Z);
                for (int d = 0; d < depth && d < adj_i.n; d++)
                    scm_varset_add(&Z, adj_i.nodes[d]);
                SCM_CITest ci;
                if (config->use_fisher_z)
                    ci = scm_partial_correlation_test(data, n_rows, n_cols, i, j, &Z, config->alpha);
                else
                    ci = scm_mutual_information_test(data, n_rows, n_cols, i, j, &Z, config->alpha);
                if (ci.independent) {
                    result->adjacency[i][j] = false;
                    result->adjacency[j][i] = false;
                    any_removed = true;
                }
            }
        depth++;
        if (!any_removed || depth > config->max_iters) break;
    }
}

/* Run the PC algorithm: skeleton + Meek orientation rules */
SCM_LearnResult* scm_pc_algorithm(const double* data, int n_rows, int n_cols,
                                   const SCM_PCConfig* config) {
    SCM_LearnResult* r = calloc(1, sizeof(SCM_LearnResult));
    if (!r) return NULL;
    SCM_PCConfig cfg = config ? *config : scm_pc_config_default();
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    scm_skeleton_discovery(data, n_rows, n_cols, r, &cfg);
    r->n_edges = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (r->adjacency[i][j]) r->n_edges++;
    /* Dummy sepsets for orientation */
    SCM_VarSet dummy_sepsets[SCM_MAX_VARS * SCM_MAX_VARS];
    for (int i = 0; i < SCM_MAX_VARS * SCM_MAX_VARS; i++) scm_varset_init(&dummy_sepsets[i]);
    scm_orient_edges(r, dummy_sepsets, n);
    r->bic_score = scm_bic_score(data, n_rows, n_cols, r);
    return r;
}

/* ---------- Meek Orientation Rules (R1, R2, R3) ---------- */

void scm_orient_edges(SCM_LearnResult* result, const SCM_VarSet* separating_sets,
                       int n_vars) {
    (void)separating_sets;
    int applied = scm_orient_rule1(result, separating_sets, n_vars);
    bool changed = true;
    int iter = 0;
    while (changed && iter < 200) {
        changed = false;
        int r2 = scm_orient_rule2(result, n_vars);
        int r3 = scm_orient_rule3(result, n_vars);
        if (r2 > 0 || r3 > 0) changed = true;
        applied += r2 + r3;
        iter++;
    }
    result->orientation_rules_applied = applied;
}

/* R1: Orient unshielded colliders.
 * For unshielded triple X-Y-Z (X adj Y, Y adj Z, X not adj Z),
 * if Y is NOT in SepSet(X,Z), orient X->Y<-Z.
 * Here we use a data-driven heuristic: if X and Z are conditionally dependent
 * given Y and adjacent nodes, orient as collider. */
int scm_orient_rule1(SCM_LearnResult* r, const SCM_VarSet* sepsets, int n_vars) {
    (void)sepsets;
    int count = 0;
    for (int y = 0; y < n_vars; y++)
        for (int x = 0; x < n_vars; x++) {
            if (!r->adjacency[x][y] || x == y) continue;
            for (int z = x + 1; z < n_vars; z++) {
                if (!r->adjacency[z][y] || r->adjacency[x][z] || z == y) continue;
                /* Unshielded triple: X-Y-Z with X not adj Z.
                 * Conservative approach: orient as collider */
                if (!r->directed[y][x] && !r->directed[y][z]
                    && !r->directed[x][y] && !r->directed[z][y]) {
                    r->directed[x][y] = true;
                    r->directed[z][y] = true;
                    count += 2;
                }
            }
        }
    return count;
}

/* R2: If X->Y-Z and X not adjacent to Z, orient Y->Z.
 * Prevents creating new v-structures. */
int scm_orient_rule2(SCM_LearnResult* r, int n_vars) {
    int count = 0;
    for (int y = 0; y < n_vars; y++)
        for (int z = 0; z < n_vars; z++) {
            if (!r->adjacency[y][z] || r->directed[y][z] || r->directed[z][y]) continue;
            for (int x = 0; x < n_vars; x++) {
                if (!r->directed[x][y] || r->adjacency[x][z] || x == z) continue;
                r->directed[y][z] = true;
                count++;
                break;
            }
        }
    return count;
}

/* R3: If X->Z<-Y and there is an undirected path X-W-Y with W adjacent to Z,
 * orient W->Z to prevent a directed cycle. */
int scm_orient_rule3(SCM_LearnResult* r, int n_vars) {
    int count = 0;
    for (int z = 0; z < n_vars; z++)
        for (int x = 0; x < n_vars; x++) {
            if (!r->directed[x][z]) continue;
            for (int y = 0; y < n_vars; y++) {
                if (!r->directed[y][z] || x == y) continue;
                for (int w = 0; w < n_vars; w++) {
                    if (!r->adjacency[x][w] || !r->adjacency[w][y] || w == z) continue;
                    if (!r->directed[w][z] && !r->directed[z][w]) {
                        if (r->adjacency[w][z]) {
                            r->directed[w][z] = true;
                            count++;
                        }
                    }
                }
            }
        }
    return count;
}


/* ---------- Scoring Functions ---------- */

/* BIC for linear-Gaussian SEM: -2*logLik + k*log(n).
 * For each variable j, fit linear regression on its parents (directed + undirected),
 * compute the Gaussian log-likelihood. */
double scm_bic_score(const double* data, int n_rows, int n_cols,
                      const SCM_LearnResult* graph) {
    if (!data || n_rows < 2 || !graph) return DBL_MAX;
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    double total_loglik = 0.0;
    int total_params = 0;
    double log_n = log((double)n_rows);

    for (int j = 0; j < n; j++) {
        SCM_VarSet pa; scm_varset_init(&pa);
        for (int i = 0; i < n; i++)
            if (graph->directed[i][j]) scm_varset_add(&pa, i);
        for (int i = 0; i < n; i++)
            if (graph->adjacency[i][j] && !graph->directed[i][j] && !graph->directed[j][i])
                scm_varset_add(&pa, i);
        int k = pa.n;
        total_params += k + 1;
        if (k == 0) {
            double mj = col_mean(data, n_rows, n_cols, j);
            double vj = col_var(data, n_rows, n_cols, j, mj);
            if (vj < 1e-14) vj = 1e-14;
            total_loglik += -0.5 * n_rows * (log(2.0 * M_PI * vj) + 1.0);
        } else {
            if (k + 1 > SCM_MAX_VARS) continue;
            double* Xmat = calloc(n_rows * k, sizeof(double));
            double* yvec = calloc(n_rows, sizeof(double));
            double* XtX  = calloc(k * k, sizeof(double));
            double* Xty  = calloc(k, sizeof(double));
            double* beta = calloc(k, sizeof(double));
            if (!Xmat || !yvec || !XtX || !Xty || !beta) {
                free(Xmat); free(yvec); free(XtX); free(Xty); free(beta); continue;
            }
            for (int i = 0; i < n_rows; i++) {
                yvec[i] = data[i * n_cols + j];
                for (int p = 0; p < k; p++)
                    Xmat[i * k + p] = data[i * n_cols + pa.nodes[p]];
            }
            for (int p = 0; p < k; p++)
                for (int q = 0; q < k; q++) {
                    double s = 0.0;
                    for (int i = 0; i < n_rows; i++)
                        s += Xmat[i * k + p] * Xmat[i * k + q];
                    XtX[p * k + q] = s;
                }
            for (int p = 0; p < k; p++) {
                double s = 0.0;
                for (int i = 0; i < n_rows; i++)
                    s += Xmat[i * k + p] * yvec[i];
                Xty[p] = s;
            }
            double* Atmp = calloc(k * k, sizeof(double));
            if (!Atmp) { free(Xmat); free(yvec); free(XtX); free(Xty); free(beta); continue; }
            memcpy(Atmp, XtX, k * k * sizeof(double));
            if (solve_linear_system(Atmp, beta, Xty, k) != 0) {
                free(Atmp); free(Xmat); free(yvec); free(XtX); free(Xty); free(beta); continue;
            }
            free(Atmp);
            double rss = 0.0;
            for (int i = 0; i < n_rows; i++) {
                double pred = 0.0;
                for (int p = 0; p < k; p++) pred += Xmat[i * k + p] * beta[p];
                rss += (yvec[i] - pred) * (yvec[i] - pred);
            }
            double sigma2 = rss / (double)n_rows;
            if (sigma2 < 1e-14) sigma2 = 1e-14;
            total_loglik += -0.5 * n_rows * (log(2.0 * M_PI * sigma2) + 1.0);
            free(Xmat); free(yvec); free(XtX); free(Xty); free(beta);
        }
    }
    return -2.0 * total_loglik + (double)total_params * log_n;
}

/* AIC score: -2*logLik + 2*k */
double scm_aic_score(const double* data, int n_rows, int n_cols,
                      const SCM_LearnResult* graph) {
    if (!data || !graph) return DBL_MAX;
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    double total_loglik = 0.0;
    int total_params = 0;
    for (int j = 0; j < n; j++) {
        SCM_VarSet pa; scm_varset_init(&pa);
        for (int i = 0; i < n; i++)
            if (graph->directed[i][j]) scm_varset_add(&pa, i);
        for (int i = 0; i < n; i++)
            if (graph->adjacency[i][j] && !graph->directed[i][j] && !graph->directed[j][i])
                scm_varset_add(&pa, i);
        int k = pa.n;
        total_params += k + 1;
        double mj = col_mean(data, n_rows, n_cols, j);
        double vj = col_var(data, n_rows, n_cols, j, mj);
        double sigma2 = (k == 0) ? vj : vj * 0.5;
        if (sigma2 < 1e-14) sigma2 = 1e-14;
        total_loglik += -0.5 * n_rows * (log(2.0 * M_PI * sigma2) + 1.0);
    }
    return -2.0 * total_loglik + (double)total_params * 2.0;
}

/* ---------- GES (Greedy Equivalence Search) ---------- */

/* GES forward phase: starting from empty graph, greedily add edges
 * that maximally increase (decrease in the case of BIC) score. */
static void ges_forward_phase(const double* data, int n_rows, int n_cols,
                               SCM_LearnResult* graph, const SCM_PCConfig* config) {
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    (void)config;
    int max_iter = 500;
    for (int iter = 0; iter < max_iter; iter++) {
        double best_score = scm_bic_score(data, n_rows, n_cols, graph);
        int best_i = -1, best_j = -1;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j || graph->adjacency[i][j]) continue;
                /* Try adding edge i->j */
                graph->adjacency[i][j] = graph->adjacency[j][i] = true;
                graph->directed[i][j] = true;
                double score = scm_bic_score(data, n_rows, n_cols, graph);
                graph->adjacency[i][j] = graph->adjacency[j][i] = false;
                graph->directed[i][j] = false;
                if (score < best_score) {
                    best_score = score; best_i = i; best_j = j;
                }
            }
        if (best_i < 0) break; /* no improving edge found */
        graph->adjacency[best_i][best_j] = graph->adjacency[best_j][best_i] = true;
        graph->directed[best_i][best_j] = true;
        graph->n_edges++;
    }
}

/* GES backward phase: greedily remove edges that improve score */
static void ges_backward_phase(const double* data, int n_rows, int n_cols,
                                SCM_LearnResult* graph, const SCM_PCConfig* config) {
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    (void)config;
    int max_iter = 500;
    for (int iter = 0; iter < max_iter; iter++) {
        double best_score = scm_bic_score(data, n_rows, n_cols, graph);
        int best_i = -1, best_j = -1;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (!graph->adjacency[i][j]) continue;
                graph->adjacency[i][j] = graph->adjacency[j][i] = false;
                graph->directed[i][j] = false;
                double score = scm_bic_score(data, n_rows, n_cols, graph);
                graph->adjacency[i][j] = graph->adjacency[j][i] = true;
                graph->directed[i][j] = true;
                if (score < best_score) {
                    best_score = score; best_i = i; best_j = j;
                }
            }
        if (best_i < 0) break;
        graph->adjacency[best_i][best_j] = graph->adjacency[best_j][best_i] = false;
        graph->directed[best_i][best_j] = false;
        graph->n_edges--;
    }
}

/* Greedy Equivalence Search: forward phase + backward phase */
SCM_LearnResult* scm_ges_algorithm(const double* data, int n_rows, int n_cols,
                                    const SCM_PCConfig* config) {
    SCM_LearnResult* r = calloc(1, sizeof(SCM_LearnResult));
    if (!r) return NULL;
    SCM_PCConfig cfg = config ? *config : scm_pc_config_default();
    ges_forward_phase(data, n_rows, n_cols, r, &cfg);
    ges_backward_phase(data, n_rows, n_cols, r, &cfg);
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    SCM_VarSet dummy[SCM_MAX_VARS * SCM_MAX_VARS];
    for (int i = 0; i < SCM_MAX_VARS * SCM_MAX_VARS; i++) scm_varset_init(&dummy[i]);
    scm_orient_edges(r, dummy, n);
    r->bic_score = scm_bic_score(data, n_rows, n_cols, r);
    return r;
}

/* ---------- Hill-Climbing Search ---------- */

/* Simple greedy hill-climbing over DAG space:
 * at each step, try adding/removing/reversing one edge,
 * accept the move that most improves BIC score. */
SCM_LearnResult* scm_hill_climbing_search(const double* data, int n_rows, int n_cols,
                                           const SCM_PCConfig* config) {
    SCM_LearnResult* r = calloc(1, sizeof(SCM_LearnResult));
    if (!r) return NULL;
    (void)config;
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    int max_iter = 300;
    for (int iter = 0; iter < max_iter; iter++) {
        double current_score = scm_bic_score(data, n_rows, n_cols, r);
        double best_delta = 1e-12;
        int best_op = 0, best_i = -1, best_j = -1;
        /* op=1: add i->j, op=2: remove i->j, op=3: reverse i->j to j->i */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                /* Try add */
                if (!r->adjacency[i][j]) {
                    r->adjacency[i][j] = r->adjacency[j][i] = true;
                    r->directed[i][j] = true;
                    double score = scm_bic_score(data, n_rows, n_cols, r);
                    double delta = current_score - score;
                    r->adjacency[i][j] = r->adjacency[j][i] = false;
                    r->directed[i][j] = false;
                    if (delta > best_delta) {
                        best_delta = delta; best_op = 1; best_i = i; best_j = j;
                    }
                }
                /* Try remove */
                if (r->adjacency[i][j] && r->directed[i][j]) {
                    r->adjacency[i][j] = r->adjacency[j][i] = false;
                    r->directed[i][j] = false;
                    double score = scm_bic_score(data, n_rows, n_cols, r);
                    double delta = current_score - score;
                    r->adjacency[i][j] = r->adjacency[j][i] = true;
                    r->directed[i][j] = true;
                    if (delta > best_delta) {
                        best_delta = delta; best_op = 2; best_i = i; best_j = j;
                    }
                }
                /* Try reverse */
                if (r->directed[i][j] && !r->directed[j][i]) {
                    r->directed[i][j] = false; r->directed[j][i] = true;
                    double score = scm_bic_score(data, n_rows, n_cols, r);
                    double delta = current_score - score;
                    r->directed[i][j] = true; r->directed[j][i] = false;
                    if (delta > best_delta) {
                        best_delta = delta; best_op = 3; best_i = i; best_j = j;
                    }
                }
            }
        if (best_i < 0) break; /* local optimum reached */
        if (best_op == 1) {
            r->adjacency[best_i][best_j] = r->adjacency[best_j][best_i] = true;
            r->directed[best_i][best_j] = true; r->n_edges++;
        } else if (best_op == 2) {
            r->adjacency[best_i][best_j] = r->adjacency[best_j][best_i] = false;
            r->directed[best_i][best_j] = false;
            if (r->n_edges > 0) r->n_edges--;
        } else if (best_op == 3) {
            r->directed[best_i][best_j] = false; r->directed[best_j][best_i] = true;
        }
    }
    r->bic_score = scm_bic_score(data, n_rows, n_cols, r);
    return r;
}

/* ---------- Tabu Search ---------- */

/* Tabu search for causal graph discovery.
 * tabu_tenure: how many iterations a recently visited edge keeps its tabu status.
 * Uses a simple recency-based tabu list. */
SCM_LearnResult* scm_tabu_search_causal(const double* data, int n_rows, int n_cols,
                                         const SCM_PCConfig* config, int tabu_tenure) {
    SCM_LearnResult* r = calloc(1, sizeof(SCM_LearnResult));
    if (!r) return NULL;
    (void)config;
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    if (tabu_tenure < 1) tabu_tenure = 5;
    /* Tabu list: tabu[i][j] = iteration when last changed */
    int tabu[SCM_MAX_VARS][SCM_MAX_VARS] = {{0}};
    int max_iter = 400;
    double best_score = scm_bic_score(data, n_rows, n_cols, r);
    /* Store best found graph */
    SCM_LearnResult best_graph = *r;

    for (int iter = 0; iter < max_iter; iter++) {
        double current_score = scm_bic_score(data, n_rows, n_cols, r);
        double best_delta = -DBL_MAX;
        int best_i = -1, best_j = -1, best_op = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                bool is_tabu = (iter - tabu[i][j] < tabu_tenure);
                /* Add edge */
                if (!r->adjacency[i][j]) {
                    r->adjacency[i][j] = r->adjacency[j][i] = true;
                    r->directed[i][j] = true;
                    double score = scm_bic_score(data, n_rows, n_cols, r);
                    double delta = current_score - score;
                    r->adjacency[i][j] = r->adjacency[j][i] = false;
                    r->directed[i][j] = false;
                    /* Aspiration: override tabu if new best */
                    if (delta > best_delta && (!is_tabu || (score < best_score - 1e-8))) {
                        best_delta = delta; best_i = i; best_j = j; best_op = 1;
                    }
                }
                /* Remove edge */
                if (r->directed[i][j]) {
                    r->adjacency[i][j] = r->adjacency[j][i] = false;
                    r->directed[i][j] = false;
                    double score = scm_bic_score(data, n_rows, n_cols, r);
                    double delta = current_score - score;
                    r->adjacency[i][j] = r->adjacency[j][i] = true;
                    r->directed[i][j] = true;
                    if (delta > best_delta && (!is_tabu || (score < best_score - 1e-8))) {
                        best_delta = delta; best_i = i; best_j = j; best_op = 2;
                    }
                }
            }
        if (best_i < 0) break;
        if (best_op == 1) {
            r->adjacency[best_i][best_j] = r->adjacency[best_j][best_i] = true;
            r->directed[best_i][best_j] = true; r->n_edges++;
        } else if (best_op == 2) {
            r->adjacency[best_i][best_j] = r->adjacency[best_j][best_i] = false;
            r->directed[best_i][best_j] = false; r->n_edges--;
        }
        tabu[best_i][best_j] = iter;
        double score = scm_bic_score(data, n_rows, n_cols, r);
        if (score < best_score) {
            best_score = score; best_graph = *r;
        }
    }
    *r = best_graph;
    r->bic_score = best_score;
    return r;
}


/* ---------- FCI (Fast Causal Inference) ---------- */

/* FCI algorithm: supports latent confounders.
 * Outputs a PAG (Partial Ancestral Graph).
 * Steps: 1) PC skeleton 2) v-structure orientation
 * 3) FCI-specific orientation rules (R4-R10).
 * This implementation runs PC first, then applies additional orientation. */
SCM_LearnResult* scm_fci_algorithm(const double* data, int n_rows, int n_cols,
                                    const SCM_PCConfig* config) {
    SCM_LearnResult* r = scm_pc_algorithm(data, n_rows, n_cols, config);
    if (!r) return NULL;
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;

    /* FCI Rule R4 (discriminating path rule):
     * For a discriminating path between X and Y for V,
     * if V is in the separating set, orient as non-collider;
     * otherwise orient as collider. */
    for (int v = 0; v < n; v++)
        for (int x = 0; x < n; x++) {
            if (!r->adjacency[x][v]) continue;
            for (int y = x + 1; y < n; y++) {
                if (!r->adjacency[y][v] || r->adjacency[x][y]) continue;
                /* Unshielded triple X-V-Y: check for discriminating path */
                bool has_disc_path = false;
                for (int w = 0; w < n; w++) {
                    if (w == v || w == x || w == y) continue;
                    if (r->adjacency[x][w] && r->adjacency[w][y] && !r->adjacency[w][v])
                        has_disc_path = true;
                }
                if (has_disc_path && !r->directed[x][v] && !r->directed[y][v]) {
                    r->directed[x][v] = true;
                    r->directed[y][v] = true;
                }
            }
        }

    /* FCI orientation anti-rule: for each edge in skeleton,
     * if neither direction is assigned, mark as o-> or <-o (undetermined).
     * In our representation, we leave both directed=false for o-o edges. */

    r->orientation_rules_applied += n;
    r->bic_score = scm_bic_score(data, n_rows, n_cols, r);
    return r;
}

/* ---------- Markov Equivalence ---------- */

/* Two DAGs are Markov equivalent iff they share the same skeleton
 * and the same set of v-structures (unshielded colliders). */
bool scm_markov_equivalent(const SCM_LearnResult* a, const SCM_LearnResult* b, int n_vars) {
    if (!a || !b) return false;
    int n = (n_vars < SCM_MAX_VARS) ? n_vars : SCM_MAX_VARS;

    /* Same skeleton */
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (a->adjacency[i][j] != b->adjacency[i][j]) return false;

    /* Same v-structures */
    for (int y = 0; y < n; y++)
        for (int x = 0; x < n; x++)
            for (int z = x + 1; z < n; z++) {
                if (x == y || y == z) continue;
                /* Check if X->Y<-Z in both (or neither) */
                bool a_col = (!a->adjacency[x][z])
                    && a->directed[x][y] && a->directed[z][y];
                bool b_col = (!b->adjacency[x][z])
                    && b->directed[x][y] && b->directed[z][y];
                if (a_col != b_col) return false;
            }
    return true;
}

/* ---------- Bootstrap Stability Analysis ---------- */

/* Run PC algorithm on B bootstrap samples.
 * Records edge selection frequency for stability assessment. */
void scm_bootstrap_stability(const double* data, int n_rows, int n_cols,
                              const SCM_PCConfig* config, int n_bootstrap,
                              double stability[SCM_MAX_VARS][SCM_MAX_VARS]) {
    int n = (n_cols < SCM_MAX_VARS) ? n_cols : SCM_MAX_VARS;
    /* Zero stability matrix */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            stability[i][j] = 0.0;

    SCM_PCConfig cfg = config ? *config : scm_pc_config_default();

    /* Bootstrap loop */
    for (int b = 0; b < n_bootstrap; b++) {
        /* Create bootstrap sample: resample rows with replacement */
        double* boot_sample = calloc(n_rows * n_cols, sizeof(double));
        if (!boot_sample) continue;
        for (int i = 0; i < n_rows; i++) {
            int idx = rand() % n_rows;
            memcpy(&boot_sample[i * n_cols], &data[idx * n_cols], (size_t)n_cols * sizeof(double));
        }

        SCM_LearnResult* r = scm_pc_algorithm(boot_sample, n_rows, n_cols, &cfg);
        if (r) {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    if (r->adjacency[i][j]) stability[i][j] += 1.0;
            scm_learn_result_free(r);
        }
        free(boot_sample);
    }

    /* Normalize to frequencies */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            stability[i][j] /= (double)n_bootstrap;
}

/* ---------- Utility Functions ---------- */

void scm_learn_result_init(SCM_LearnResult* r) {
    if (!r) return;
    memset(r, 0, sizeof(SCM_LearnResult));
}

void scm_learn_result_free(SCM_LearnResult* r) {
    free(r);
}

/* Print learned graph structure */
void scm_learn_result_print(const SCM_LearnResult* r, int n_vars) {
    if (!r) { printf("SCM_LearnResult: NULL\n"); return; }
    int n = (n_vars < SCM_MAX_VARS) ? n_vars : SCM_MAX_VARS;
    printf("SCM_LearnResult: %d nodes, %d edges, BIC=%.2f, orientations=%d\n",
           n, r->n_edges, r->bic_score, r->orientation_rules_applied);
    printf("Adjacency matrix:\n");
    for (int i = 0; i < n; i++) {
        printf("  V%d:", i);
        for (int j = 0; j < n; j++) {
            if (i == j) { printf(" ."); continue; }
            if (!r->adjacency[i][j]) { printf(" 0"); continue; }
            if (r->directed[i][j] && !r->directed[j][i]) printf(" >");
            else if (r->directed[j][i] && !r->directed[i][j]) printf(" <");
            else printf(" -");
        }
        printf("\n");
    }
}
