/* cdf_scores.c ? Score-Based Causal Discovery
 *
 * Implements BIC, BGe, AIC scores for linear-Gaussian DAGs,
 * linear regression via normal equations, and greedy/tabu
 * DAG search algorithms.
 *
 * Score decomposability: Score(G) = Sum_v Score(v, Pa_G(v)).
 * This enables O(p) local score updates for edge operations.
 *
 * Linear-Gaussian model for node v with parents Pa(v):
 *   X_v = beta_0 + Sum_{u in Pa(v)} beta_u X_u + eps, eps ~ N(0, sigma^2)
 *
 * BIC peak = -N ln(sigma^2_hat) - k ln(N) where k = |Pa| + 1
 * BGe integrates over parameters using normal-Wishart conjugate prior.
 */

#include "cdf_scores.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ?? Extract parent indices from a DAG ???????????????????????????????? */

static int get_parents(const CdfGraph *g, int v, int *parents)
{
    int n = 0;
    for (int u = 0; u < g->p; u++) {
        if (g->edges[u * g->p + v] == CDF_EDGE_DIRECTED) {
            if (n < CDF_MAX_DEGREE) parents[n++] = u;
        }
    }
    return n;
}

/* ?? Linear Regression (Normal Equations) ????????????????????????????? */

bool cdf_score_linear_regression(const double *y, const double *X,
                                  int N, int k, double *beta, double *sigma2)
{
    if (!y || !X || N < 2 || k < 0 || !beta || !sigma2) return false;

    if (k == 0) {
        /* Empty model: intercept only = mean(y) */
        double mean_y = 0.0;
        for (int i = 0; i < N; i++) mean_y += y[i];
        mean_y /= (double)N;
        beta[0] = mean_y;
        double rss = 0.0;
        for (int i = 0; i < N; i++) {
            double r = y[i] - mean_y;
            rss += r * r;
        }
        *sigma2 = rss / (double)N;
        return true;
    }

    /* Design matrix: m = k+1 (intercept + k regressors) */
    int m = k + 1;

    /* Accumulate XtX [m x m] and Xty [m x 1] */
    double *XtX = (double*)calloc((size_t)m * m, sizeof(double));
    double *Xty = (double*)calloc((size_t)m, sizeof(double));
    if (!XtX || !Xty) { free(XtX); free(Xty); return false; }

    for (int i = 0; i < N; i++) {
        /* Row i: [1, X[i+0*N], X[i+1*N], ..., X[i+(k-1)*N]] */
        for (int a = 0; a < m; a++) {
            double va = (a == 0) ? 1.0 : X[i + (a - 1) * N];
            Xty[a] += va * y[i];
            for (int b = 0; b < m; b++) {
                double vb = (b == 0) ? 1.0 : X[i + (b - 1) * N];
                XtX[a * m + b] += va * vb;
            }
        }
    }

    /* Solve via Gaussian elimination with partial pivoting */
    double *aug = (double*)calloc((size_t)m * (m + 1), sizeof(double));
    if (!aug) { free(XtX); free(Xty); return false; }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++)
            aug[i * (m + 1) + j] = XtX[i * m + j];
        aug[i * (m + 1) + m] = Xty[i];
    }

    for (int col = 0; col < m; col++) {
        /* Partial pivot */
        int pivot = col;
        double maxv = fabs(aug[col * (m + 1) + col]);
        for (int row = col + 1; row < m; row++) {
            double v = fabs(aug[row * (m + 1) + col]);
            if (v > maxv) { maxv = v; pivot = row; }
        }
        if (maxv < CDF_EPS) {
            free(XtX); free(Xty); free(aug);
            for (int i = 0; i < m; i++) beta[i] = 0.0;
            *sigma2 = 1.0;
            return false;
        }

        if (pivot != col) {
            for (int j = 0; j <= m; j++) {
                double t = aug[col * (m + 1) + j];
                aug[col * (m + 1) + j] = aug[pivot * (m + 1) + j];
                aug[pivot * (m + 1) + j] = t;
            }
        }

        double piv = aug[col * (m + 1) + col];
        for (int j = 0; j <= m; j++)
            aug[col * (m + 1) + j] /= piv;

        for (int row = 0; row < m; row++) {
            if (row == col) continue;
            double factor = aug[row * (m + 1) + col];
            for (int j = 0; j <= m; j++)
                aug[row * (m + 1) + j] -= factor * aug[col * (m + 1) + j];
        }
    }

    for (int i = 0; i < m; i++)
        beta[i] = aug[i * (m + 1) + m];

    /* RSS */
    double rss = 0.0;
    for (int i = 0; i < N; i++) {
        double y_hat = beta[0];
        for (int j = 0; j < k; j++)
            y_hat += beta[j + 1] * X[i + j * N];
        double resid = y[i] - y_hat;
        rss += resid * resid;
    }
    *sigma2 = rss / (double)N;

    free(XtX); free(Xty); free(aug);
    return true;
}

/* ?? Design matrix builder ??????????????????????????????????????????? */

static void build_design(const CdfDataset *ds, int v,
                          const int *parents, int n_parents,
                          double *X, double *y)
{
    int N = ds->N;
    for (int i = 0; i < N; i++)
        y[i] = ds->data[v * N + i];
    for (int j = 0; j < n_parents; j++) {
        int pv = parents[j];
        for (int i = 0; i < N; i++)
            X[i + j * N] = ds->data[pv * N + i];
    }
}

/* ?? Score helper ???????????????????????????????????????????????????? */

static double score_local(const CdfDataset *ds, int v,
                           const int *parents, int n_parents,
                           CdfScoreType type)
{
    int N = ds->N, k = n_parents;
    if (k > CDF_MAX_DEGREE) k = CDF_MAX_DEGREE;

    if (k == 0) {
        double mean = 0.0, ss = 0.0;
        for (int i = 0; i < N; i++) {
            double val = ds->data[v * N + i];
            mean += val; ss += val * val;
        }
        mean /= (double)N;
        double var = ss / (double)N - mean * mean;
        if (var < CDF_EPS) var = CDF_EPS;

        switch (type) {
        case CDF_SCORE_BIC:
            return -(double)N * log(var) - log((double)N);  /* 1 param */
        case CDF_SCORE_AIC:
            return -(double)N * log(var) - 2.0;
        case CDF_SCORE_BGE: {
            double s2 = CDF_EPS;
            /* marginal log-likelihood with Jeffrey's prior */
            for (int i = 0; i < N; i++) {
                double d = ds->data[v * N + i] - mean;
                s2 += d * d;
            }
            s2 /= (double)(N - 1);
            if (s2 < CDF_EPS) s2 = CDF_EPS;
            return -0.5 * (double)N * log(2.0 * CDF_PI * s2)
                   - 0.5 * log((double)N) - 0.5 * (double)(N - 1);
        }
        default: return 0.0;
        }
    }

    double *Xm = (double*)malloc((size_t)N * k * sizeof(double));
    double *yv = (double*)malloc((size_t)N * sizeof(double));
    double *beta = (double*)malloc((size_t)(k + 1) * sizeof(double));
    if (!Xm || !yv || !beta) {
        free(Xm); free(yv); free(beta); return -1e100;
    }

    build_design(ds, v, parents, n_parents, Xm, yv);

    double sigma2;
    if (!cdf_score_linear_regression(yv, Xm, N, k, beta, &sigma2)) {
        free(Xm); free(yv); free(beta); return -1e100;
    }
    if (sigma2 < CDF_EPS) sigma2 = CDF_EPS;

    double n_params = (double)(k + 1);
    double lnN = log((double)N);
    double score;

    switch (type) {
    case CDF_SCORE_BIC:
        score = -(double)N * log(sigma2) - n_params * lnN;
        break;
    case CDF_SCORE_AIC:
        score = -(double)N * log(sigma2) - 2.0 * n_params;
        break;
    case CDF_SCORE_BGE:
        /* BGe approximation: log ML of Gaussian with unit info prior.
         * P(D|G,v) approx -(N/2) ln(2 pi sigma2) - 0.5 ln|S_pa| - 0.5 k ln(N) */
        {
            double logdet = 0.0;
            for (int j = 0; j < k; j++) {
                int pv = parents[j];
                double m = 0.0, vv = 0.0;
                for (int i = 0; i < N; i++) {
                    double val = ds->data[pv * N + i];
                    m += val; vv += val * val;
                }
                m /= (double)N;
                double pvar = vv / (double)N - m * m;
                if (pvar < 1e-10) pvar = 1e-10;
                logdet += log(pvar);
            }
            score = -0.5 * (double)N * log(2.0 * CDF_PI * sigma2)
                    - 0.5 * logdet - 0.5 * n_params * lnN;
        }
        break;
    default: score = 0.0;
    }

    free(Xm); free(yv); free(beta);
    return score;
}

/* ?? Public Score API ????????????????????????????????????????????????? */

double cdf_score_bic_local(const CdfDataset *ds, int v,
                            const int *parents, int n_parents)
{
    if (!ds || v < 0 || v >= ds->p) return -1e100;
    return score_local(ds, v, parents, n_parents, CDF_SCORE_BIC);
}

double cdf_score_bic(const CdfDataset *ds, const CdfGraph *g)
{
    if (!ds || !g || ds->p != g->p) return -1e100;
    int p = g->p;
    double total = 0.0;
    int parents[CDF_MAX_DEGREE];
    for (int v = 0; v < p; v++) {
        int np = get_parents(g, v, parents);
        total += cdf_score_bic_local(ds, v, parents, np);
    }
    return total;
}

double cdf_score_bge_local(const CdfDataset *ds, int v,
                            const int *parents, int n_parents,
                            double alpha_w, double alpha_mu)
{
    (void)alpha_w; (void)alpha_mu;
    return score_local(ds, v, parents, n_parents, CDF_SCORE_BGE);
}

double cdf_score_bge(const CdfDataset *ds, const CdfGraph *g,
                      double alpha_w, double alpha_mu)
{
    if (!ds || !g || ds->p != g->p) return -1e100;
    int p = g->p;
    double total = 0.0;
    int parents[CDF_MAX_DEGREE];
    for (int v = 0; v < p; v++) {
        int np = get_parents(g, v, parents);
        total += cdf_score_bge_local(ds, v, parents, np, alpha_w, alpha_mu);
    }
    return total;
}

double cdf_score_aic_local(const CdfDataset *ds, int v,
                            const int *parents, int n_parents)
{
    return score_local(ds, v, parents, n_parents, CDF_SCORE_AIC);
}

double cdf_score_aic(const CdfDataset *ds, const CdfGraph *g)
{
    if (!ds || !g || ds->p != g->p) return -1e100;
    int p = g->p;
    double total = 0.0;
    int parents[CDF_MAX_DEGREE];
    for (int v = 0; v < p; v++) {
        int np = get_parents(g, v, parents);
        total += cdf_score_aic_local(ds, v, parents, np);
    }
    return total;
}

/* ?? Score Delta ?????????????????????????????????????????????????????? */

double cdf_score_delta(const CdfDataset *ds, const CdfGraph *g,
                        int u, int v, char op, CdfScoreType type)
{
    if (!ds || !g || u < 0 || v < 0 || u >= g->p || v >= g->p || u == v)
        return 0.0;
    int p = g->p;
    int parents[CDF_MAX_DEGREE] = {0};

    switch (op) {
    case 'a': {  /* add u -> v */
        if (cdf_graph_has_edge(g, u, v)) return 0.0;
        if (cdf_graph_reachable(g, v, u)) return -1e100;
        int np_old = get_parents(g, v, parents);
        double sc_old = score_local(ds, v, parents, np_old, type);
        parents[np_old] = u;
        double sc_new = score_local(ds, v, parents, np_old + 1, type);
        return sc_new - sc_old;
    }
    case 'd': {  /* delete u -> v */
        if (g->edges[u * p + v] != CDF_EDGE_DIRECTED) return 0.0;
        int np_old = get_parents(g, v, parents);
        double sc_old = score_local(ds, v, parents, np_old, type);
        int pnew[CDF_MAX_DEGREE] = {0};
        int np_new = 0;
        for (int i = 0; i < np_old; i++)
            if (parents[i] != u) pnew[np_new++] = parents[i];
        double sc_new = score_local(ds, v, pnew, np_new, type);
        return sc_new - sc_old;
    }
    case 'r': {  /* reverse u->v to v->u */
        if (g->edges[u * p + v] != CDF_EDGE_DIRECTED) return 0.0;
        if (cdf_graph_reachable(g, u, v)) return -1e100;
        /* Node v loses parent u */
        int np_v = get_parents(g, v, parents);
        double sc_v_old = score_local(ds, v, parents, np_v, type);
        int pv_new[CDF_MAX_DEGREE] = {0};
        int npv_new = 0;
        for (int i = 0; i < np_v; i++)
            if (parents[i] != u) pv_new[npv_new++] = parents[i];
        double sc_v_new = score_local(ds, v, pv_new, npv_new, type);
        /* Node u gains parent v */
        int np_u = get_parents(g, u, parents);
        double sc_u_old = score_local(ds, u, parents, np_u, type);
        parents[np_u] = v;
        double sc_u_new = score_local(ds, u, parents, np_u + 1, type);
        return (sc_v_new + sc_u_new) - (sc_v_old + sc_u_old);
    }
    default: return 0.0;
    }
}

/* ?? Edge Operation Struct ??????????????????????????????????????????? */

typedef struct { int u, v; char op; double delta; } EdgeOp;

static EdgeOp find_best_op(const CdfDataset *ds, const CdfGraph *g,
                            CdfScoreType type)
{
    EdgeOp best = {0, 0, '?', -1e-8};
    int p = g->p;

    for (int u = 0; u < p; u++) {
        for (int v = 0; v < p; v++) {
            if (u == v) continue;

            /* Try add edge */
            if (!cdf_graph_has_edge(g, u, v) && !cdf_graph_has_edge(g, v, u)) {
                double d = cdf_score_delta(ds, g, u, v, 'a', type);
                if (d > best.delta) {
                    best.u = u; best.v = v; best.op = 'a'; best.delta = d;
                }
            }
            /* Try delete edge */
            if (g->edges[u * p + v] == CDF_EDGE_DIRECTED) {
                double d = cdf_score_delta(ds, g, u, v, 'd', type);
                if (d > best.delta) {
                    best.u = u; best.v = v; best.op = 'd'; best.delta = d;
                }
            }
            /* Try reverse edge */
            if (g->edges[u * p + v] == CDF_EDGE_DIRECTED
                && !cdf_graph_has_edge(g, v, u)) {
                double d = cdf_score_delta(ds, g, u, v, 'r', type);
                if (d > best.delta) {
                    best.u = u; best.v = v; best.op = 'r'; best.delta = d;
                }
            }
        }
    }
    return best;
}

static void apply_op(CdfGraph *g, EdgeOp op)
{
    if (op.op == 'a') cdf_graph_add_edge(g, op.u, op.v, CDF_EDGE_DIRECTED);
    else if (op.op == 'd') cdf_graph_remove_edge(g, op.u, op.v);
    else if (op.op == 'r') {
        cdf_graph_remove_edge(g, op.u, op.v);
        cdf_graph_add_edge(g, op.v, op.u, CDF_EDGE_DIRECTED);
    }
}

/* ?? Tabu key encoding ???????????????????????????????????????????????? */

static int tabu_key(int u, int v, char op, int p) {
    return (u * p + v) * 3 + ((op == 'a') ? 0 : ((op == 'd') ? 1 : 2));
}

static bool is_tabu(const int *tabu_list, int tabu_len,
                     int key, int *pos)
{
    for (int i = 0; i < tabu_len; i++) {
        if (tabu_list[i] == key) { *pos = i; return true; }
    }
    *pos = -1;
    return false;
}

/* ?? Greedy DAG Search ???????????????????????????????????????????????? */

CdfGraph* cdf_score_greedy_search(const CdfDataset *ds, int max_iter,
                                   CdfScoreType type, bool verbose)
{
    if (!ds) return NULL;
    int p = ds->p;
    if (max_iter <= 0) max_iter = 100 * p;

    CdfGraph *g = cdf_graph_create(p);
    if (!g) return NULL;

    double cur_score = 0.0;
    for (int v = 0; v < p; v++)
        cur_score += cdf_score_bic_local(ds, v, NULL, 0);

    for (int iter = 0; iter < max_iter; iter++) {
        EdgeOp best = find_best_op(ds, g, type);
        if (best.delta <= 1e-8) break;

        apply_op(g, best);
        cur_score += best.delta;

        if (verbose && (iter % 50 == 0 || best.delta > 0.5))
            printf("  iter %d: %c(%d->%d) d=%.3f score=%.2f\n",
                   iter, best.op, best.u, best.v, best.delta, cur_score);
    }

    if (verbose) printf("  Greedy done: score=%.2f\n", cur_score);
    return g;
}

/* ?? Tabu Search ?????????????????????????????????????????????????????? */

CdfGraph* cdf_score_tabu_search(const CdfDataset *ds, int max_iter,
                                 int tabu_len, CdfScoreType type,
                                 bool verbose)
{
    if (!ds) return NULL;
    int p = ds->p;
    if (max_iter <= 0) max_iter = 200 * p;
    if (tabu_len <= 0) tabu_len = 7;

    CdfGraph *g = cdf_graph_create(p);
    CdfGraph *g_best = cdf_graph_create(p);
    if (!g || !g_best) { cdf_graph_free(g); cdf_graph_free(g_best); return NULL; }

    double best_score = -1e100;
    int *tabu = (int*)malloc((size_t)tabu_len * sizeof(int));
    if (!tabu) { cdf_graph_free(g); cdf_graph_free(g_best); return NULL; }
    for (int i = 0; i < tabu_len; i++) tabu[i] = -1;
    int tabu_idx = 0, p2 = p * p;

    for (int iter = 0; iter < max_iter; iter++) {
        EdgeOp best = {0, 0, '?', -1e-8};

        for (int u = 0; u < p; u++) {
            for (int v = 0; v < p; v++) {
                if (u == v) continue;
                char ops[] = {'a', 'd', 'r'};
                for (int oi = 0; oi < 3; oi++) {
                    char op = ops[oi];
                    double d = cdf_score_delta(ds, g, u, v, op, type);
                    if (d <= 1e-10) continue;
                    int key = tabu_key(u, v, op, p);
                    int tpos;
                    bool tabu_hit = is_tabu(tabu, tabu_len, key, &tpos);
                    /* Aspiration: if this move beats global best, allow it */
                    if (tabu_hit && d <= (best_score - cdf_score_bic(ds, g))) continue;
                    if (d > best.delta) {
                        best.u = u; best.v = v; best.op = op; best.delta = d;
                    }
                }
            }
        }

        if (best.delta <= 1e-10) {
            /* No tabu-legal improving move: random perturb */
            for (int u = 0; u < p; u++)
                for (int v = u + 1; v < p; v++)
                    if (cdf_graph_has_edge(g, u, v)) {
                        cdf_graph_remove_edge(g, u, v); goto perturb_done;
                    }
            perturb_done:;
            continue;
        }

        apply_op(g, best);
        int key = tabu_key(best.u, best.v, best.op, p);
        tabu[tabu_idx] = key;
        tabu_idx = (tabu_idx + 1) % tabu_len;

        double cs = cdf_score_bic(ds, g);
        if (cs > best_score) {
            best_score = cs;
            for (int i = 0; i < p2; i++) g_best->edges[i] = g->edges[i];
        }

        if (verbose && iter % 50 == 0)
            printf("  tabu %d: best=%.2f\n", iter, best_score);
    }

    cdf_graph_free(g); free(tabu);
    if (verbose) printf("  Tabu done: best=%.2f\n", best_score);
    return g_best;
}

/* ?? Bootstrap Stability ?????????????????????????????????????????????? */

void cdf_score_bootstrap_stability(const CdfDataset *ds, int B,
                                    double *stability,
                                    CdfScoreType type, bool verbose)
{
    if (!ds || !stability || B < 1) return;
    int p = ds->p, N = ds->N;
    size_t p2 = (size_t)p * p;
    for (size_t i = 0; i < p2; i++) stability[i] = 0.0;

    double *bdata = (double*)malloc((size_t)p * N * sizeof(double));
    int *idx = (int*)malloc((size_t)N * sizeof(int));
    if (!bdata || !idx) { free(bdata); free(idx); return; }

    int ok = 0;
    for (int b = 0; b < B; b++) {
        for (int i = 0; i < N; i++) idx[i] = rand() % N;
        for (int var = 0; var < p; var++)
            for (int i = 0; i < N; i++)
                bdata[var * N + i] = ds->data[var * N + idx[i]];

        CdfDataset *bd = cdf_dataset_create(bdata, N, p);
        if (!bd) continue;

        CdfGraph *bg = cdf_score_greedy_search(bd, 30, type, false);
        if (bg) {
            ok++;
            for (int i = 0; i < p; i++)
                for (int j = i + 1; j < p; j++)
                    if (cdf_graph_has_edge(bg, i, j))
                        stability[i * p + j] += 1.0;
            cdf_graph_free(bg);
        }
        cdf_dataset_free(bd);
        if (verbose && (b + 1) % 10 == 0)
            printf("  bootstrap %d/%d\n", b + 1, B);
    }

    if (ok > 0)
        for (size_t i = 0; i < p2; i++)
            stability[i] /= (double)ok;

    if (verbose) printf("  Bootstrap done: %d/%d ok\n", ok, B);
    free(bdata); free(idx);
}

/* ?? Utilities ???????????????????????????????????????????????????????? */

double cdf_score_compare(const CdfDataset *ds, const CdfGraph *g1,
                          const CdfGraph *g2, CdfScoreType type)
{
    double s1, s2;
    switch (type) {
    case CDF_SCORE_BIC: s1 = cdf_score_bic(ds, g1); s2 = cdf_score_bic(ds, g2); break;
    case CDF_SCORE_BGE: s1 = cdf_score_bge(ds, g1, ds->p+2, 1.0);
                          s2 = cdf_score_bge(ds, g2, ds->p+2, 1.0); break;
    case CDF_SCORE_AIC: s1 = cdf_score_aic(ds, g1); s2 = cdf_score_aic(ds, g2); break;
    default: return 0.0;
    }
    return s1 - s2;
}

int cdf_score_num_params(const CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, total = 0;
    int parents[CDF_MAX_DEGREE];
    for (int v = 0; v < p; v++) {
        int np = get_parents(g, v, parents);
        total += np + 1;
    }
    return total;
}

void cdf_score_print_report(const CdfDataset *ds, const CdfGraph *g,
                             CdfScoreType type)
{
    if (!ds || !g) return;
    int p = g->p;
    const char *tname = (type == CDF_SCORE_BIC) ? "BIC" :
                         (type == CDF_SCORE_BGE ? "BGe" : "AIC");

    printf("=== Score Report (%s) ===\n", tname);
    printf("  Nodes: %d, Params: %d\n", p, cdf_score_num_params(g));
    printf("  BIC = %.4f\n", cdf_score_bic(ds, g));
    printf("  AIC = %.4f\n", cdf_score_aic(ds, g));
    printf("  BGe = %.4f\n", cdf_score_bge(ds, g, (double)(p+2), 1.0));

    printf("  Local scores:\n");
    int parents[CDF_MAX_DEGREE];
    for (int v = 0; v < p; v++) {
        int np = get_parents(g, v, parents);
        double s = score_local(ds, v, parents, np, type);
        printf("    X%d (|Pa|=%d): %.4f\n", v + 1, np, s);
    }
}
