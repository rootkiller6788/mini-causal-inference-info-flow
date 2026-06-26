/* cdf_ida.c -- IDA: Intervention-calculus when the DAG is Absent
 *
 * Implements IDA (Maathuis, Kalisch, Buhlmann, 2009; 2010) for
 * estimating causal effects from observational data when the
 * true DAG is unknown.
 *
 * Given a learned CPDAG from PC, IDA:
 *  1. Enumerates all valid parent sets of Y consistent with CPDAG
 *  2. For each, estimates causal effect: Y ~ X + Pa(Y)
 *  3. Returns multiset of possible effects as bounds
 */

#include "cdf_ida.h"
#include "cdf_citest.h"
#include "cdf_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

static bool is_par(const CdfGraph *cp, int y, int u) {
    return cp->edges[u * cp->p + y] == CDF_EDGE_DIRECTED;
}
static bool is_und(const CdfGraph *cp, int y, int u) {
    return cp->edges[y * cp->p + u] == CDF_EDGE_UNDIRECTED;
}

/* Enumerate valid parent sets from a CPDAG for node y. */
void cdf_ida_enumerate_parent_sets(const CdfGraph *cp, int y,
                                    int **parent_sets, int *n_sets,
                                    int max_sets)
{
    if (!cp || !parent_sets || !n_sets || max_sets < 1) return;
    int p = cp->p;
    *n_sets = 0;

    int mandatory[CDF_MAX_DEGREE], n_mand = 0;
    int optional[CDF_MAX_DEGREE], n_opt = 0;

    for (int u = 0; u < p; u++) {
        if (u == y) continue;
        if (is_par(cp, y, u)) {
            if (n_mand < CDF_MAX_DEGREE) mandatory[n_mand++] = u;
        } else if (is_und(cp, y, u)) {
            if (n_opt < CDF_MAX_DEGREE) optional[n_opt++] = u;
        }
    }

    int n_combos = 1 << n_opt;
    for (int mask = 0; mask < n_combos && *n_sets < max_sets; mask++) {
        int np = n_mand;
        for (int i = 0; i < n_mand; i++)
            parent_sets[*n_sets][i] = mandatory[i];
        for (int b = 0; b < n_opt; b++)
            if (mask & (1 << b))
                parent_sets[*n_sets][np++] = optional[b];
        if (np <= CDF_MAX_DEGREE) {
            parent_sets[*n_sets][np] = -1;
            (*n_sets)++;
        }
    }
}

/* Causal effect via linear regression Y ~ beta0 + beta_x*X + sum_j beta_j*Z_j */
double cdf_ida_effect_estimate(const CdfDataset *ds, int x, int y,
                                const int *Z, int nZ)
{
    if (!ds || x < 0 || y < 0 || x >= ds->p || y >= ds->p) return 0.0;
    if (nZ < 0) return 0.0;
    int N = ds->N;
    int k = nZ + 1;
    if (k < 1) return 0.0;

    double *Xm = (double*)malloc((size_t)N * k * sizeof(double));
    double *Yv = (double*)malloc((size_t)N * sizeof(double));
    if (!Xm || !Yv) { free(Xm); free(Yv); return 0.0; }

    for (int i = 0; i < N; i++) Xm[i] = ds->data[x * N + i];
    for (int j = 0; j < nZ; j++) {
        int zj = Z[j];
        if (zj < 0 || zj >= ds->p) continue;
        for (int i = 0; i < N; i++)
            Xm[i + (j + 1) * N] = ds->data[zj * N + i];
    }
    for (int i = 0; i < N; i++) Yv[i] = ds->data[y * N + i];

    int m = k + 1;
    double *XtX = (double*)calloc((size_t)m * m, sizeof(double));
    double *Xty = (double*)calloc((size_t)m, sizeof(double));
    if (!XtX || !Xty) { free(Xm); free(Yv); free(XtX); free(Xty); return 0.0; }

    for (int i = 0; i < N; i++) {
        double row[64];
        row[0] = 1.0;
        row[1] = Xm[i];
        for (int j = 1; j < k; j++) row[j + 1] = Xm[i + j * N];
        for (int a = 0; a < m; a++)
            for (int b = 0; b < m; b++)
                XtX[a * m + b] += row[a] * row[b];
        for (int a = 0; a < m; a++)
            Xty[a] += row[a] * Yv[i];
    }

    double *aug = (double*)calloc((size_t)m * (m + 1), sizeof(double));
    if (!aug) { free(Xm); free(Yv); free(XtX); free(Xty); return 0.0; }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++)
            aug[i * (m + 1) + j] = XtX[i * m + j];
        aug[i * (m + 1) + m] = Xty[i];
    }

    for (int col = 0; col < m; col++) {
        int pivot = col;
        double maxv = fabs(aug[col * (m + 1) + col]);
        for (int r = col + 1; r < m; r++) {
            double v = fabs(aug[r * (m + 1) + col]);
            if (v > maxv) { maxv = v; pivot = r; }
        }
        if (maxv < 1e-12) {
            free(Xm); free(Yv); free(XtX); free(Xty); free(aug); return 0.0;
        }
        if (pivot != col) {
            for (int j = 0; j <= m; j++) {
                double t = aug[col * (m + 1) + j];
                aug[col * (m + 1) + j] = aug[pivot * (m + 1) + j];
                aug[pivot * (m + 1) + j] = t;
            }
        }
        double piv = aug[col * (m + 1) + col];
        for (int j = 0; j <= m; j++) aug[col * (m + 1) + j] /= piv;
        for (int r = 0; r < m; r++) {
            if (r == col) continue;
            double f = aug[r * (m + 1) + col];
            for (int j = 0; j <= m; j++)
                aug[r * (m + 1) + j] -= f * aug[col * (m + 1) + j];
        }
    }

    double beta_x = aug[1 * (m + 1) + m];
    free(Xm); free(Yv); free(XtX); free(Xty); free(aug);
    return beta_x;
}

/* Back-door criterion check */
bool cdf_ida_is_backdoor_admissible(const CdfGraph *g, int x, int y,
                                     const int *Z, int nZ)
{
    if (!g) return false;
    for (int i = 0; i < nZ; i++)
        if (cdf_graph_reachable(g, x, Z[i])) return false;
    int p = g->p;
    for (int u = 0; u < p; u++) {
        if (u == x || u == y) continue;
        if (g->edges[u * p + x] != CDF_EDGE_DIRECTED) continue;
        if (!cdf_graph_d_separated(g, u, y, Z, nZ)) return false;
    }
    return true;
}

/* Front-door criterion check */
bool cdf_ida_is_frontdoor_admissible(const CdfGraph *g, int x, int y,
                                      const int *Z, int nZ)
{
    if (!g || nZ < 1) return false;
    for (int i = 0; i < nZ; i++) {
        if (!cdf_graph_reachable(g, x, Z[i])) return false;
        if (!cdf_graph_reachable(g, Z[i], y)) return false;
    }
    if (!cdf_graph_d_separated(g, x, Z[0], NULL, 0)) return false;
    int cx[1] = {x};
    for (int i = 0; i < nZ; i++)
        if (!cdf_graph_d_separated(g, Z[i], y, cx, 1)) return false;
    return true;
}

/* Full IDA estimation */
CdfIDAResult* cdf_ida_estimate(const CdfGraph *cp, const CdfDataset *ds,
                                int x, int y)
{
    if (!cp || !ds || x < 0 || y < 0 || x >= cp->p || y >= cp->p) return NULL;
    if (x == y) return NULL;

    int p = cp->p;
    CdfIDAResult *res = (CdfIDAResult*)calloc(1, sizeof(CdfIDAResult));
    if (!res) return NULL;

    int max_sets = (p <= 10) ? (1 << p) : 1024;
    int **psets = (int**)malloc((size_t)max_sets * sizeof(int*));
    if (!psets) { free(res); return NULL; }
    for (int i = 0; i < max_sets; i++) {
        psets[i] = (int*)malloc((size_t)(p + 1) * sizeof(int));
        if (!psets[i]) {
            for (int j = 0; j < i; j++) free(psets[j]);
            free(psets); free(res); return NULL;
        }
    }

    int n_sets = 0;
    cdf_ida_enumerate_parent_sets(cp, y, psets, &n_sets, max_sets);
    if (n_sets == 0) {
        for (int i = 0; i < max_sets; i++) free(psets[i]);
        free(psets); free(res); return NULL;
    }

    res->effects = (CdfIDAEffect*)calloc((size_t)n_sets, sizeof(CdfIDAEffect));
    if (!res->effects) {
        for (int i = 0; i < max_sets; i++) free(psets[i]);
        free(psets); free(res); return NULL;
    }
    res->n_effects = n_sets;
    res->min_effect = 1e100; res->max_effect = -1e100;
    res->all_same_sign = true;
    int sign_ref = 0;

    for (int s = 0; s < n_sets; s++) {
        int pa[CDF_MAX_DEGREE], np = 0;
        for (int i = 0; i < p && psets[s][i] >= 0; i++)
            if (psets[s][i] != x && np < CDF_MAX_DEGREE)
                pa[np++] = psets[s][i];
        double eff = cdf_ida_effect_estimate(ds, x, y, pa, np);
        res->effects[s].x = x; res->effects[s].y = y;
        res->effects[s].effect_coef = eff;
        res->effects[s].n_parents = np;
        res->effects[s].parents_used = (int*)malloc((size_t)(np+1)*sizeof(int));
        if (res->effects[s].parents_used) {
            for (int i = 0; i < np; i++)
                res->effects[s].parents_used[i] = pa[i];
            res->effects[s].parents_used[np] = -1;
        }
        if (eff < res->min_effect) res->min_effect = eff;
        if (eff > res->max_effect) res->max_effect = eff;
        int cs = (eff > 1e-9) ? 1 : ((eff < -1e-9) ? -1 : 0);
        if (s == 0) sign_ref = cs;
        else if (cs != sign_ref) res->all_same_sign = false;
    }

    for (int i = 0; i < n_sets; i++)
        for (int j = i + 1; j < n_sets; j++)
            if (res->effects[i].effect_coef > res->effects[j].effect_coef) {
                CdfIDAEffect t = res->effects[i];
                res->effects[i] = res->effects[j]; res->effects[j] = t;
            }
    if (n_sets % 2 == 1)
        res->median_effect = res->effects[n_sets/2].effect_coef;
    else
        res->median_effect = 0.5*(res->effects[n_sets/2-1].effect_coef
                                 + res->effects[n_sets/2].effect_coef);

    for (int i = 0; i < max_sets; i++) free(psets[i]);
    free(psets);
    return res;
}

void cdf_ida_result_free(CdfIDAResult *res)
{
    if (!res) return;
    if (res->effects) {
        for (int i = 0; i < res->n_effects; i++)
            free(res->effects[i].parents_used);
        free(res->effects);
    }
    free(res);
}

void cdf_ida_print(const CdfIDAResult *res)
{
    if (!res) return;
    printf("=== IDA: X%d -> X%d ===\n",
           res->effects[0].x + 1, res->effects[0].y + 1);
    printf("  Sets: %d  Range: [%.4f, %.4f]  Median: %.4f\n",
           res->n_effects, res->min_effect, res->max_effect,
           res->median_effect);
    printf("  Sign consistent: %s\n", res->all_same_sign ? "YES" : "NO");
}
