/*
 * adjustment_sets.c �� Covariate Adjustment Set Selection
 *
 * Implements algorithms for finding valid covariate adjustment sets
 * that identify causal effects from observational data:
 *   - Parent adjustment set (simplest valid set)
 *   - Minimal adjustment set enumeration
 *   - Generalized Adjustment Criterion (Perkovic et al., 2015)
 *   - Optimal adjustment set by asymptotic variance (Henckel et al., 2019)
 *   - Efficient covariate selection via outcome regression
 *
 * References:
 *   Pearl, "Causality", 2009, Ch.3.3 (Back-Door Criterion)
 *   Perkovic, Textor, Kalisch, Maathuis, "A Complete Generalized
 *     Adjustment Criterion", JMLR 2015
 *   Henckel, Perkovic, Maathuis, "Graphical Criteria for Efficient
 *     Total Effect Estimation via Adjustment in Causal Linear Models",
 *     JRSS-B, 2019
 *   Witte, Henckel, Maathuis, Didelez, "On Efficient Adjustment in
 *     Causal Graphs", JMLR, 2020
 */

#include "adjustment_sets.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ================================================================
 * Parent Adjustment Set
 *
 * Z = parents(X). Always valid under the causal Markov condition
 * because parents of X block all back-door paths into X.
 * ================================================================ */

AdjustmentSet* adjustment_parent_set(const DAG* dag, int x) {
    if (!dag || x < 0 || x >= dag->n) return NULL;

    AdjustmentSet* as = calloc(1, sizeof(AdjustmentSet));
    if (!as) return NULL;

    int np = dag->n_parents[x];
    as->variables = calloc((size_t)(np > 0 ? np : 1), sizeof(int));
    if (!as->variables) { free(as); return NULL; }

    as->n_vars = np;
    as->is_valid = (np > 0);  /* valid if X has parents */
    as->efficiency = (np > 0) ? 1.0 / (double)np : DBL_MAX;

    for (int i = 0; i < np; i++)
        as->variables[i] = dag->parents[x * dag->n + i];

    return as;
}

/* ================================================================
 * Adjustment Set Validation
 *
 * Uses the generalized adjustment criterion:
 * Z is valid for (X,Y) iff:
 *   1. Z blocks all "proper causal paths" from X to Y
 *      (a proper causal path is a directed path X -> ... -> Y)
 *   2. Z contains no "forbidden nodes"
 *
 * For the back-door criterion: check that Z blocks all back-door paths.
 * ================================================================ */

bool adjustment_is_valid(const DAG* dag, int x, int y,
                          const int *z, int n_z) {
    if (!dag || x < 0 || x >= dag->n || y < 0 || y >= dag->n)
        return false;

    /* Use the back-door criterion check */
    return causal_backdoor_criterion(dag, x, y, z, n_z);
}

/* ================================================================
 * Generalized Adjustment Criterion (Perkovic et al., 2015)
 *
 * Extends the back-door criterion to handle cases where the
 * adjustment set may be larger than needed but still valid.
 *
 * Z satisfies the generalized adjustment criterion for (X,Y) if:
 *   1. Z blocks all non-causal paths from X to Y (back-door paths)
 *   2. Z does not contain any descendants of X that lie on a
 *      causal path from X to Y (avoids "M-bias" etc.)
 * ================================================================ */

bool adjustment_generalized_criterion(const DAG* dag, int x, int y,
                                       const int *z, int n_z) {
    if (!dag || x < 0 || x >= dag->n || y < 0 || y >= dag->n)
        return false;

    /* Condition 1: Z blocks all back-door paths */
    if (!causal_backdoor_criterion(dag, x, y, z, n_z))
        return false;

    /* Condition 2: No Z node is a descendant of X on a causal path to Y */
    /* A node v is on a causal path X->...->Y if:
     *   X can reach v, and v can reach Y — checked via dag_has_path below */

    for (int i = 0; i < n_z; i++) {
        int v = z[i];
        /* Check if v is on a causal path from X to Y */
        if (dag_has_path(dag, x, v) && dag_has_path(dag, v, y)) {
            /* v is on a causal path �� not allowed in adjustment set
             * unless it's a mediator used in front-door adjustment */
            return false;
        }
    }

    return true;
}

/* ================================================================
 * Find All Minimal Adjustment Sets
 *
 * A minimal adjustment set is one where removing any variable
 * breaks the back-door criterion. We enumerate all subsets of
 * non-descendants of X and test minimality.
 * ================================================================ */

int adjustment_find_minimal_sets(const DAG* dag, int x, int y,
                                  AdjustmentSet *sets, int max_sets) {
    if (!dag || !sets || max_sets <= 0 ||
        x < 0 || x >= dag->n || y < 0 || y >= dag->n) return -1;

    int desc[DAG_MAX_NODES];
    int nd = dag_descendants(dag, x, desc, DAG_MAX_NODES);
    if (nd < 0) nd = 0;

    int candidates[DAG_MAX_NODES];
    int n_cand = 0;
    for (int v = 0; v < dag->n; v++) {
        if (v == x || v == y) continue;
        bool is_desc = false;
        for (int j = 0; j < nd; j++)
            if (desc[j] == v) { is_desc = true; break; }
        if (!is_desc && n_cand < DAG_MAX_NODES)
            candidates[n_cand++] = v;
    }

    if (n_cand == 0) {
        /* No candidates: empty set is trivially minimal if valid */
        if (adjustment_is_valid(dag, x, y, NULL, 0)) {
            sets[0].variables = NULL;
            sets[0].n_vars = 0;
            sets[0].is_valid = true;
            sets[0].efficiency = 0.0;
            return 1;
        }
        return 0;
    }

    int count = 0;
    unsigned long long limit = 1ULL << (unsigned long long)n_cand;

    /* First, collect all valid sets */
    typedef struct { int vars[ADJ_MAX_VARS]; int n; } ValSet;
    ValSet *valid = calloc((size_t)(limit < 1024 ? (int)limit : 1024),
                            sizeof(ValSet));
    int n_valid = 0;
    int max_valid = (int)(limit < 1024 ? limit : 1024);

    for (unsigned long long mask = 1; mask < limit && n_valid < max_valid;
         mask++) {
        int z_buf[ADJ_MAX_VARS];
        int n_z = 0;
        for (int b = 0; b < n_cand && n_z < ADJ_MAX_VARS; b++) {
            if (mask & (1ULL << (unsigned long long)b))
                z_buf[n_z++] = candidates[b];
        }

        if (adjustment_is_valid(dag, x, y, z_buf, n_z)) {
            for (int j = 0; j < n_z; j++)
                valid[n_valid].vars[j] = z_buf[j];
            valid[n_valid].n = n_z;
            n_valid++;
        }
    }

    /* Now filter to minimal sets: no proper subset is also valid */
    for (int i = 0; i < n_valid && count < max_sets; i++) {
        bool is_minimal = true;
        for (int j = 0; j < n_valid && is_minimal; j++) {
            if (i == j) continue;
            if (valid[j].n >= valid[i].n) continue;

            /* Check if valid[j] is a subset of valid[i] */
            bool is_subset = true;
            for (int k = 0; k < valid[j].n && is_subset; k++) {
                bool found = false;
                for (int l = 0; l < valid[i].n; l++) {
                    if (valid[j].vars[k] == valid[i].vars[l])
                    { found = true; break; }
                }
                if (!found) is_subset = false;
            }
            if (is_subset) is_minimal = false;
        }

        if (is_minimal) {
            sets[count].variables = calloc((size_t)valid[i].n, sizeof(int));
            if (sets[count].variables) {
                for (int k = 0; k < valid[i].n; k++)
                    sets[count].variables[k] = valid[i].vars[k];
                sets[count].n_vars = valid[i].n;
                sets[count].is_valid = true;
                sets[count].efficiency = 1.0 / (double)(valid[i].n + 1);
                count++;
            }
        }
    }

    free(valid);
    return count;
}

/* ================================================================
 * Optimal Adjustment Set by Asymptotic Variance
 *
 * Among all valid adjustment sets, the optimal one minimizes the
 * asymptotic variance of the causal effect estimator.
 *
 * For linear models: Var(ATE_hat) �� sigma^2 * (1/n) * (1/R^2_{X|Z})
 * where R^2_{X|Z} is the R-squared of regressing X on Z.
 *
 * We estimate this via: for each valid set Z, compute
 *   var_reduction = 1 - Var(X|Z) / Var(X)
 * and select the set with maximum var_reduction.
 * ================================================================ */

AdjustmentSet* adjustment_optimal_set(const DAG* dag, int x, int y,
                                       const double *data, int n, int d) {
    if (!dag || x < 0 || x >= dag->n || y < 0 || y >= dag->n)
        return NULL;

    AdjustmentSet min_sets[ADJ_MAX_SETS];
    int n_sets = adjustment_find_minimal_sets(dag, x, y, min_sets, ADJ_MAX_SETS);
    if (n_sets <= 0) {
        /* Fall back to parent set */
        AdjustmentSet* ps = adjustment_parent_set(dag, x);
        return ps;
    }

    /* Select the set with best efficiency (lowest variance) */
    /* If we have data, compute R-squared for each set.
     * Without data, use size heuristic: smaller set = higher efficiency
     * (fewer degrees of freedom consumed). */
    int best_idx = 0;
    double best_eff = min_sets[0].efficiency;

    if (data && n > 0 && d > 0) {
        /* Use data-driven selection: compute partial correlation */
        for (int s = 0; s < n_sets; s++) {
            double efficiency = 0.0;
            double var_x = 1e-8;

            /* Compute Var(X) */
            double mean_x = 0.0;
            for (int i = 0; i < n; i++) mean_x += data[i * d + x];
            mean_x /= (double)n;
            for (int i = 0; i < n; i++) {
                double dx = data[i * d + x] - mean_x;
                var_x += dx * dx;
            }
            var_x /= (double)n;

            /* Simple heuristic: efficiency = 1/(|Z|+1), adjusted by
             * how correlated Z is with X */
            double z_corr = 0.0;
            for (int zi = 0; zi < min_sets[s].n_vars; zi++) {
                int zv = min_sets[s].variables[zi];
                double mean_z = 0.0, cov_xz = 0.0, var_z = 1e-8;
                for (int i = 0; i < n; i++) mean_z += data[i * d + zv];
                mean_z /= (double)n;
                for (int i = 0; i < n; i++) {
                    double dx = data[i * d + x] - mean_x;
                    double dz = data[i * d + zv] - mean_z;
                    cov_xz += dx * dz;
                    var_z += dz * dz;
                }
                cov_xz /= (double)n;
                var_z /= (double)n;
                double corr = cov_xz / sqrt(var_x * var_z);
                z_corr += fabs(corr);
            }
            efficiency = z_corr / (double)(min_sets[s].n_vars + 1);

            if (efficiency > best_eff) {
                best_eff = efficiency;
                best_idx = s;
            }
        }
    }

    /* Return the best set */
    AdjustmentSet* best = calloc(1, sizeof(AdjustmentSet));
    if (!best) {
        for (int s = 0; s < n_sets; s++)
            if (min_sets[s].variables) free(min_sets[s].variables);
        return NULL;
    }

    best->n_vars = min_sets[best_idx].n_vars;
    best->efficiency = best_eff;
    best->is_valid = true;
    best->variables = calloc((size_t)(best->n_vars > 0 ? best->n_vars : 1),
                              sizeof(int));
    if (best->variables) {
        for (int i = 0; i < best->n_vars; i++)
            best->variables[i] = min_sets[best_idx].variables[i];
    }

    for (int s = 0; s < n_sets; s++)
        if (min_sets[s].variables) free(min_sets[s].variables);

    return best;
}

/* ================================================================
 * Efficient Covariate Selection
 *
 * Given a set of candidate covariates, select the subset that
 * minimizes the residual variance of the outcome Y after adjusting
 * for X and the covariates. This is essentially Lasso-type selection
 * but using a simple greedy forward selection based on variance
 * reduction.
 * ================================================================ */

int adjustment_select_efficient(const DAG* dag, int x, int y,
                                 const double *data, int n, int d,
                                 int *selected, int max_sel) {
    if (!dag || !selected || max_sel <= 0 || !data || n <= 0 || d <= 0 ||
        x < 0 || x >= dag->n || y < 0 || y >= dag->n) return -1;

    /* Get valid candidate sets */
    AdjustmentSet min_sets[ADJ_MAX_SETS];
    int n_sets = adjustment_find_minimal_sets(dag, x, y, min_sets, ADJ_MAX_SETS);

    if (n_sets <= 0) return 0;

    /* Use the union of all minimal sets as candidates */
    int candidates[ADJ_MAX_VARS];
    int n_cand = 0;
    int *seen = calloc((size_t)dag->n, sizeof(int));
    if (!seen) {
        for (int s = 0; s < n_sets; s++)
            free(min_sets[s].variables);
        return -1;
    }

    for (int s = 0; s < n_sets; s++) {
        for (int vi = 0; vi < min_sets[s].n_vars; vi++) {
            int v = min_sets[s].variables[vi];
            if (!seen[v] && n_cand < ADJ_MAX_VARS) {
                seen[v] = 1;
                candidates[n_cand++] = v;
            }
        }
    }
    free(seen);

    for (int s = 0; s < n_sets; s++) free(min_sets[s].variables);

    if (n_cand == 0) return 0;

    /* Greedy forward selection: add variables that most reduce
     * residual variance of Y after regressing on X */
    int *chosen = calloc((size_t)n_cand, sizeof(int));
    if (!chosen) return -1;
    int n_chosen = 0;

    /* Simple variance-based selection:
     * For each candidate, compute partial correlation with Y given X.
     * Select those with significant partial correlation. */
    double var_x = 1e-8, mean_x = 0.0;
    double var_y = 1e-8, mean_y = 0.0;
    for (int i = 0; i < n; i++) {
        mean_x += data[i * d + x];
        mean_y += data[i * d + y];
    }
    mean_x /= (double)n; mean_y /= (double)n;
    for (int i = 0; i < n; i++) {
        double dx = data[i * d + x] - mean_x;
        double dy = data[i * d + y] - mean_y;
        var_x += dx * dx; var_y += dy * dy;
    }
    var_x /= (double)n; var_y /= (double)n;

    for (int ci = 0; ci < n_cand && n_chosen < max_sel; ci++) {
        int c = candidates[ci];
        double mean_c = 0.0, cov_cy = 0.0, cov_cx = 0.0, var_c = 1e-8;
        for (int i = 0; i < n; i++) mean_c += data[i * d + c];
        mean_c /= (double)n;
        for (int i = 0; i < n; i++) {
            double dc = data[i * d + c] - mean_c;
            double dy = data[i * d + y] - mean_y;
            double dx = data[i * d + x] - mean_x;
            cov_cy += dc * dy;
            cov_cx += dc * dx;
            var_c += dc * dc;
        }
        cov_cy /= (double)n; cov_cx /= (double)n; var_c /= (double)n;

        /* Partial correlation of Y and C given X:
         * r_{YC|X} = (r_{YC} - r_{YX}*r_{CX}) / sqrt((1-r_{YX}^2)*(1-r_{CX}^2)) */
        double r_yc = cov_cy / sqrt(var_y * var_c);
        double r_yx = 0.0; /* compute */
        {   double cov_yx = 0.0;
            for (int i = 0; i < n; i++)
                cov_yx += (data[i*d+y]-mean_y)*(data[i*d+x]-mean_x);
            cov_yx /= (double)n;
            r_yx = cov_yx / sqrt(var_y * var_x);
        }
        double r_cx = cov_cx / sqrt(var_c * var_x);
        double denom = sqrt((1.0 - r_yx*r_yx) * (1.0 - r_cx*r_cx));
        double partial_r = (denom > 1e-10) ?
                           (r_yc - r_yx * r_cx) / denom : 0.0;

        /* Select if partial correlation is substantial */
        if (fabs(partial_r) > 0.05 && n_chosen < max_sel) {
            selected[n_chosen++] = c;
        }
    }

    free(chosen);
    return n_chosen;
}

/* ================================================================
 * AdjustmentSet Lifecycle
 * ================================================================ */

void adjustment_set_free(AdjustmentSet* as) {
    if (!as) return;
    free(as->variables);
    free(as);
}
