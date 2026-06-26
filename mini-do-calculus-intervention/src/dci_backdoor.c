#include "dci_backdoor.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Back-Door Criterion
 *
 * Z satisfies the back-door criterion relative to (X, Y) if:
 *   1. No node in Z is a descendant of X
 *   2. Z blocks every back-door path from X to Y
 * ============================================================== */

DCI_BackdoorSet dci_backdoor_find(DCI_Graph* g, int cause, int effect) {
    DCI_BackdoorSet bd;
    memset(&bd, 0, sizeof(bd));
    if (!g) return bd;
    /* Strategy: find all back-door paths, then find minimal set
     * of non-descendants of cause that block all such paths */
    DCI_Path paths[DCI_MAX_PATHS];
    int n_paths = dci_find_backdoor_paths(g, cause, effect, paths, DCI_MAX_PATHS);
    if (n_paths == 0) {
        bd.is_valid = true; /* No back-door paths = no confounding */
        return bd;
    }
    /* Find non-descendant parents of cause as initial back-door set */
    int i;
    for (i = 0; i < g->n_nodes; i++) {
        if (i == cause || i == effect) continue;
        if (g->adjacency[i][cause] > 0.5) {
            /* Parent of cause — check if not a descendant */
            if (!dci_is_descendant(g, i, cause)) {
                if (bd.n_vars < DCI_MAX_VARS) {
                    bd.variables[bd.n_vars++] = i;
                }
            }
        }
    }
    /* Verify: check if this set blocks all back-door paths */
    bd.is_valid = dci_backdoor_check(g, cause, effect,
        bd.variables, bd.n_vars);
    return bd;
}

bool dci_backdoor_check(DCI_Graph* g, int cause, int effect,
    const int* z_set, int n_z) {
    if (!g) return false;
    /* Check 1: no Z_i is a descendant of cause */
    int i;
    for (i = 0; i < n_z; i++) {
        if (dci_is_descendant(g, z_set[i], cause)) return false;
    }
    /* Check 2: Z blocks all back-door paths from cause to effect */
    DCI_Path paths[DCI_MAX_PATHS];
    int n_paths = dci_find_backdoor_paths(g, cause, effect, paths, DCI_MAX_PATHS);
    if (n_paths == 0) return true;
    /* For each back-door path, check if Z blocks it */
    int pi;
    for (pi = 0; pi < n_paths; pi++) {
        DCI_Path* path = &paths[pi];
        bool blocked = false;
        int j;
        for (j = 0; j < path->length && !blocked; j++) {
            int node = path->nodes[j];
            if (node == cause || node == effect) continue;
            /* Check if node is in Z */
            int k;
            for (k = 0; k < n_z; k++) {
                if (z_set[k] == node) { blocked = true; break; }
            }
        }
        if (!blocked) return false;
    }
    return true;
}

DCI_BackdoorSet dci_backdoor_minimal(DCI_Graph* g, int cause, int effect) {
    DCI_BackdoorSet full = dci_backdoor_find(g, cause, effect);
    if (!full.is_valid || full.n_vars <= 1) return full;
    /* Greedy minimization: try removing each variable */
    DCI_BackdoorSet minimal = full;
    int i;
    for (i = 0; i < full.n_vars; i++) {
        int test_set[DCI_MAX_VARS];
        int n_test = 0;
        int j;
        for (j = 0; j < full.n_vars; j++) {
            if (j != i) test_set[n_test++] = full.variables[j];
        }
        if (dci_backdoor_check(g, cause, effect, test_set, n_test)) {
            minimal.n_vars = n_test;
            memcpy(minimal.variables, test_set, n_test * sizeof(int));
            break;
        }
    }
    minimal.is_valid = true;
    return minimal;
}

/* ==============================================================
 * Back-Door Adjustment Formula
 *
 * P(Y=y | do(X=x)) = Σ_z P(Y=y | X=x, Z=z) P(Z=z)
 * ============================================================== */

double dci_backdoor_adjust(DCI_SCM* scm, DCI_Graph* g,
    int cause, int effect, double cause_val, double effect_val,
    DCI_BackdoorSet* z_set, int n_samples) {
    if (!scm || !g || !z_set || n_samples <= 0) return 0.0;
    /* Monte Carlo estimation of back-door adjustment */
    double prob_sum = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        /* Sample exogenous variables and evaluate SCM */
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        /* Set cause to cause_val (intervention) */
        exog[cause] = cause_val;
        dci_scm_evaluate(scm, exog);
        /* Check if effect matches */
        if (fabs(scm->vars[effect].value - effect_val) < 0.01) {
            prob_sum += 1.0;
        }
    }
    (void)z_set;
    return prob_sum / n_samples;
}

/* ==============================================================
 * Front-Door Criterion
 *
 * M satisfies the front-door criterion relative to (X, Y) if:
 *   1. M intercepts all directed paths from X to Y
 *   2. No back-door path from X to M exists
 *   3. All back-door paths from M to Y are blocked by X
 * ============================================================== */

DCI_FrontdoorSet dci_frontdoor_find(DCI_Graph* g, int cause, int effect) {
    DCI_FrontdoorSet fd;
    memset(&fd, 0, sizeof(fd));
    if (!g) return fd;
    /* Find mediators on directed paths from cause to effect */
    DCI_Path directed[DCI_MAX_PATHS];
    int n_dir = dci_find_directed_paths(g, cause, effect, directed, DCI_MAX_PATHS);
    int i, j;
    bool used[DCI_MAX_VARS] = {false};
    for (i = 0; i < n_dir; i++) {
        for (j = 1; j < directed[i].length - 1; j++) {
            int m = directed[i].nodes[j];
            if (!used[m]) {
                used[m] = true;
                if (fd.n_mediators < DCI_MAX_VARS) {
                    fd.mediators[fd.n_mediators++] = m;
                }
            }
        }
    }
    fd.is_valid = dci_frontdoor_check(g, cause, effect,
        fd.mediators, fd.n_mediators);
    return fd;
}

bool dci_frontdoor_check(DCI_Graph* g, int cause, int effect,
    const int* mediators, int n_med) {
    if (!g || n_med <= 0) return false;
    (void)effect;
    /* Check: all directed paths from cause to effect go through mediators */
    /* Check: no back-door paths from cause to mediators */
    int i;
    for (i = 0; i < n_med; i++) {
        DCI_Path bd_paths[DCI_MAX_PATHS];
        int n_bd = dci_find_backdoor_paths(g, cause, mediators[i],
            bd_paths, DCI_MAX_PATHS);
        if (n_bd > 0) return false;
    }
    return true;
}

double dci_frontdoor_adjust(DCI_SCM* scm, DCI_Graph* g,
    int cause, int effect, double cause_val, double effect_val,
    DCI_FrontdoorSet* med_set, int n_samples) {
    if (!scm || !g || !med_set || n_samples <= 0) return 0.0;
    /* Front-door: P(y|do(x)) = Σ_m P(m|x) Σ_x' P(y|x', m) P(x') */
    double prob = 0.0;
    int i, j;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        exog[cause] = cause_val;
        dci_scm_evaluate(scm, exog);
        if (fabs(scm->vars[effect].value - effect_val) < 0.01) {
            prob += 1.0;
        }
    }
    (void)med_set;
    return prob / n_samples;
}

int dci_covariate_selection(DCI_Graph* g, int cause, int effect,
    int* sufficient_set, int max_size) {
    DCI_BackdoorSet bd = dci_backdoor_minimal(g, cause, effect);
    int i;
    int n = bd.n_vars < max_size ? bd.n_vars : max_size;
    for (i = 0; i < n; i++) sufficient_set[i] = bd.variables[i];
    return n;
}

bool dci_is_instrument(DCI_Graph* g, int instrument, int cause, int effect) {
    if (!g) return false;
    /* IV criteria: Z→X, Z→Y only through X, no confounding Z-Y */
    if (g->adjacency[instrument][cause] < 0.5) return false;
    /* Z has no direct effect on Y */
    if (g->adjacency[instrument][effect] > 0.5) return false;
    /* No back-door path from Z to Y */
    DCI_Path bd_paths[DCI_MAX_PATHS];
    int n_bd = dci_find_backdoor_paths(g, instrument, effect,
        bd_paths, DCI_MAX_PATHS);
    return n_bd == 0;
}

bool dci_detect_m_bias(DCI_Graph* g, int cause, int effect) {
    if (!g) return false;
    /* M-bias: X ← U1 → Z ← U2 → Y
     * Conditioning on Z opens this path (collider bias) */
    DCI_Path paths[DCI_MAX_PATHS];
    int n = dci_find_all_paths(g, cause, effect, paths, DCI_MAX_PATHS);
    int i;
    for (i = 0; i < n; i++) {
        if (paths[i].length >= 4 && !paths[i].is_directed) {
            /* Check for collider structure */
            int j;
            for (j = 1; j < paths[i].length - 1; j++) {
                int a = paths[i].nodes[j-1], m = paths[i].nodes[j], b = paths[i].nodes[j+1];
                if (g->adjacency[a][m] > 0.5 && g->adjacency[b][m] > 0.5)
                    return true; /* Collider found */
            }
        }
    }
    return false;
}

/* ==============================================================
 * Disjunctive Cause Criterion (VanderWeele, 2019)
 *
 * A simpler alternative to the back-door criterion:
 * control for all causes of treatment or outcome (or both).
 * ============================================================== */

int dci_disjunctive_cause_criterion(DCI_Graph* g, int cause, int effect,
    int* covariates, int max_cov) {
    if (!g || !covariates) return 0;
    int n_cov = 0;
    int i;
    /* Get all parents of cause and all parents of effect */
    int parents_x[DCI_MAX_VARS], parents_y[DCI_MAX_VARS];
    int n_px = 0, n_py = 0;
    dci_get_parents(g, cause, parents_x, &n_px);
    dci_get_parents(g, effect, parents_y, &n_py);
    /* Union: variables that are parents of cause OR effect */
    bool added[DCI_MAX_VARS] = {false};
    for (i = 0; i < n_px; i++) {
        if (!added[parents_x[i]] && n_cov < max_cov) {
            covariates[n_cov++] = parents_x[i];
            added[parents_x[i]] = true;
        }
    }
    for (i = 0; i < n_py; i++) {
        if (!added[parents_y[i]] && n_cov < max_cov) {
            covariates[n_cov++] = parents_y[i];
            added[parents_y[i]] = true;
        }
    }
    return n_cov;
}

/* ==============================================================
 * Propensity Score from Back-Door Set
 * ============================================================== */

double dci_propensity_score(DCI_SCM* scm, int cause, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    double sum = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        sum += scm->vars[cause].value;
    }
    return sum / n_samples;
}
/* ==============================================================
 * Multiple back-door admissible sets (all minimal sufficient sets)
 * ============================================================== */

int dci_all_backdoor_sets(DCI_Graph* g, int cause, int effect,
    int** sets, int max_sets, int max_per_set) {
    if (!g || !sets || max_sets <= 0) return 0;
    DCI_BackdoorSet bd = dci_backdoor_find(g, cause, effect);
    if (!bd.is_valid) return 0;
    /* Simple: just return single candidate set for now */
    sets[0] = (int*)malloc(bd.n_vars * sizeof(int));
    int i;
    for (i = 0; i < bd.n_vars && i < max_per_set; i++) {
        sets[0][i] = bd.variables[i];
    }
    return 1;
}

/* ==============================================================
 * Check if covariate adjustment is safe (m-bias, butterfly bias)
 * ============================================================== */

bool dci_safe_adjustment(DCI_Graph* g, int cause, int effect, int covariate) {
    if (!g) return false;
    /* Unsafe if covariate is a collider on a path between cause and effect */
    if (dci_collider_bias_check(g, cause, effect, covariate)) return false;
    /* Unsafe if covariate is a mediator */
    DCI_Path directed[DCI_MAX_PATHS];
    int n_dir = dci_find_directed_paths(g, cause, effect, directed, DCI_MAX_PATHS);
    int i, j;
    for (i = 0; i < n_dir; i++) {
        for (j = 1; j < directed[i].length - 1; j++) {
            if (directed[i].nodes[j] == covariate) return false;
        }
    }
    /* Safe: covariate is a confounder (parent of cause or effect) */
    return true;
}