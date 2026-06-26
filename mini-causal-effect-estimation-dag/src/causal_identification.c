/*
 * causal_identification.c ！ Causal Effect Identification
 *
 * Implements Pearl's causal identification criteria:
 *   - Back-door criterion (Pearl, 1993)
 *   - Front-door criterion (Pearl, 1995)
 *   - do-calculus Rule 2 (Pearl, 1995)
 *   - Complete identification via back-door/front-door search
 *
 * A causal effect P(y | do(x)) is identifiable from observational data
 * iff there exists a set of covariates Z satisfying the back-door
 * criterion, or a mediator M satisfying the front-door criterion,
 * or the do-calculus rules can transform it into an observable expression.
 *
 * References:
 *   Pearl, "Causality", 2nd ed, 2009, Ch.3 (The Effects of Interventions)
 *   Pearl, "Causal diagrams for empirical research", Biometrika, 1995
 *   Shpitser & Pearl, "Complete identification of causal effects", UAI 2006
 */

#include "causal_identification.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ================================================================
 * Back-Door Criterion
 *
 * A set Z satisfies the back-door criterion for (X, Y) if:
 *   1. No node in Z is a descendant of X
 *   2. Z blocks every back-door path from X to Y
 *
 * Equivalently: (Y ? X | Z) in the graph G_{X} (graph with outgoing
 * edges from X removed).
 * ================================================================ */

bool causal_backdoor_criterion(const DAG* dag, int x, int y,
                                const int *z, int n_z) {
    if (!dag || x < 0 || x >= dag->n || y < 0 || y >= dag->n)
        return false;

    /* Condition 1: No node in Z is a descendant of X */
    int desc[DAG_MAX_NODES];
    int nd = dag_descendants(dag, x, desc, DAG_MAX_NODES);
    if (nd < 0) nd = 0;

    /* X itself is a descendant of X, but that's fine ！ we check for Z nodes */
    for (int i = 0; i < n_z; i++) {
        if (z[i] == x) return false;  /* X cannot be in its own adjustment set */
        for (int j = 0; j < nd; j++) {
            if (desc[j] == z[i]) return false;
        }
    }

    /* Condition 2: In G_{X} (remove outgoing edges from X), Y d-sep X | Z */
    DAG* g_x = dag_do_intervention(dag, x);
    if (!g_x) {
        /* Fallback: create graph without outgoing edges from x */
        g_x = dag_create(dag->n);
        if (!g_x) return false;
        for (int i = 0; i < dag->n; i++)
            for (int j = 0; j < dag->n; j++)
                if (dag->adj[i * dag->n + j])
                    dag_add_edge(g_x, i, j);
        /* Remove outgoing from x */
        for (int j = 0; j < dag->n; j++)
            dag_remove_edge(g_x, x, j);
    }

    /* Actually, the back-door check uses the original graph but checks
     * if Z blocks all back-door paths. We verify using d-separation
     * in the original graph: Y ? X | Z in the absence of confounding
     * (i.e., check all back-door paths are blocked). */

    /* Find all back-door paths and verify each is blocked by Z */
    Path* bd_paths[DAG_MAX_PATHS];
    int n_bd = dag_find_backdoor_paths(dag, x, y, bd_paths, DAG_MAX_PATHS);

    if (n_bd < 0) { dag_free(g_x); return false; }

    /* If no back-door paths exist, back-door criterion is trivially satisfied
     * (as long as no Z node is a descendant of X). */
    bool criterion_holds = true;

    for (int p = 0; p < n_bd && criterion_holds; p++) {
        /* Check if this back-door path is blocked by Z */
        /* A path is blocked if any non-collider node on it is in Z */
        bool path_blocked = false;
        for (int k = 0; k < bd_paths[p]->length && !path_blocked; k++) {
            int node = bd_paths[p]->nodes[k];
            /* Skip X and Y ！ they are endpoints */
            if (node == x || node == y) continue;

            /* Check if node is a collider on this path */
            bool is_coll = false;
            if (k > 0 && k < bd_paths[p]->length - 1) {
                int a = bd_paths[p]->nodes[k - 1];
                int b = node;
                int c = bd_paths[p]->nodes[k + 1];
                is_coll = (dag_has_edge(dag, a, b) && dag_has_edge(dag, c, b));
            }

            /* Non-collider in Z blocks the path */
            if (!is_coll) {
                for (int zi = 0; zi < n_z; zi++) {
                    if (z[zi] == node) { path_blocked = true; break; }
                }
            }
        }

        if (!path_blocked) {
            criterion_holds = false;
        }

        /* Free individual path */
        if (bd_paths[p]) { free(bd_paths[p]->nodes); free(bd_paths[p]); }
    }

    dag_free(g_x);
    return criterion_holds;
}

/* ================================================================
 * Front-Door Criterion
 *
 * A variable M satisfies the front-door criterion for (X, Y) if:
 *   1. M intercepts all directed paths from X to Y
 *   2. There is no unblocked back-door path from X to M
 *   3. All back-door paths from M to Y are blocked by X
 *
 * Then: P(y | do(x)) = sum_m P(m | x) * sum_{x'} P(y | x', m) * P(x')
 * ================================================================ */

bool causal_frontdoor_criterion(const DAG* dag, int x, int y, int m) {
    if (!dag || x < 0 || x >= dag->n || y < 0 || y >= dag->n ||
        m < 0 || m >= dag->n || m == x || m == y)
        return false;

    /* Condition 1: All directed paths from X to Y go through M */
    /* Check: is there a directed path from X to Y that does NOT pass through M? */
    /* First check if a path from X to Y exists */
    if (!dag_has_path(dag, x, y)) return false;

    /* Build graph without M to check if X can still reach Y */
    DAG* g_no_m = dag_create(dag->n);
    if (!g_no_m) return false;
    for (int i = 0; i < dag->n; i++)
        for (int j = 0; j < dag->n; j++)
            if (dag->adj[i * dag->n + j] && i != m && j != m)
                dag_add_edge(g_no_m, i, j);

    bool path_without_m = dag_has_path(g_no_m, x, y);
    dag_free(g_no_m);
    if (path_without_m) return false;  /* M doesn't intercept all paths */

    /* Condition 2: No unblocked back-door path from X to M */
    /* Z = empty set ★ check if there are any unblocked back-door paths */
    Path* bd_xm[DAG_MAX_PATHS];
    int n_bd_xm = dag_find_backdoor_paths(dag, x, m, bd_xm, DAG_MAX_PATHS);
    if (n_bd_xm > 0) {
        for (int i = 0; i < n_bd_xm; i++) {
            if (bd_xm[i]->nodes) free(bd_xm[i]->nodes);
            free(bd_xm[i]);
        }
        return false;  /* There are open back-door paths from X to M */
    }

    /* Condition 3: All back-door paths from M to Y are blocked by X */
    Path* bd_my[DAG_MAX_PATHS];
    int n_bd_my = dag_find_backdoor_paths(dag, m, y, bd_my, DAG_MAX_PATHS);

    if (n_bd_my < 0) return false;

    bool all_blocked_by_x = true;
    for (int p = 0; p < n_bd_my && all_blocked_by_x; p++) {
        bool path_blocked = false;
        for (int k = 0; k < bd_my[p]->length && !path_blocked; k++) {
            int node = bd_my[p]->nodes[k];
            if (node == m || node == y) continue;
            if (node == x) {
                /* Check if X is a non-collider on this path */
                bool is_coll = false;
                if (k > 0 && k < bd_my[p]->length - 1) {
                    int a = bd_my[p]->nodes[k - 1];
                    int c = bd_my[p]->nodes[k + 1];
                    is_coll = (dag_has_edge(dag, a, node) &&
                               dag_has_edge(dag, c, node));
                }
                if (!is_coll) path_blocked = true;
            }
        }
        if (!path_blocked) all_blocked_by_x = false;
        if (bd_my[p]->nodes) free(bd_my[p]->nodes);
        free(bd_my[p]);
    }

    return all_blocked_by_x;
}

/* ================================================================
 * do-Calculus Rule 2 (Action/Observation Exchange)
 *
 * P(y | do(x), z, w) = P(y | x, z, w)
 * if Y is d-separated from X by {Z, W} in the graph G_{X}
 * (graph with incoming edges to X removed).
 *
 * This is the key rule that allows replacing do(x) with conditioning
 * on x when the right independence holds.
 * ================================================================ */

bool do_calculus_rule2(const DAG* dag, int x, int y,
                        const int *z, int n_z, const int *w, int n_w) {
    if (!dag || x < 0 || x >= dag->n || y < 0 || y >= dag->n)
        return false;

    /* Create the graph G_{X_} with incoming edges to X removed */
    DAG* g_x_underline = dag_do_intervention(dag, x);
    if (!g_x_underline) return false;

    /* Actually, do_calculus_rule2 uses G_{X_} (remove outgoing from X).
     * Wait ！ Rule 2: P(y|do(x),z) = P(y|x,z) if Y ? X | Z in G_{X_}
     * where G_{X_} removes edges OUT of X.
     *
     * Let me re-read: Pearl's Rule 2 states:
     * P(y | do(x), z) = P(y | x, z) if (Y ? X | Z) in G_{\underline{X}}
     * where G_{\underline{X}} is the graph with arrows emerging from X deleted.
     *
     * So we need: DAG with outgoing edges from X removed. */

    /* Rebuild g_x_underline: remove outgoing edges from x */
    dag_free(g_x_underline);
    g_x_underline = dag_create(dag->n);
    if (!g_x_underline) return false;
    for (int i = 0; i < dag->n; i++)
        for (int j = 0; j < dag->n; j++)
            if (dag->adj[i * dag->n + j])
                dag_add_edge(g_x_underline, i, j);

    /* Remove outgoing edges from x */
    for (int j = 0; j < dag->n; j++)
        dag_remove_edge(g_x_underline, x, j);

    /* Combine Z and W into a single conditioning set */
    int combined[DAG_MAX_NODES];
    int n_comb = 0;
    for (int i = 0; i < n_z && n_comb < DAG_MAX_NODES; i++)
        combined[n_comb++] = z[i];
    for (int i = 0; i < n_w && n_comb < DAG_MAX_NODES; i++)
        combined[n_comb++] = w[i];

    /* Test d-separation: Y ? X | (Z, W) in G_{\underline{X}} */
    DSepResult* ds = dag_d_separation(g_x_underline, x, y, combined, n_comb);
    dag_free(g_x_underline);

    bool result = false;
    if (ds) {
        result = ds->is_separated;
        if (ds->blocking_path) {
            if (ds->blocking_path->nodes) free(ds->blocking_path->nodes);
            free(ds->blocking_path);
        }
        free(ds);
    }
    return result;
}

/* ================================================================
 * Causal Identify: determine if effect of X on Y is identifiable
 *
 * The algorithm tries:
 *   1. Check front-door criterion
 *   2. Find a back-door adjustment set
 *   3. If neither works, the effect may still be identifiable via
 *      do-calculus (we return "not identified" for complex cases)
 * ================================================================ */

IdentificationResult* causal_identify(const DAG* dag, int x, int y) {
    if (!dag || x < 0 || x >= dag->n || y < 0 || y >= dag->n)
        return NULL;

    IdentificationResult* r = calloc(1, sizeof(IdentificationResult));
    if (!r) return NULL;

    r->adjustment_set = calloc((size_t)DAG_MAX_NODES, sizeof(int));
    if (!r->adjustment_set) { free(r); return NULL; }

    /* Strategy 1: Search for front-door mediator */
    for (int m = 0; m < dag->n; m++) {
        if (m == x || m == y) continue;
        if (causal_frontdoor_criterion(dag, x, y, m)) {
            r->identifiable = true;
            r->method = ID_METHOD_FRONTDOOR;
            /* Front-door: no direct adjustment set, use M as mediator */
            r->adjustment_set[0] = m;
            r->n_adjust = 1;
            r->explanation = calloc(256, sizeof(char));
            if (r->explanation)
                snprintf(r->explanation, 256,
                         "Front-door criterion: mediator V%d intercepts all "
                         "causal paths from V%d to V%d", m, x, y);
            return r;
        }
    }

    /* Strategy 2: Find back-door adjustment set */
    int adj_found = causal_find_backdoor_adjustment(dag, x, y,
                        r->adjustment_set, DAG_MAX_NODES);
    if (adj_found >= 0) {
        r->identifiable = true;
        r->method = ID_METHOD_BACKDOOR;
        r->n_adjust = adj_found;
        r->explanation = calloc(256, sizeof(char));
        if (r->explanation)
            snprintf(r->explanation, 256,
                     "Back-door criterion: adjustment set of size %d found", adj_found);
        return r;
    }

    /* Strategy 3: Try do-calculus Rule 2 with parents */
    /* (Simplified ！ full do-calculus requires the ID algorithm) */
    int parents_x[DAG_MAX_NODES];
    int np_x = 0;
    for (int p = 0; p < dag->n; p++)
        if (dag_has_edge(dag, p, x) && np_x < DAG_MAX_NODES)
            parents_x[np_x++] = p;

    if (np_x > 0) {
        /* Check Rule 2: Y ? X | parents(X) in G_{\underline{X}} */
        if (do_calculus_rule2(dag, x, y, parents_x, np_x, NULL, 0)) {
            r->identifiable = true;
            r->method = ID_METHOD_DO_CALC;
            r->n_adjust = np_x;
            for (int i = 0; i < np_x; i++) r->adjustment_set[i] = parents_x[i];
            r->explanation = calloc(256, sizeof(char));
            if (r->explanation)
                snprintf(r->explanation, 256,
                         "do-calculus Rule 2: effect identified via parents of V%d", x);
            return r;
        }
    }

    /* Not identifiable via back-door or front-door */
    r->identifiable = false;
    r->method = ID_METHOD_NONE;
    r->n_adjust = 0;
    r->explanation = calloc(256, sizeof(char));
    if (r->explanation)
        snprintf(r->explanation, 256,
                 "Effect of V%d on V%d is not identifiable via "
                 "back-door or front-door from this DAG", x, y);
    return r;
}

void ident_result_free(IdentificationResult* r) {
    if (!r) return;
    free(r->adjustment_set);
    free(r->explanation);
    free(r);
}

/* ================================================================
 * Find Back-Door Adjustment Set
 *
 * Enumerate subsets of non-descendants of X and test if they
 * satisfy the back-door criterion. Return the smallest valid set.
 * This is an exhaustive search over subsets of potential confounders.
 * ================================================================ */

int causal_find_backdoor_adjustment(const DAG* dag, int x, int y,
                                     int *adjustment_set, int max_n) {
    if (!dag || !adjustment_set || max_n <= 0 ||
        x < 0 || x >= dag->n || y < 0 || y >= dag->n) return -1;

    /* Build list of candidate variables: all non-descendants of X */
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

    /* Search subsets of candidates by increasing size */
    for (int sz = 0; sz <= n_cand && sz <= max_n; sz++) {
        /* Generate all subsets of size sz */
        /* Use combinatorial enumeration via bitmask */
        unsigned long long limit = 1ULL << (unsigned long long)n_cand;
        for (unsigned long long mask = 0; mask < limit; mask++) {
            /* Count bits */
            int bits = 0;
            unsigned long long m = mask;
            while (m) { bits += (int)(m & 1); m >>= 1; }
            if (bits != sz) continue;

            /* Build the candidate set */
            int z_buf[DAG_MAX_NODES];
            int n_z = 0;
            for (int b = 0; b < n_cand && n_z < DAG_MAX_NODES; b++) {
                if (mask & (1ULL << (unsigned long long)b))
                    z_buf[n_z++] = candidates[b];
            }

            /* Test if this set satisfies back-door criterion */
            if (causal_backdoor_criterion(dag, x, y, z_buf, n_z)) {
                int copy_n = (n_z < max_n) ? n_z : max_n;
                for (int i = 0; i < copy_n; i++)
                    adjustment_set[i] = z_buf[i];
                return copy_n;
            }
        }
    }

    return -1;  /* No valid back-door adjustment set found */
}
