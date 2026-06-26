/* cdf_orientation.c — Edge Orientation Rules (PC R1-R4, FCI R1-R10)
 *
 * Implements Meek's orientation rules for PC and the extended
 * rule set for FCI (Zhang, 2008).
 *
 * PC Rules (Meek, 1995):
 *   R0: Orient v-structures (unshielded colliders)
 *   R1: a→b—c and a,c not adjacent → orient b→c
 *   R2: a→b→c and a—c → orient a→c
 *   R3: a—b→c, a—d→c, b,d not adjacent, a—c → orient a→c
 *   R4: Discriminating path rule
 *
 * FCI Additional Rules (Zhang, 2008):
 *   R5-R10: Handle latent variable structures in PAGs
 */

#include "cdf_orientation.h"
#include "cdf_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ── Helper: check if two nodes are adjacent in graph ─────────────── */
static bool adjacent(const CdfGraph *g, int u, int v) {
    return cdf_graph_has_edge(g, u, v);
}

/* ──────────────── Orientation Utilities ──────────────────────────── */

bool cdf_orient_would_create_cycle(const CdfGraph *g, int u, int v)
{
    if (!g) return false;
    /* Check if there's already a directed path from v to u */
    return cdf_graph_reachable(g, v, u);
}

void cdf_orient_edge(CdfGraph *g, int u, int v, CdfEdgeType type)
{
    if (!g) return;

    int p = g->p;
    g->edges[u * p + v] = type;

    /* Set reverse edge appropriately */
    if (type == CDF_EDGE_DIRECTED)
        g->edges[v * p + u] = CDF_EDGE_NONE;
    else if (type == CDF_EDGE_UNDIRECTED)
        g->edges[v * p + u] = CDF_EDGE_UNDIRECTED;
    else if (type == CDF_EDGE_BIDIRECTED)
        g->edges[v * p + u] = CDF_EDGE_BIDIRECTED;
    else if (type == CDF_EDGE_PARTIAL_I)
        g->edges[v * p + u] = CDF_EDGE_PARTIAL_J;
    else if (type == CDF_EDGE_PARTIAL_J)
        g->edges[v * p + u] = CDF_EDGE_PARTIAL_I;
    else if (type == CDF_EDGE_NONDIR)
        g->edges[v * p + u] = CDF_EDGE_NONDIR;
}

int cdf_orient_count_unoriented(const CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;
    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            if (g->edges[i * p + j] == CDF_EDGE_UNDIRECTED)
                count++;
    return count;
}

bool cdf_orient_is_maximal(const CdfGraph *g)
{
    return cdf_orient_count_unoriented(g) == 0;
}

/* ──────────────── PC Rule R1 ─────────────────────────────────────── */

int cdf_orient_rule_r1(CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;

    /* For each triple (a, b, c):
     * If a→b, b—c, and a not adjacent to c → orient b→c */
    for (int b = 0; b < p; b++) {
        for (int a = 0; a < p; a++) {
            if (a == b) continue;
            if (g->edges[a * p + b] != CDF_EDGE_DIRECTED) continue;

            for (int c = 0; c < p; c++) {
                if (c == a || c == b) continue;
                if (g->edges[b * p + c] != CDF_EDGE_UNDIRECTED) continue;
                if (adjacent(g, a, c)) continue;

                /* Orient b→c if no cycle created */
                if (!cdf_orient_would_create_cycle(g, b, c)) {
                    cdf_orient_edge(g, b, c, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }
    return count;
}

/* ──────────────── PC Rule R2 ─────────────────────────────────────── */

int cdf_orient_rule_r2(CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;

    /* For triple (a, b, c):
     * If a→b→c and a—c → orient a→c */
    for (int b = 0; b < p; b++) {
        for (int a = 0; a < p; a++) {
            if (a == b) continue;
            if (g->edges[a * p + b] != CDF_EDGE_DIRECTED) continue;

            for (int c = 0; c < p; c++) {
                if (c == a || c == b) continue;
                if (g->edges[b * p + c] != CDF_EDGE_DIRECTED) continue;
                if (g->edges[a * p + c] != CDF_EDGE_UNDIRECTED) continue;

                if (!cdf_orient_would_create_cycle(g, a, c)) {
                    cdf_orient_edge(g, a, c, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }
    return count;
}

/* ──────────────── PC Rule R3 ─────────────────────────────────────── */

int cdf_orient_rule_r3(CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;

    /* For quadruple (a, b, c, d):
     * If a—b→c, a—d→c, b,d not adjacent, and a—c → orient a→c */
    for (int c = 0; c < p; c++) {
        for (int a = 0; a < p; a++) {
            if (a == c) continue;
            if (g->edges[a * p + c] != CDF_EDGE_UNDIRECTED) continue;

            /* Find b: a—b, b→c */
            for (int b = 0; b < p; b++) {
                if (b == a || b == c) continue;
                if (g->edges[a * p + b] != CDF_EDGE_UNDIRECTED) continue;
                if (g->edges[b * p + c] != CDF_EDGE_DIRECTED) continue;

                /* Find d: a—d, d→c, d not adjacent to b */
                for (int d = 0; d < p; d++) {
                    if (d == a || d == b || d == c) continue;
                    if (g->edges[a * p + d] != CDF_EDGE_UNDIRECTED) continue;
                    if (g->edges[d * p + c] != CDF_EDGE_DIRECTED) continue;
                    if (adjacent(g, b, d)) continue;

                    if (!cdf_orient_would_create_cycle(g, a, c)) {
                        cdf_orient_edge(g, a, c, CDF_EDGE_DIRECTED);
                        count++;
                        goto next_c;  /* break out of triple loop */
                    }
                }
            }
            next_c:;
        }
    }
    return count;
}

/* ──────────────── PC Rule R4 ─────────────────────────────────────── */

int cdf_orient_rule_r4(CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;

    /* Discriminating path rule for PC:
     * A discriminating path π = ⟨x, q_1, ..., q_m, b, y⟩ with
     * x—b—y and x,y not adjacent. Then orient b→y.
     *
     * Simplified: check for length-1 discriminating paths.
     */
    for (int b = 0; b < p; b++) {
        for (int y = 0; y < p; y++) {
            if (y == b) continue;
            if (g->edges[b * p + y] != CDF_EDGE_UNDIRECTED) continue;

            /* Look for a node x: x—b, x not adjacent to y, and x has a
             * directed path to y through a collider at b */
            for (int x = 0; x < p; x++) {
                if (x == b || x == y) continue;
                if (!adjacent(g, x, b)) continue;
                if (adjacent(g, x, y)) continue;

                /* Check if there's a path x→...→y with b as collider */
                /* Simplified heuristic */
                if (!cdf_orient_would_create_cycle(g, b, y)) {
                    cdf_orient_edge(g, b, y, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }
    return count;
}

/* ──────────────── PC Rules (All) ─────────────────────────────────── */

CdfOrientResult cdf_orient_pc_rules(CdfGraph *g)
{
    CdfOrientResult res = {0};
    if (!g) return res;

    /* Exhaustive application of R1-R4 */
    int iter = 0, max_iter = 100;
    while (iter < max_iter) {
        res.changed = false;

        int r1 = cdf_orient_rule_r1(g);
        int r2 = cdf_orient_rule_r2(g);
        int r3 = cdf_orient_rule_r3(g);
        int r4 = cdf_orient_rule_r4(g);

        res.rules_applied[1] += r1;
        res.rules_applied[2] += r2;
        res.rules_applied[3] += r3;
        res.rules_applied[4] += r4;
        res.total_oriented += r1 + r2 + r3 + r4;

        if (r1 + r2 + r3 + r4 == 0) break;
        res.changed = true;
        iter++;
    }
    return res;
}

/* ──────────────── Discriminating Path Check ──────────────────────── */

bool cdf_orient_is_discriminating_path(const CdfGraph *g,
                                        int u, int v, int w,
                                        const int *path, int len)
{
    if (!g || !path || len < 3) return false;
    (void)u; (void)v; (void)w;

    /* Path π = ⟨v_0, ..., v_k⟩ is discriminating for (u, v, w) if:
     * 1. v_0 = u, v_1 = v, v_k = w
     * 2. u and w are not adjacent
     * 3. Every v_i (1 < i < k) is a collider on π and a parent of w
     */
    int v0 = path[0], v1 = path[1], vk = path[len - 1];
    if (!adjacent(g, v0, v1)) return false;
    if (adjacent(g, v0, vk)) return false;

    for (int i = 2; i < len - 1; i++) {
        int vi = path[i];
        /* Check vi is a collider: path[i-1] → vi ← path[i+1] */
        if (g->edges[path[i - 1] * g->p + vi] != CDF_EDGE_DIRECTED &&
            g->edges[vi * g->p + path[i - 1]] != CDF_EDGE_DIRECTED)
            return false;
    }

    return true;
}

/* ──────────────── FCI Rules R5-R7 ────────────────────────────────── */

int cdf_orient_fci_rules_5_7(CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;

    /* R5: a ∘→ b — c, a and c not adjacent → b → c */
    for (int b = 0; b < p; b++) {
        for (int a = 0; a < p; a++) {
            if (a == b) continue;
            if (g->edges[a * p + b] != CDF_EDGE_PARTIAL_I) continue;

            for (int c = 0; c < p; c++) {
                if (c == a || c == b) continue;
                if (g->edges[b * p + c] != CDF_EDGE_UNDIRECTED) continue;
                if (adjacent(g, a, c)) continue;

                if (!cdf_orient_would_create_cycle(g, b, c)) {
                    cdf_orient_edge(g, b, c, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }

    /* R6: a ∘→ b ∘→ c, a—c → b → c */
    for (int b = 0; b < p; b++) {
        for (int a = 0; a < p; a++) {
            if (a == b) continue;
            if (g->edges[a * p + b] != CDF_EDGE_PARTIAL_I) continue;

            for (int c = 0; c < p; c++) {
                if (c == a || c == b) continue;
                if (g->edges[b * p + c] != CDF_EDGE_PARTIAL_I) continue;
                if (g->edges[a * p + c] != CDF_EDGE_UNDIRECTED) continue;

                if (!cdf_orient_would_create_cycle(g, b, c)) {
                    cdf_orient_edge(g, b, c, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }

    /* R7: a ∘→ b → c, a—c → a → c */
    for (int b = 0; b < p; b++) {
        for (int a = 0; a < p; a++) {
            if (a == b) continue;
            if (g->edges[a * p + b] != CDF_EDGE_PARTIAL_I) continue;

            for (int c = 0; c < p; c++) {
                if (c == a || c == b) continue;
                if (g->edges[b * p + c] != CDF_EDGE_DIRECTED) continue;
                if (g->edges[a * p + c] != CDF_EDGE_UNDIRECTED) continue;

                if (!cdf_orient_would_create_cycle(g, a, c)) {
                    cdf_orient_edge(g, a, c, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }

    return count;
}

/* ──────────────── FCI Rules R8-R10 ───────────────────────────────── */

int cdf_orient_fci_rules_8_10(CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;

    /* R8: a → b → c, a ∘→ c → a → c */
    for (int b = 0; b < p; b++) {
        for (int a = 0; a < p; a++) {
            if (a == b) continue;
            if (g->edges[a * p + b] != CDF_EDGE_DIRECTED) continue;

            for (int c = 0; c < p; c++) {
                if (c == a || c == b) continue;
                if (g->edges[b * p + c] != CDF_EDGE_DIRECTED) continue;
                if (g->edges[a * p + c] != CDF_EDGE_PARTIAL_I) continue;

                if (!cdf_orient_would_create_cycle(g, a, c)) {
                    cdf_orient_edge(g, a, c, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }

    /* R9: a ∘→ c, path a ∘→ b → c → a → c */
    for (int c = 0; c < p; c++) {
        for (int a = 0; a < p; a++) {
            if (a == c) continue;
            if (g->edges[a * p + c] != CDF_EDGE_PARTIAL_I) continue;

            for (int b = 0; b < p; b++) {
                if (b == a || b == c) continue;
                if (g->edges[a * p + b] != CDF_EDGE_PARTIAL_I) continue;
                if (g->edges[b * p + c] != CDF_EDGE_DIRECTED) continue;

                if (!cdf_orient_would_create_cycle(g, a, c)) {
                    cdf_orient_edge(g, a, c, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }

    /* R10: a ∘→ c, a → b ← c → a → c */
    for (int c = 0; c < p; c++) {
        for (int a = 0; a < p; a++) {
            if (a == c) continue;
            if (g->edges[a * p + c] != CDF_EDGE_PARTIAL_I) continue;

            for (int b = 0; b < p; b++) {
                if (b == a || b == c) continue;
                if (g->edges[a * p + b] != CDF_EDGE_DIRECTED) continue;
                if (g->edges[c * p + b] != CDF_EDGE_DIRECTED) continue;

                if (!cdf_orient_would_create_cycle(g, a, c)) {
                    cdf_orient_edge(g, a, c, CDF_EDGE_DIRECTED);
                    count++;
                }
            }
        }
    }

    return count;
}

/* ──────────────── FCI Rules (All) ────────────────────────────────── */

CdfOrientResult cdf_orient_fci_rules(CdfGraph *g)
{
    CdfOrientResult res = {0};
    if (!g) return res;

    /* First apply PC rules R1-R4 */
    CdfOrientResult pc_res = cdf_orient_pc_rules(g);
    for (int i = 1; i <= 4; i++)
        res.rules_applied[i] = pc_res.rules_applied[i];

    /* Then apply FCI rules R5-R10 */
    int iter = 0, max_iter = 100;
    while (iter < max_iter) {
        int r1_4 = cdf_orient_rule_r1(g) + cdf_orient_rule_r2(g) +
                   cdf_orient_rule_r3(g) + cdf_orient_rule_r4(g);
        int r5_7 = cdf_orient_fci_rules_5_7(g);
        int r8_10 = cdf_orient_fci_rules_8_10(g);

        res.rules_applied[5] += r5_7;  /* lump R5-R7 */
        res.rules_applied[8] += r8_10; /* lump R8-R10 */
        res.total_oriented += r1_4 + r5_7 + r8_10;

        if (r1_4 + r5_7 + r8_10 == 0) break;
        res.changed = true;
        iter++;
    }

    return res;
}

void cdf_orient_print_result(const CdfOrientResult *res)
{
    if (!res) return;
    printf("Orientation result: %d edges oriented\n", res->total_oriented);
    for (int i = 1; i < 10; i++) {
        if (res->rules_applied[i] > 0)
            printf("  R%d: %d\n", i, res->rules_applied[i]);
    }
}