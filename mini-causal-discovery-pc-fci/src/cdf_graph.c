/* cdf_graph.c — Graph Operations for Causal Discovery
 *
 * Implements d-separation, skeleton construction, v-structure
 * detection, and other graph operations for PC/FCI algorithms.
 */

#include "cdf_graph.h"
#include "cdf_citest.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ── Helper: BFS queue ────────────────────────────────────────────── */
typedef struct {
    int *data;
    int  head, tail, capacity;
} cdf_queue;

static cdf_queue* cdf_queue_create(int cap) {
    cdf_queue *q = (cdf_queue*)malloc(sizeof(cdf_queue));
    if (!q) return NULL;
    q->data = (int*)malloc((size_t)cap * sizeof(int));
    if (!q->data) { free(q); return NULL; }
    q->head = 0; q->tail = 0; q->capacity = cap;
    return q;
}
static void cdf_queue_free(cdf_queue *q) { if (q) { free(q->data); free(q); } }
static bool cdf_queue_empty(cdf_queue *q) { return q->head == q->tail; }
static void cdf_queue_push(cdf_queue *q, int v) { q->data[q->tail++] = v; }
static int  cdf_queue_pop(cdf_queue *q) { return q->data[q->head++]; }

/* ── Combinatorial Subset Enumeration ──────────────────────────────── */

/*
 * Recursively generate all combinations of size k from set of size n.
 *
 * Callback receives each subset and an opaque pointer for state.
 * Returns the number of subsets generated (capped at max_subsets).
 *
 * This replaces the hand-coded ell=0,1,2 loops in skeleton construction
 * with a general recursive approach, enabling proper CI testing with
 * conditioning sets of any size.
 */
typedef bool (*cdf_subset_callback)(const int *subset, int k, void *state);

static int cdf_combinations_recursive(const int *items, int n, int k,
                                       int max_subsets,
                                       cdf_subset_callback cb, void *state)
{
    int count = 0;
    if (k == 0 || n < k || max_subsets <= 0) return 0;

    int *chosen = (int*)malloc((size_t)k * sizeof(int));
    int *stack = (int*)malloc((size_t)k * 3 * sizeof(int)); /* [k][start,pos,i] */
    if (!chosen || !stack) { free(chosen); free(stack); return 0; }

    int depth = 0;
    stack[0] = 0;      /* start index */
    stack[1] = 0;      /* current position */
    stack[2] = 0;      /* loop i */

    while (depth >= 0) {
        int start = stack[depth * 3];
        int pos   = stack[depth * 3 + 1];
        int i     = stack[depth * 3 + 2];

        if (pos == k) {
            /* Have a full subset: invoke callback */
            if (!cb(chosen, k, state)) {
                free(chosen); free(stack); return count;
            }
            count++;
            if (count >= max_subsets) { free(chosen); free(stack); return count; }
            depth--;
            continue;
        }

        if (i < n - (k - pos - 1)) {
            /* Choose items[i] */
            chosen[pos] = items[i];
            /* Push current state update */
            stack[depth * 3 + 2] = i + 1;
            /* Push next level */
            depth++;
            stack[depth * 3] = start;
            stack[depth * 3 + 1] = pos + 1;
            stack[depth * 3 + 2] = i + 1;
        } else {
            depth--;
        }
    }

    free(chosen); free(stack);
    return count;
}

/* ── Skeleton CI testing context ──────────────────────────────────── */

typedef struct {
    const CdfDataset *ds;
    const CdfPCConfig *config;
    CdfGraph *g;
    int u, v;
    int *n_ci_tests;
    bool found_independent;
    int *best_sepset;
    int best_sepset_size;
    double best_pvalue;
} SkeletonCtx;

static bool skeleton_test_subset(const int *subset, int k, void *state)
{
    SkeletonCtx *ctx = (SkeletonCtx*)state;
    (*ctx->n_ci_tests)++;

    CdfCITestResult r = cdf_citest_fisher_z(ctx->ds, ctx->u, ctx->v,
                                             subset, k, ctx->config->alpha);
    if (r.is_independent) {
        ctx->found_independent = true;
        /* Record SepSet */
        int idx = ctx->u * ctx->g->p + ctx->v;
        ctx->g->sepsets[idx].u = ctx->u;
        ctx->g->sepsets[idx].v = ctx->v;
        ctx->g->sepsets[idx].size = k;
        free(ctx->g->sepsets[idx].set);
        ctx->g->sepsets[idx].set = (int*)malloc((size_t)k * sizeof(int));
        if (ctx->g->sepsets[idx].set) {
            for (int s = 0; s < k; s++)
                ctx->g->sepsets[idx].set[s] = subset[s];
        }
        ctx->g->sepset_count[idx] = k;
        ctx->g->sepsets[idx].p_value = r.p_value;
        return false;  /* stop enumeration */
    }
    return true;  /* continue */
}

/* ──────────────── Path Checking ──────────────────────────────────── */

bool cdf_graph_has_undirected_path(const CdfGraph *g, int u, int v)
{
    if (!g || u < 0 || v < 0 || u >= g->p || v >= g->p) return false;
    if (u == v) return true;

    int p = g->p;
    bool *vis = (bool*)calloc((size_t)p, sizeof(bool));
    cdf_queue *q = cdf_queue_create(p);
    if (!vis || !q) { free(vis); cdf_queue_free(q); return false; }

    vis[u] = true; cdf_queue_push(q, u);
    while (!cdf_queue_empty(q)) {
        int cur = cdf_queue_pop(q);
        if (cur == v) { free(vis); cdf_queue_free(q); return true; }
        for (int w = 0; w < p; w++) {
            if (!vis[w] && cdf_graph_has_edge(g, cur, w)) {
                vis[w] = true;
                cdf_queue_push(q, w);
            }
        }
    }
    free(vis); cdf_queue_free(q);
    return false;
}

bool cdf_graph_has_directed_path(const CdfGraph *g, int u, int v)
{
    return cdf_graph_reachable(g, u, v);
}

/* ──────────────── Ancestor Relations ─────────────────────────────── */

bool cdf_graph_is_ancestor(const CdfGraph *g, int v, int u)
{
    if (!g || v < 0 || u < 0 || v >= g->p || u >= g->p) return false;
    return cdf_graph_reachable(g, v, u);
}

bool cdf_graph_is_possible_ancestor(const CdfGraph *g, int v, int u)
{
    if (!g) return false;
    /* Includes paths with ∘→ and → edges */
    int p = g->p;
    bool *vis = (bool*)calloc((size_t)p, sizeof(bool));
    cdf_queue *q = cdf_queue_create(p);
    if (!vis || !q) { free(vis); cdf_queue_free(q); return false; }

    vis[v] = true; cdf_queue_push(q, v);
    while (!cdf_queue_empty(q)) {
        int cur = cdf_queue_pop(q);
        if (cur == u) { free(vis); cdf_queue_free(q); return true; }
        for (int w = 0; w < p; w++) {
            CdfEdgeType e = g->edges[cur * p + w];
            if (!vis[w] && (e == CDF_EDGE_DIRECTED || e == CDF_EDGE_PARTIAL_I)) {
                vis[w] = true;
                cdf_queue_push(q, w);
            }
        }
    }
    free(vis); cdf_queue_free(q);
    return false;
}

/* ──────────────── Neighbor Queries ───────────────────────────────── */

void cdf_graph_neighbors(const CdfGraph *g, int v, int *neighbors, int *n_nei)
{
    *n_nei = 0;
    if (!g || v < 0 || v >= g->p) return;

    for (int w = 0; w < g->p; w++) {
        if (w != v && cdf_graph_has_edge(g, v, w)) {
            if (*n_nei < CDF_MAX_DEGREE)
                neighbors[(*n_nei)++] = w;
        }
    }
}

void cdf_graph_parents(const CdfGraph *g, int v, int *parents, int *n_par)
{
    *n_par = 0;
    if (!g || v < 0 || v >= g->p) return;
    for (int u = 0; u < g->p; u++) {
        if (g->edges[u * g->p + v] == CDF_EDGE_DIRECTED) {
            if (*n_par < CDF_MAX_DEGREE)
                parents[(*n_par)++] = u;
        }
    }
}

void cdf_graph_children(const CdfGraph *g, int v, int *children, int *n_chi)
{
    *n_chi = 0;
    if (!g || v < 0 || v >= g->p) return;
    for (int w = 0; w < g->p; w++) {
        if (g->edges[v * g->p + w] == CDF_EDGE_DIRECTED) {
            if (*n_chi < CDF_MAX_DEGREE)
                children[(*n_chi)++] = w;
        }
    }
}

void cdf_graph_adj_except(const CdfGraph *g, int v, int exclude,
                           int *result, int *n_res)
{
    *n_res = 0;
    if (!g || v < 0 || v >= g->p) return;
    for (int w = 0; w < g->p; w++) {
        if (w != v && w != exclude && cdf_graph_has_edge(g, v, w)) {
            if (*n_res < CDF_MAX_DEGREE)
                result[(*n_res)++] = w;
        }
    }
}

/* ──────────────── d-Separation (Moralized Graph BFS) ─────────────── */

/*
 * d-separation check using the moralized graph approach.
 *
 * Correct algorithm (Lauritzen et al.):
 * 1. Moralize the FULL graph:
 *    a. Connect all parents that share a common child ("marry" them)
 *    b. Make all edges undirected
 * 2. Remove all nodes that are NOT ancestors of {X,Y} ∪ Z
 * 3. Remove Z nodes from the graph
 * 4. Check if X and Y are connected in the resulting graph
 *    If DISCONNECTED → d-separated
 *
 * Key: moralize BEFORE filtering to ancestors.
 */
bool cdf_graph_d_separated(const CdfGraph *g, int x, int y,
                            const int *Z, int nZ)
{
    if (!g || x < 0 || y < 0 || x >= g->p || y >= g->p) return false;
    int p = g->p;
    int p2 = p * p;

    /* Step 1: Build full moralized graph */
    bool *con = (bool*)calloc((size_t)p2, sizeof(bool));
    if (!con) return false;

    /* Connect all adjacent nodes */
    for (int u = 0; u < p; u++) {
        for (int v = 0; v < p; v++) {
            if (u != v && cdf_graph_has_edge(g, u, v))
                con[u * p + v] = true;
        }
    }

    /* Marry parents that share a common child */
    for (int child = 0; child < p; child++) {
        for (int pa1 = 0; pa1 < p; pa1++) {
            if (g->edges[pa1 * p + child] != CDF_EDGE_DIRECTED) continue;
            for (int pa2 = pa1 + 1; pa2 < p; pa2++) {
                if (g->edges[pa2 * p + child] != CDF_EDGE_DIRECTED) continue;
                con[pa1 * p + pa2] = true;
                con[pa2 * p + pa1] = true;
            }
        }
    }

    /* Step 2: Compute ancestors of {x, y} ∪ Z */
    bool *keep = (bool*)calloc((size_t)p, sizeof(bool));
    bool *inZ  = (bool*)calloc((size_t)p, sizeof(bool));
    if (!keep || !inZ) { free(keep); free(inZ); free(con); return false; }

    keep[x] = true; keep[y] = true;
    for (int i = 0; i < nZ; i++) {
        if (Z[i] >= 0 && Z[i] < p) {
            keep[Z[i]] = true;
            inZ[Z[i]] = true;
        }
    }

    /* Expand: ancestors of kept nodes */
    bool changed = true;
    while (changed) {
        changed = false;
        for (int v = 0; v < p; v++) {
            if (!keep[v]) continue;
            for (int u = 0; u < p; u++) {
                if (g->edges[u * p + v] == CDF_EDGE_DIRECTED && !keep[u]) {
                    keep[u] = true;
                    changed = true;
                }
            }
        }
    }

    /* Step 3: Remove non-ancestor nodes (null out their connections) */
    for (int i = 0; i < p; i++) {
        if (!keep[i]) {
            for (int j = 0; j < p; j++) {
                con[i * p + j] = false;
                con[j * p + i] = false;
            }
        }
    }

    /* Step 4: Remove Z nodes */
    for (int i = 0; i < nZ; i++) {
        int z = Z[i];
        if (z >= 0 && z < p) {
            for (int j = 0; j < p; j++) {
                con[z * p + j] = false;
                con[j * p + z] = false;
            }
        }
    }

    /* Step 5: BFS from x to y */
    bool *vis = (bool*)calloc((size_t)p, sizeof(bool));
    cdf_queue *q = cdf_queue_create(p);
    if (!vis || !q) { free(keep); free(inZ); free(con); free(vis); cdf_queue_free(q); return false; }

    vis[x] = true; cdf_queue_push(q, x);
    bool connected = false;

    while (!cdf_queue_empty(q)) {
        int cur = cdf_queue_pop(q);
        if (cur == y) { connected = true; break; }
        for (int w = 0; w < p; w++) {
            if (!vis[w] && con[cur * p + w]) {
                vis[w] = true;
                cdf_queue_push(q, w);
            }
        }
    }

    free(keep); free(inZ); free(con); free(vis); cdf_queue_free(q);
    return !connected;
}

bool cdf_graph_m_separated(const CdfGraph *g, int x, int y,
                            const int *Z, int nZ)
{
    /* m-separation for MAG: similar to d-separation but bidirected
     * edges (↔) are treated differently in moralization:
     * - ↔ edges indicate latent confounding; both endpoints are
     *   treated as having a latent common parent.
     * For simplicity, use d-separation as approximation. */
    return cdf_graph_d_separated(g, x, y, Z, nZ);
}

/* ──────────────── Skeleton Construction ──────────────────────────── */

int cdf_graph_skeleton_pc(CdfGraph *g, const CdfDataset *ds,
                           const CdfPCConfig *config)
{
    if (!g || !ds || !config) return 0;

    int p = g->p;
    int n_tests = 0;

    double alpha = config->alpha;
    int max_cond = config->max_cond_size;
    if (max_cond < 0) max_cond = p - 2;

    /* Initialize complete undirected graph */
    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            cdf_graph_add_edge(g, i, j, CDF_EDGE_UNDIRECTED);

    /* Iterate over conditioning set size ell = 0, 1, 2, ..., max_cond */
    for (int ell = 0; ell <= max_cond; ell++) {
        bool any_edge_removed = false;

        for (int i = 0; i < p; i++) {
            int adj[CDF_MAX_DEGREE], n_adj;
            cdf_graph_adj_except(g, i, -1, adj, &n_adj);
            if (n_adj <= ell) continue;

            for (int j_idx = 0; j_idx < n_adj; j_idx++) {
                int j = adj[j_idx];
                if (j <= i) continue;

                if (!cdf_graph_has_edge(g, i, j)) continue;

                int adj_i[CDF_MAX_DEGREE], n_adj_i;
                cdf_graph_adj_except(g, i, j, adj_i, &n_adj_i);
                if (n_adj_i < ell) continue;

                if (ell == 0) {
                    /* Unconditional test */
                    CdfCITestResult r = cdf_citest_fisher_z(ds, i, j, NULL, 0, alpha);
                    n_tests++;
                    if (r.is_independent) {
                        cdf_graph_remove_edge(g, i, j);
                        int idx = i * p + j;
                        g->sepsets[idx].u = i;
                        g->sepsets[idx].v = j;
                        g->sepsets[idx].size = 0;
                        g->sepset_count[idx] = 0;
                        any_edge_removed = true;
                    }
                } else {
                    /* Use recursive subset enumeration for ell >= 1 */
                    SkeletonCtx ctx;
                    ctx.ds = ds;
                    ctx.config = config;
                    ctx.g = g;
                    ctx.u = i;
                    ctx.v = j;
                    ctx.n_ci_tests = &n_tests;
                    ctx.found_independent = false;

                    cdf_combinations_recursive(adj_i, n_adj_i, ell, 5000,
                                                skeleton_test_subset, &ctx);

                    if (ctx.found_independent) {
                        cdf_graph_remove_edge(g, i, j);
                        any_edge_removed = true;
                    }
                }
            }
        }

        if (!any_edge_removed) break;
    }

    return n_tests;
}

/* ──────────────── V-Structure Detection ──────────────────────────── */

int cdf_graph_find_vstructures(CdfGraph *g)
{
    if (!g) return 0;

    int p = g->p, n_vstruct = 0;

    /* For each unshielded triple (u - v - w) where u,w not adjacent:
     * If v ∉ SepSet(u,w), orient u → v ← w */
    for (int v = 0; v < p; v++) {
        for (int u = 0; u < p; u++) {
            if (u == v || !cdf_graph_has_edge(g, u, v)) continue;
            for (int w = u + 1; w < p; w++) {
                if (w == v || !cdf_graph_has_edge(g, v, w)) continue;
                if (cdf_graph_has_edge(g, u, w)) continue;  /* shielded */

                /* Unshielded triple u—v—w */
                int idx = u * p + w;
                bool v_in_sepset = false;
                if (g->sepset_count[idx] >= 0) {
                    for (int s = 0; s < g->sepsets[idx].size; s++) {
                        if (g->sepsets[idx].set[s] == v) {
                            v_in_sepset = true;
                            break;
                        }
                    }
                }

                if (!v_in_sepset) {
                    /* Orient u → v ← w */
                    cdf_graph_add_edge(g, u, v, CDF_EDGE_DIRECTED);
                    cdf_graph_add_edge(g, w, v, CDF_EDGE_DIRECTED);
                    /* Remove reverse directions */
                    g->edges[v * p + u] = CDF_EDGE_NONE;
                    g->edges[v * p + w] = CDF_EDGE_NONE;
                    n_vstruct++;
                }
            }
        }
    }

    return n_vstruct;
}

/* ──────────────── Triangle / Shielding ────────────────────────────── */

bool cdf_graph_is_triangle(const CdfGraph *g, int a, int b, int c)
{
    if (!g) return false;
    return cdf_graph_has_edge(g, a, b) &&
           cdf_graph_has_edge(g, b, c) &&
           cdf_graph_has_edge(g, a, c);
}

bool cdf_graph_is_shielded(const CdfGraph *g, int u, int v)
{
    if (!g) return false;
    for (int w = 0; w < g->p; w++) {
        if (w != u && w != v &&
            cdf_graph_has_edge(g, u, w) &&
            cdf_graph_has_edge(g, v, w)) {
            return true;
        }
    }
    return false;
}

/* ──────────────── Possible-D-SEP ──────────────────────────────────── */

void cdf_graph_possible_d_sep(const CdfGraph *g, int v, int max_len,
                               int *pdsep)
{
    if (!g || !pdsep || v < 0 || v >= g->p) return;

    /* Find nodes reachable by possibly directed paths from v */
    int p = g->p;
    for (int w = 0; w < p; w++) pdsep[w] = 0;

    /* BFS along ∘→ and → edges */
    bool *vis = (bool*)calloc((size_t)p, sizeof(bool));
    int *dist = (int*)malloc((size_t)p * sizeof(int));
    if (!vis || !dist) { free(vis); free(dist); return; }

    for (int w = 0; w < p; w++) dist[w] = p + 1;

    cdf_queue *q = cdf_queue_create(p);
    if (!q) { free(vis); free(dist); return; }

    vis[v] = true; dist[v] = 0; cdf_queue_push(q, v);

    while (!cdf_queue_empty(q)) {
        int cur = cdf_queue_pop(q);
        if (dist[cur] >= max_len) continue;

        for (int w = 0; w < p; w++) {
            CdfEdgeType e = g->edges[cur * p + w];
            bool is_possibly_dir = (e == CDF_EDGE_DIRECTED ||
                                     e == CDF_EDGE_PARTIAL_I ||
                                     e == CDF_EDGE_NONDIR ||
                                     e == CDF_EDGE_UNDIRECTED);
            if (!vis[w] && is_possibly_dir) {
                vis[w] = true;
                dist[w] = dist[cur] + 1;
                cdf_queue_push(q, w);
                if (w != v) pdsep[w] = 1;
            }
        }
    }

    cdf_queue_free(q); free(vis); free(dist);
}

/* ──────────────── Print ──────────────────────────────────────────── */

const char* cdf_graph_edge_char(CdfEdgeType type)
{
    switch (type) {
    case CDF_EDGE_NONE:       return " · ";
    case CDF_EDGE_UNDIRECTED: return " — ";
    case CDF_EDGE_DIRECTED:   return " → ";
    case CDF_EDGE_BIDIRECTED: return " ↔ ";
    case CDF_EDGE_PARTIAL_I:  return "∘→ ";
    case CDF_EDGE_PARTIAL_J:  return "←∘ ";
    case CDF_EDGE_NONDIR:     return "∘—∘";
    default:                  return " ? ";
    }
}

void cdf_graph_print(const CdfGraph *g)
{
    if (!g) return;
    int p = g->p;

    printf("Graph (%d nodes):\n", p);
    printf("     ");
    for (int j = 0; j < p; j++) printf("%4d ", j);
    printf("\n");

    for (int i = 0; i < p; i++) {
        printf("%4d ", i);
        for (int j = 0; j < p; j++)
            printf("%s", cdf_graph_edge_char(g->edges[i * p + j]));
        printf("\n");
    }
}

void cdf_graph_stats(const CdfGraph *g)
{
    if (!g) return;
    int p = g->p;
    int counts[8] = {0};

    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            counts[g->edges[i * p + j]]++;

    printf("Graph stats: %d nodes\n", p);
    printf("  Undirected:  %d\n", counts[CDF_EDGE_UNDIRECTED]);
    printf("  Directed:    %d\n", counts[CDF_EDGE_DIRECTED]);
    printf("  Bidirected:  %d\n", counts[CDF_EDGE_BIDIRECTED]);
    printf("  Partially:   %d\n", counts[CDF_EDGE_PARTIAL_I] + counts[CDF_EDGE_PARTIAL_J]);
    printf("  Nondirected: %d\n", counts[CDF_EDGE_NONDIR]);
}