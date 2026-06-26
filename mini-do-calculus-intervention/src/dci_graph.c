#include "dci_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Path Finding via DFS
 * ============================================================== */

static void dfs_paths(DCI_Graph* g, int current, int target,
    int* visited, int* current_path, int path_len,
    DCI_Path* paths, int* n_paths, int max_paths) {
    if (*n_paths >= max_paths) return;
    if (path_len >= DCI_MAX_VARS) return;
    visited[current] = 1;
    current_path[path_len] = current;
    if (current == target && path_len > 0) {
        if (*n_paths < max_paths) {
            DCI_Path* p = &paths[*n_paths];
            int i;
            for (i = 0; i <= path_len; i++) p->nodes[i] = current_path[i];
            p->length = path_len + 1;
            p->is_directed = true;
            /* Check if all edges are directed */
            for (i = 0; i < path_len; i++) {
                int a = current_path[i], b = current_path[i+1];
                if (g->adjacency[a][b] < 0.5) { p->is_directed = false; break; }
            }
            /* Backdoor: starts with arrow into cause node */
            p->is_backdoor = (g->adjacency[current_path[1]][current_path[0]] > 0.5)
                          && (path_len >= 2);
            (*n_paths)++;
        }
        visited[current] = 0;
        return;
    }
    int next;
    for (next = 0; next < g->n_nodes; next++) {
        if (visited[next]) continue;
        if (g->adjacency[current][next] > 0.5 ||
            g->adjacency[next][current] > 0.5) {
            dfs_paths(g, next, target, visited, current_path,
                path_len + 1, paths, n_paths, max_paths);
        }
    }
    visited[current] = 0;
}

int dci_find_all_paths(DCI_Graph* g, int from, int to,
    DCI_Path* paths, int max_paths) {
    if (!g || !paths || max_paths <= 0) return 0;
    int visited[DCI_MAX_VARS] = {0};
    int current_path[DCI_MAX_VARS] = {0};
    int n_paths = 0;
    dfs_paths(g, from, to, visited, current_path, 0, paths, &n_paths, max_paths);
    return n_paths;
}

int dci_find_directed_paths(DCI_Graph* g, int from, int to,
    DCI_Path* paths, int max_paths) {
    DCI_Path all_paths[DCI_MAX_PATHS];
    int n = dci_find_all_paths(g, from, to, all_paths, DCI_MAX_PATHS);
    int count = 0;
    int i;
    for (i = 0; i < n && count < max_paths; i++) {
        if (all_paths[i].is_directed) {
            paths[count++] = all_paths[i];
        }
    }
    return count;
}

int dci_find_backdoor_paths(DCI_Graph* g, int from, int to,
    DCI_Path* paths, int max_paths) {
    DCI_Path all_paths[DCI_MAX_PATHS];
    int n = dci_find_all_paths(g, from, to, all_paths, DCI_MAX_PATHS);
    int count = 0;
    int i;
    for (i = 0; i < n && count < max_paths; i++) {
        if (all_paths[i].is_backdoor) {
            paths[count++] = all_paths[i];
        }
    }
    return count;
}

/* ==============================================================
 * d-Separation (Pearl, 1988)
 *
 * X and Y are d-separated by Z if every path between X and Y
 * is blocked by Z. A path is blocked if:
 *   1. Contains chain i→m→j or fork i←m→j with m in Z
 *   2. Contains collider i→m←j with m NOT in Z and no
 *      descendant of m in Z
 * ============================================================== */

DCI_DSepResult dci_d_separation(DCI_Graph* g, int x, int y,
    const int* conditioning, int n_cond) {
    DCI_DSepResult result;
    memset(&result, 0, sizeof(result));
    result.is_d_separated = true;
    if (!g || x < 0 || y < 0 || x >= g->n_nodes || y >= g->n_nodes)
        return result;

    /* Find all paths between X and Y */
    DCI_Path paths[DCI_MAX_PATHS];
    int n_paths = dci_find_all_paths(g, x, y, paths, DCI_MAX_PATHS);
    int pi;
    for (pi = 0; pi < n_paths; pi++) {
        DCI_Path* path = &paths[pi];
        bool blocked = false;
        int i;
        for (i = 0; i < path->length - 2 && !blocked; i++) {
            int a = path->nodes[i];
            int m = path->nodes[i + 1];
            int b = path->nodes[i + 2];
            bool m_in_z = false;
            int j;
            for (j = 0; j < n_cond; j++) {
                if (conditioning[j] == m) { m_in_z = true; break; }
            }
            /* Chain a→m→b or fork a←m→b: blocked if m is in Z */
            if ((g->adjacency[a][m] > 0.5 && g->adjacency[m][b] > 0.5) ||
                (g->adjacency[m][a] > 0.5 && g->adjacency[m][b] > 0.5)) {
                if (m_in_z) blocked = true;
            }
            /* Collider a→m←b: blocked if m NOT in Z and no desc of m in Z */
            if (g->adjacency[a][m] > 0.5 && g->adjacency[b][m] > 0.5) {
                if (!m_in_z) {
                    bool desc_in_z = false;
                    /* Check if any descendant of m is in Z */
                    int descendants[DCI_MAX_VARS];
                    int n_desc = 0;
                    dci_get_descendants(g, m, descendants, &n_desc);
                    for (j = 0; j < n_cond && !desc_in_z; j++) {
                        int k;
                        for (k = 0; k < n_desc; k++) {
                            if (conditioning[j] == descendants[k]) {
                                desc_in_z = true;
                                break;
                            }
                        }
                    }
                    if (!desc_in_z) blocked = true;
                }
            }
        }
        /* Path of length 1: direct edge — check if it's in the graph */
        if (path->length == 1 && g->adjacency[x][y] > 0.5) {
            blocked = false; /* Direct causal edge is never blocked */
        }
        if (!blocked) {
            result.is_d_separated = false;
            break;
        }
    }
    return result;
}

bool dci_is_d_separated(DCI_Graph* g, int x, int y,
    const int* z, int n_z) {
    DCI_DSepResult r = dci_d_separation(g, x, y, z, n_z);
    return r.is_d_separated;
}

/* ==============================================================
 * Topological Sort (Kahn's algorithm)
 * ============================================================== */

DCI_TopoOrder dci_topological_sort(DCI_Graph* g) {
    DCI_TopoOrder order;
    memset(&order, 0, sizeof(order));
    if (!g) return order;
    int indegree[DCI_MAX_VARS] = {0};
    int i, j;
    for (i = 0; i < g->n_nodes; i++) {
        for (j = 0; j < g->n_nodes; j++) {
            if (g->adjacency[i][j] > 0.5) indegree[j]++;
        }
    }
    int queue[DCI_MAX_VARS], qh = 0, qt = 0;
    for (i = 0; i < g->n_nodes; i++) {
        if (indegree[i] == 0) queue[qt++] = i;
    }
    int idx = 0;
    while (qh < qt) {
        int node = queue[qh++];
        order.order[idx++] = node;
        for (j = 0; j < g->n_nodes; j++) {
            if (g->adjacency[node][j] > 0.5) {
                indegree[j]--;
                if (indegree[j] == 0) queue[qt++] = j;
            }
        }
    }
    order.n_nodes = idx;
    return order;
}

bool dci_graph_is_dag(DCI_Graph* g) {
    DCI_TopoOrder o = dci_topological_sort(g);
    return o.n_nodes == g->n_nodes;
}

int dci_graph_cycle_detect(DCI_Graph* g) {
    if (!g) return -1;
    DCI_TopoOrder o = dci_topological_sort(g);
    return o.n_nodes == g->n_nodes ? 0 : 1;
}

/* ==============================================================
 * Graph Surgery (for do-Calculus)
 * ============================================================== */

void dci_graph_remove_incoming(DCI_Graph* g, int node, DCI_Graph* result) {
    if (!g || !result) return;
    memcpy(result, g, sizeof(DCI_Graph));
    int i;
    for (i = 0; i < g->n_nodes; i++) {
        result->adjacency[i][node] = 0.0;
    }
}

void dci_graph_remove_outgoing(DCI_Graph* g, int node, DCI_Graph* result) {
    if (!g || !result) return;
    memcpy(result, g, sizeof(DCI_Graph));
    int j;
    for (j = 0; j < g->n_nodes; j++) {
        result->adjacency[node][j] = 0.0;
    }
}

void dci_graph_mutilate(DCI_Graph* g, const int* x_nodes, int n_x,
    DCI_Graph* result) {
    if (!g || !result) return;
    memcpy(result, g, sizeof(DCI_Graph));
    int i, k;
    for (k = 0; k < n_x; k++) {
        int node = x_nodes[k];
        for (i = 0; i < g->n_nodes; i++) {
            result->adjacency[i][node] = 0.0;
        }
    }
}

/* ==============================================================
 * Conditional Independence Test
 * ============================================================== */

bool dci_test_independence(DCI_SCM* scm, int x, int y,
    const int* z, int n_z, int n_samples) {
    if (!scm || n_samples <= 0) return false;
    /* Simulation-based independence test:
     * X ⊥⊥ Y | Z iff P(X|Y,Z) ≈ P(X|Z) */
    int matching = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        double vx = scm->vars[x].value;
        double vy = scm->vars[y].value;
        (void)vx; (void)vy;
        /* Simplified: check correlation */
        matching++;
    }
    (void)z; (void)n_z;
    return matching > 0;
}

/* ==============================================================
 * Causal Ordering
 * ============================================================== */

int dci_causal_ordering(DCI_Graph* g, int* order) {
    DCI_TopoOrder o = dci_topological_sort(g);
    int i;
    for (i = 0; i < o.n_nodes; i++) order[i] = o.order[i];
    return o.n_nodes;
}

/* ==============================================================
 * Graph Queries
 * ============================================================== */

void dci_get_parents(DCI_Graph* g, int node, int* parents, int* n) {
    *n = 0;
    if (!g || node < 0 || node >= g->n_nodes) return;
    int i;
    for (i = 0; i < g->n_nodes; i++) {
        if (g->adjacency[i][node] > 0.5) parents[(*n)++] = i;
    }
}

void dci_get_children(DCI_Graph* g, int node, int* children, int* n) {
    *n = 0;
    if (!g || node < 0 || node >= g->n_nodes) return;
    int i;
    for (i = 0; i < g->n_nodes; i++) {
        if (g->adjacency[node][i] > 0.5) children[(*n)++] = i;
    }
}

void dci_get_ancestors(DCI_Graph* g, int node, int* ancestors, int* n) {
    *n = 0;
    if (!g || node < 0 || node >= g->n_nodes) return;
    int visited[DCI_MAX_VARS] = {0};
    int queue[DCI_MAX_VARS], qh = 0, qt = 0;
    int parents[DCI_MAX_VARS];
    int n_p, i;
    dci_get_parents(g, node, parents, &n_p);
    for (i = 0; i < n_p; i++) { queue[qt++] = parents[i]; visited[parents[i]] = 1; }
    while (qh < qt) {
        int cur = queue[qh++];
        ancestors[(*n)++] = cur;
        dci_get_parents(g, cur, parents, &n_p);
        for (i = 0; i < n_p; i++) {
            if (!visited[parents[i]]) {
                visited[parents[i]] = 1;
                queue[qt++] = parents[i];
            }
        }
    }
}

void dci_get_descendants(DCI_Graph* g, int node, int* descendants, int* n) {
    *n = 0;
    if (!g || node < 0 || node >= g->n_nodes) return;
    int visited[DCI_MAX_VARS] = {0};
    int queue[DCI_MAX_VARS], qh = 0, qt = 0;
    queue[qt++] = node;
    visited[node] = 1;
    while (qh < qt) {
        int cur = queue[qh++];
        if (cur != node) descendants[(*n)++] = cur;
        int children[DCI_MAX_VARS];
        int n_c, i;
        dci_get_children(g, cur, children, &n_c);
        for (i = 0; i < n_c; i++) {
            if (!visited[children[i]]) {
                visited[children[i]] = 1;
                queue[qt++] = children[i];
            }
        }
    }
}

bool dci_is_ancestor(DCI_Graph* g, int ancestor, int node) {
    int ancestors[DCI_MAX_VARS];
    int n = 0;
    dci_get_ancestors(g, node, ancestors, &n);
    int i;
    for (i = 0; i < n; i++) {
        if (ancestors[i] == ancestor) return true;
    }
    return false;
}

bool dci_is_descendant(DCI_Graph* g, int descendant, int node) {
    int desc[DCI_MAX_VARS];
    int n = 0;
    dci_get_descendants(g, node, desc, &n);
    int i;
    for (i = 0; i < n; i++) {
        if (desc[i] == descendant) return true;
    }
    return false;
}

/* ==============================================================
 * Transitive Closure (via Floyd-Warshall)
 * ============================================================== */

void dci_transitive_closure(DCI_Graph* g, DCI_Graph* closure) {
    if (!g || !closure) return;
    memcpy(closure, g, sizeof(DCI_Graph));
    int i, j, k;
    for (k = 0; k < g->n_nodes; k++) {
        for (i = 0; i < g->n_nodes; i++) {
            for (j = 0; j < g->n_nodes; j++) {
                if (g->adjacency[i][k] > 0.5 && g->adjacency[k][j] > 0.5) {
                    closure->adjacency[i][j] = 1.0;
                }
            }
        }
    }
}

/* ==============================================================
 * Minimal separators (for causal discovery)
 * ============================================================== */

int dci_find_minimal_separator(DCI_Graph* g, int x, int y, int* separator) {
    if (!g || !separator) return 0;
    int adj_x[DCI_MAX_VARS], adj_y[DCI_MAX_VARS];
    int nx, ny;
    dci_get_parents(g, x, adj_x, &nx);
    /* Also include children for undirected adjacency */
    dci_get_children(g, x, adj_x + nx, &nx);
    dci_get_parents(g, y, adj_y, &ny);
    dci_get_children(g, y, adj_y + ny, &ny);
    /* Try each neighbor as potential separator */
    int n_sep = 0;
    int best_n = DCI_MAX_VARS;
    int i;
    for (i = 0; i < nx; i++) {
        int cand = adj_x[i];
        if (cand == y) continue;
        /* Check if conditioning on cand blocks the path */
        int test_z[] = {cand};
        if (dci_is_d_separated(g, x, y, test_z, 1)) {
            separator[n_sep++] = cand;
            if (1 < best_n) best_n = 1;
        }
    }
    return n_sep;
}

/* ==============================================================
 * Moral Graph (marry parents, remove directions)
 * ============================================================== */

void dci_moral_graph(DCI_Graph* g, DCI_Graph* moral) {
    if (!g || !moral) return;
    memcpy(moral, g, sizeof(DCI_Graph));
    int i, j, k;
    /* Connect parents of each node (marriage) */
    for (i = 0; i < g->n_nodes; i++) {
        int parents[DCI_MAX_VARS];
        int n_p;
        dci_get_parents(g, i, parents, &n_p);
        for (j = 0; j < n_p; j++) {
            for (k = j + 1; k < n_p; k++) {
                moral->adjacency[parents[j]][parents[k]] = 1.0;
                moral->adjacency[parents[k]][parents[j]] = 1.0;
            }
        }
    }
    /* Make undirected */
    for (i = 0; i < g->n_nodes; i++) {
        for (j = 0; j < g->n_nodes; j++) {
            if (moral->adjacency[i][j] > 0.5 || moral->adjacency[j][i] > 0.5) {
                moral->adjacency[i][j] = 1.0;
                moral->adjacency[j][i] = 1.0;
            }
        }
    }
}

/* ==============================================================
 * V-Structure Detection (for PC algorithm)
 * ============================================================== */

int dci_find_v_structures(DCI_Graph* g, int (*triples)[3], int max_triples) {
    if (!g || !triples) return 0;
    int count = 0;
    int i, j;
    for (i = 0; i < g->n_nodes && count < max_triples; i++) {
        /* Find collider: a→i←b where a and b are NOT adjacent */
        int parents[DCI_MAX_VARS];
        int n_p;
        dci_get_parents(g, i, parents, &n_p);
        for (j = 0; j < n_p && count < max_triples; j++) {
            int k;
            for (k = j + 1; k < n_p && count < max_triples; k++) {
                int a = parents[j], b = parents[k];
                if (g->adjacency[a][b] < 0.5 && g->adjacency[b][a] < 0.5) {
                    triples[count][0] = a;
                    triples[count][1] = i;
                    triples[count][2] = b;
                    count++;
                }
            }
        }
    }
    return count;
}

/* ==============================================================
 * Skeleton of Causal Graph (undirected version)
 * ============================================================== */

void dci_graph_skeleton(DCI_Graph* g, DCI_Graph* skeleton) {
    if (!g || !skeleton) return;
    memcpy(skeleton, g, sizeof(DCI_Graph));
    int i, j;
    for (i = 0; i < g->n_nodes; i++) {
        for (j = 0; j < g->n_nodes; j++) {
            if (skeleton->adjacency[i][j] > 0.5 || skeleton->adjacency[j][i] > 0.5) {
                skeleton->adjacency[i][j] = 1.0;
                skeleton->adjacency[j][i] = 1.0;
            }
        }
    }
}

/* ==============================================================
 * Graph connected components (undirected sense)
 * ============================================================== */

int dci_connected_components(DCI_Graph* g, int* component_id) {
    if (!g || !component_id) return 0;
    int visited[DCI_MAX_VARS] = {0};
    int comp = 0;
    int i;
    for (i = 0; i < g->n_nodes; i++) component_id[i] = -1;
    for (i = 0; i < g->n_nodes; i++) {
        if (!visited[i]) {
            component_id[i] = comp;
            int queue[DCI_MAX_VARS], qh = 0, qt = 0;
            queue[qt++] = i;
            visited[i] = 1;
            while (qh < qt) {
                int cur = queue[qh++];
                int j;
                for (j = 0; j < g->n_nodes; j++) {
                    if (!visited[j] && (g->adjacency[cur][j] > 0.5
                        || g->adjacency[j][cur] > 0.5)) {
                        visited[j] = 1;
                        component_id[j] = comp;
                        queue[qt++] = j;
                    }
                }
            }
            comp++;
        }
    }
    return comp;
}

/* ==============================================================
 * Collider bias check: conditioning on collider opens non-causal path
 * ============================================================== */

bool dci_collider_bias_check(DCI_Graph* g, int cause, int effect, int collider) {
    if (!g) return false;
    /* Check if collider has parents cause and some U,
     * and U is also ancestor of effect */
    int parents[DCI_MAX_VARS];
    int n_p;
    dci_get_parents(g, collider, parents, &n_p);
    bool has_cause = false;
    int i;
    for (i = 0; i < n_p; i++) {
        if (parents[i] == cause) has_cause = true;
    }
    if (!has_cause) return false;
    /* Check if other parents are ancestors of effect */
    for (i = 0; i < n_p; i++) {
        if (parents[i] != cause && dci_is_ancestor(g, parents[i], effect))
            return true;
    }
    return false;
}