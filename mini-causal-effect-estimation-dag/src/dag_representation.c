/*
 * dag_representation.c ！ DAG Data Structures and Graph Algorithms
 *
 * Implements the fundamental DAG operations for causal inference:
 * adjacency matrix representation, topological sort (Kahn's algorithm),
 * path finding (BFS/DFS), cycle detection, d-separation (moral graph method),
 * back-door path enumeration, and graph surgery for interventions.
 *
 * References:
 *   Pearl, "Causality: Models, Reasoning, and Inference", 2nd ed, 2009, Ch.1-2
 *   Koller & Friedman, "Probabilistic Graphical Models", 2009, Ch.3
 *   Bang-Jensen & Gutin, "Digraphs: Theory, Algorithms and Applications", 2009
 */

#include "dag_representation.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ================================================================
 * DAG Lifecycle
 * ================================================================ */

DAG* dag_create(int n) {
    if (n <= 0 || n > DAG_MAX_NODES) return NULL;
    DAG* dag = calloc(1, sizeof(DAG));
    if (!dag) return NULL;
    dag->n = n;
    dag->adj        = calloc((size_t)n * n, sizeof(int));
    dag->parents    = calloc((size_t)n * n, sizeof(int));
    dag->children   = calloc((size_t)n * n, sizeof(int));
    dag->n_parents  = calloc((size_t)n, sizeof(int));
    dag->n_children = calloc((size_t)n, sizeof(int));
    dag->topo_order = calloc((size_t)n, sizeof(int));
    dag->names      = calloc((size_t)n, sizeof(char*));
    if (!dag->adj || !dag->parents || !dag->children ||
        !dag->n_parents || !dag->n_children || !dag->topo_order || !dag->names) {
        dag_free(dag);
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        dag->names[i] = calloc(DAG_NAME_LEN, sizeof(char));
        if (dag->names[i]) snprintf(dag->names[i], DAG_NAME_LEN, "V%d", i);
        dag->topo_order[i] = -1;
    }
    return dag;
}

void dag_free(DAG* dag) {
    if (!dag) return;
    if (dag->names) {
        for (int i = 0; i < dag->n; i++) free(dag->names[i]);
        free(dag->names);
    }
    free(dag->adj);
    free(dag->parents);
    free(dag->children);
    free(dag->n_parents);
    free(dag->n_children);
    free(dag->topo_order);
    free(dag);
}

/* ================================================================
 * Edge Operations
 * ================================================================ */

int dag_add_edge(DAG* dag, int i, int j) {
    if (!dag || i < 0 || i >= dag->n || j < 0 || j >= dag->n || i == j)
        return -1;
    if (dag->adj[i * dag->n + j]) return 0;

    dag->adj[i * dag->n + j] = 1;

    /* Cycle check via DFS from j back to i */
    int *visited = calloc((size_t)dag->n, sizeof(int));
    int *stack = calloc((size_t)dag->n, sizeof(int));
    if (!visited || !stack) {
        dag->adj[i * dag->n + j] = 0;
        free(visited); free(stack);
        return -1;
    }
    int top = 0;
    stack[top++] = j;
    visited[j] = 1;
    bool creates_cycle = false;
    while (top > 0) {
        int v = stack[--top];
        if (v == i) { creates_cycle = true; break; }
        for (int k = 0; k < dag->n; k++) {
            if (dag->adj[v * dag->n + k] && !visited[k]) {
                visited[k] = 1;
                stack[top++] = k;
            }
        }
    }
    free(visited); free(stack);

    if (creates_cycle) {
        dag->adj[i * dag->n + j] = 0;
        return -1;
    }

    int pc = dag->n_children[i];
    dag->children[i * dag->n + pc] = j;
    dag->n_children[i]++;
    int pp = dag->n_parents[j];
    dag->parents[j * dag->n + pp] = i;
    dag->n_parents[j]++;
    dag->topo_order[0] = -1;
    return 0;
}

int dag_remove_edge(DAG* dag, int i, int j) {
    if (!dag || i < 0 || i >= dag->n || j < 0 || j >= dag->n) return -1;
    if (!dag->adj[i * dag->n + j]) return -1;
    dag->adj[i * dag->n + j] = 0;

    int nc = dag->n_children[i];
    for (int k = 0; k < nc; k++) {
        if (dag->children[i * dag->n + k] == j) {
            dag->children[i * dag->n + k] = dag->children[i * dag->n + nc - 1];
            dag->n_children[i]--;
            break;
        }
    }
    int np = dag->n_parents[j];
    for (int k = 0; k < np; k++) {
        if (dag->parents[j * dag->n + k] == i) {
            dag->parents[j * dag->n + k] = dag->parents[j * dag->n + np - 1];
            dag->n_parents[j]--;
            break;
        }
    }
    dag->topo_order[0] = -1;
    return 0;
}

bool dag_has_edge(const DAG* dag, int i, int j) {
    if (!dag || i < 0 || i >= dag->n || j < 0 || j >= dag->n) return false;
    return dag->adj[i * dag->n + j] != 0;
}

/* ================================================================
 * Topological Sort ！ Kahn's Algorithm O(V+E)
 * ================================================================ */

int dag_topological_sort(DAG* dag) {
    if (!dag) return -1;
    int *in_degree = calloc((size_t)dag->n, sizeof(int));
    int *queue = calloc((size_t)dag->n, sizeof(int));
    if (!in_degree || !queue) { free(in_degree); free(queue); return -1; }

    for (int i = 0; i < dag->n; i++)
        for (int j = 0; j < dag->n; j++)
            if (dag->adj[i * dag->n + j]) in_degree[j]++;

    int qh = 0, qt = 0;
    for (int i = 0; i < dag->n; i++)
        if (in_degree[i] == 0) queue[qt++] = i;

    int idx = 0;
    while (qh < qt) {
        int u = queue[qh++];
        dag->topo_order[idx++] = u;
        for (int v = 0; v < dag->n; v++) {
            if (dag->adj[u * dag->n + v]) {
                in_degree[v]--;
                if (in_degree[v] == 0) queue[qt++] = v;
            }
        }
    }
    free(in_degree); free(queue);

    if (idx != dag->n) { dag->topo_order[0] = -1; return -1; }
    return 0;
}

/* ================================================================
 * Cycle Detection ！ DFS with 3-coloring O(V+E)
 * ================================================================ */

bool dag_is_acyclic(const DAG* dag) {
    if (!dag) return false;
    int *color = calloc((size_t)dag->n, sizeof(int));
    int *stack = calloc((size_t)dag->n * 2, sizeof(int));
    if (!color || !stack) { free(color); free(stack); return false; }

    bool has_cycle = false;
    for (int start = 0; start < dag->n && !has_cycle; start++) {
        if (color[start] != 0) continue;
        int top = 0;
        stack[top * 2] = start;
        stack[top * 2 + 1] = 0;
        top++;

        while (top > 0) {
            top--;
            int v = stack[top * 2];
            int state = stack[top * 2 + 1];
            if (state == 0) {
                if (color[v] == 1) { has_cycle = true; break; }
                if (color[v] == 2) continue;
                color[v] = 1;
                stack[top * 2] = v;
                stack[top * 2 + 1] = 1;
                top++;
                for (int w = dag->n - 1; w >= 0; w--) {
                    if (dag->adj[v * dag->n + w]) {
                        if (color[w] == 1) { has_cycle = true; break; }
                        if (color[w] == 0) {
                            stack[top * 2] = w;
                            stack[top * 2 + 1] = 0;
                            top++;
                        }
                    }
                }
            } else {
                color[v] = 2;
            }
        }
    }
    free(color); free(stack);
    return !has_cycle;
}

/* ================================================================
 * Path Finding ！ BFS O(V+E)
 * ================================================================ */

bool dag_has_path(const DAG* dag, int from, int to) {
    if (!dag || from < 0 || from >= dag->n || to < 0 || to >= dag->n)
        return false;
    if (from == to) return true;

    int *visited = calloc((size_t)dag->n, sizeof(int));
    int *queue = calloc((size_t)dag->n, sizeof(int));
    if (!visited || !queue) { free(visited); free(queue); return false; }

    int qh = 0, qt = 0;
    queue[qt++] = from;
    visited[from] = 1;
    bool found = false;

    while (qh < qt) {
        int v = queue[qh++];
        for (int w = 0; w < dag->n; w++) {
            if (dag->adj[v * dag->n + w] && !visited[w]) {
                if (w == to) { found = true; break; }
                visited[w] = 1;
                queue[qt++] = w;
            }
        }
        if (found) break;
    }
    free(visited); free(queue);
    return found;
}

/* ================================================================
 * Ancestor/Descendant Relations ！ Reverse/Forward BFS
 * ================================================================ */

int dag_ancestors(const DAG* dag, int node, int *ancestors, int max_n) {
    if (!dag || !ancestors || node < 0 || node >= dag->n || max_n <= 0)
        return -1;

    int *visited = calloc((size_t)dag->n, sizeof(int));
    int *queue = calloc((size_t)dag->n, sizeof(int));
    if (!visited || !queue) { free(visited); free(queue); return -1; }

    int qh = 0, qt = 0;
    queue[qt++] = node;
    visited[node] = 1;
    int count = 0;

    while (qh < qt && count < max_n) {
        int v = queue[qh++];
        if (v != node) ancestors[count++] = v;
        for (int u = 0; u < dag->n; u++) {
            if (dag->adj[u * dag->n + v] && !visited[u]) {
                visited[u] = 1;
                queue[qt++] = u;
            }
        }
    }
    free(visited); free(queue);
    return count;
}

int dag_descendants(const DAG* dag, int node, int *descendants, int max_n) {
    if (!dag || !descendants || node < 0 || node >= dag->n || max_n <= 0)
        return -1;

    int *visited = calloc((size_t)dag->n, sizeof(int));
    int *queue = calloc((size_t)dag->n, sizeof(int));
    if (!visited || !queue) { free(visited); free(queue); return -1; }

    int qh = 0, qt = 0;
    queue[qt++] = node;
    visited[node] = 1;
    int count = 0;

    while (qh < qt && count < max_n) {
        int v = queue[qh++];
        if (v != node) descendants[count++] = v;
        for (int w = 0; w < dag->n; w++) {
            if (dag->adj[v * dag->n + w] && !visited[w]) {
                visited[w] = 1;
                queue[qt++] = w;
            }
        }
    }
    free(visited); free(queue);
    return count;
}

/* ================================================================
 * Directed Path Search ！ DFS with parent tracking
 * ================================================================ */

Path* dag_find_directed_path(const DAG* dag, int from, int to) {
    if (!dag || from < 0 || from >= dag->n || to < 0 || to >= dag->n)
        return NULL;
    if (from == to) {
        Path* p = calloc(1, sizeof(Path));
        if (!p) return NULL;
        p->nodes = calloc(1, sizeof(int));
        if (!p->nodes) { free(p); return NULL; }
        p->nodes[0] = from;
        p->length = 1;
        p->is_directed = true;
        return p;
    }

    int *visited = calloc((size_t)dag->n, sizeof(int));
    int *parent = calloc((size_t)dag->n, sizeof(int));
    int *stack = calloc((size_t)dag->n, sizeof(int));
    if (!visited || !parent || !stack) {
        free(visited); free(parent); free(stack); return NULL;
    }
    for (int i = 0; i < dag->n; i++) parent[i] = -1;

    int top = 0;
    stack[top++] = from;
    visited[from] = 1;
    bool found = false;

    while (top > 0) {
        int v = stack[--top];
        if (v == to) { found = true; break; }
        for (int w = 0; w < dag->n; w++) {
            if (dag->adj[v * dag->n + w] && !visited[w]) {
                visited[w] = 1;
                parent[w] = v;
                stack[top++] = w;
            }
        }
    }

    if (!found) {
        free(visited); free(parent); free(stack);
        return NULL;
    }

    int rev[DAG_MAX_PATH_LEN];
    int len = 0;
    for (int v = to; v != -1 && len < DAG_MAX_PATH_LEN; v = parent[v])
        rev[len++] = v;

    Path* p = calloc(1, sizeof(Path));
    if (!p) { free(visited); free(parent); free(stack); return NULL; }
    p->nodes = calloc((size_t)len, sizeof(int));
    if (!p->nodes) { free(p); free(visited); free(parent); free(stack); return NULL; }
    for (int i = 0; i < len; i++) p->nodes[i] = rev[len - 1 - i];
    p->length = len;
    p->is_directed = true;

    free(visited); free(parent); free(stack);
    return p;
}

/* ================================================================
 * Back-Door Path Enumeration
 *
 * A back-door path between X and Y is any undirected path that:
 *   1. Starts with an edge INTO X (X <- ...)
 *   2. Is not blocked by a collider unconditionally
 *
 * We find all simple paths starting with X <- edge up to length limit.
 * ================================================================ */

static bool is_collider_triple(const DAG* dag, int a, int b, int c) {
    return (dag->adj[a * dag->n + b] && dag->adj[c * dag->n + b]);
}

int dag_find_backdoor_paths(const DAG* dag, int from, int to,
                              Path **paths, int max_paths) {
    if (!dag || !paths || from < 0 || from >= dag->n ||
        to < 0 || to >= dag->n || max_paths <= 0) return -1;

    int *undirected = calloc((size_t)dag->n * dag->n, sizeof(int));
    if (!undirected) return -1;
    for (int i = 0; i < dag->n; i++)
        for (int j = 0; j < dag->n; j++)
            if (dag->adj[i * dag->n + j] || dag->adj[j * dag->n + i])
                undirected[i * dag->n + j] = undirected[j * dag->n + i] = 1;

    int count = 0;
    int *visited = calloc((size_t)dag->n, sizeof(int));
    int *parent = calloc((size_t)dag->n, sizeof(int));
    int *queue = calloc((size_t)dag->n, sizeof(int));
    if (!visited || !parent || !queue) {
        free(undirected); free(visited); free(parent); free(queue);
        return count;
    }

    for (int nb = 0; nb < dag->n && count < max_paths; nb++) {
        if (!dag->adj[nb * dag->n + from]) continue;
        if (nb == to) continue;

        memset(visited, 0, (size_t)dag->n * sizeof(int));
        for (int i = 0; i < dag->n; i++) parent[i] = -1;
        visited[from] = visited[nb] = 1;
        parent[nb] = from;

        int qh = 0, qt = 0;
        queue[qt++] = nb;

        while (qh < qt && count < max_paths) {
            int v = queue[qh++];
            if (v == to) {
                int rev[DAG_MAX_PATH_LEN];
                int rl = 0, cv = to;
                while (cv != -1 && rl < DAG_MAX_PATH_LEN) {
                    rev[rl++] = cv;
                    cv = parent[cv];
                }
                bool collider_free = true;
                for (int k = 1; k < rl - 1 && collider_free; k++) {
                    if (is_collider_triple(dag, rev[k + 1], rev[k], rev[k - 1]))
                        collider_free = false;
                }
                if (collider_free && rl >= 2) {
                    Path* pp = calloc(1, sizeof(Path));
                    if (pp) {
                        pp->nodes = calloc((size_t)rl, sizeof(int));
                        if (pp->nodes) {
                            for (int k = 0; k < rl; k++)
                                pp->nodes[k] = rev[rl - 1 - k];
                            pp->length = rl;
                            pp->is_directed = false;
                            paths[count++] = pp;
                        } else { free(pp); }
                    }
                }
                continue;
            }

            for (int w = 0; w < dag->n; w++) {
                if (undirected[v * dag->n + w] && !visited[w]) {
                    visited[w] = 1;
                    parent[w] = v;
                    queue[qt++] = w;
                }
            }
        }
    }

    free(undirected); free(visited); free(parent); free(queue);
    return count;
}

/* ================================================================
 * d-Separation Criterion
 *
 * X and Y are d-separated by Z iff every undirected path between
 * X and Y is blocked by Z:
 *   - Chain i->m->j or fork i<-m->j: blocked if m in Z
 *   - Collider i->m<-j: blocked if m NOT in Z and no descendant of m in Z
 * ================================================================ */

static bool node_or_descendant_in_set(const DAG* dag, int node,
                                       const int *z, int n_z) {
    for (int i = 0; i < n_z; i++)
        if (z[i] == node) return true;
    int desc[DAG_MAX_NODES];
    int nd = dag_descendants(dag, node, desc, DAG_MAX_NODES);
    for (int i = 0; i < n_z; i++)
        for (int j = 0; j < nd; j++)
            if (desc[j] == z[i]) return true;
    return false;
}

DSepResult* dag_d_separation(const DAG* dag, int x, int y,
                              const int *z, int n_z) {
    DSepResult* result = calloc(1, sizeof(DSepResult));
    if (!result || !dag || x < 0 || x >= dag->n || y < 0 || y >= dag->n) {
        if (result) result->is_separated = false;
        return result;
    }
    result->is_separated = true;

    int *ud = calloc((size_t)dag->n * dag->n, sizeof(int));
    if (!ud) { result->is_separated = false; return result; }
    for (int i = 0; i < dag->n; i++)
        for (int j = 0; j < dag->n; j++)
            if (dag->adj[i * dag->n + j] || dag->adj[j * dag->n + i])
                ud[i * dag->n + j] = 1;

    int *vis = calloc((size_t)dag->n, sizeof(int));
    int *par = calloc((size_t)dag->n, sizeof(int));
    int *q = calloc((size_t)dag->n, sizeof(int));
    if (!vis || !par || !q) { free(ud); free(vis); free(par); free(q);
        result->is_separated = false; return result; }

    for (int nb = 0; nb < dag->n; nb++) {
        if (!ud[x * dag->n + nb]) continue;
        memset(vis, 0, (size_t)dag->n * sizeof(int));
        for (int i = 0; i < dag->n; i++) par[i] = -1;
        vis[x] = vis[nb] = 1; par[nb] = x;
        int qh = 0, qt = 0; q[qt++] = nb;

        while (qh < qt) {
            int v = q[qh++];
            if (v == y) {
                int rev[DAG_MAX_PATH_LEN];
                int rl = 0, cv = y;
                while (cv != -1 && rl < DAG_MAX_PATH_LEN) {
                    rev[rl++] = cv; cv = par[cv];
                }
                bool blocked = false;
                for (int k = 1; k < rl - 1; k++) {
                    int a = rev[k + 1], b = rev[k], c = rev[k - 1];
                    bool a2b = dag->adj[a * dag->n + b] != 0;
                    bool b2a = dag->adj[b * dag->n + a] != 0;
                    bool c2b = dag->adj[c * dag->n + b] != 0;
                    bool b2c = dag->adj[b * dag->n + c] != 0;
                    bool b_in_z = false;
                    for (int zi = 0; zi < n_z; zi++)
                        if (z[zi] == b) { b_in_z = true; break; }

                    if ((a2b && b2c) || (b2a && b2c)) {
                        if (b_in_z) { blocked = true; break; }
                    } else if (a2b && c2b) {
                        if (!b_in_z && !node_or_descendant_in_set(dag, b, z, n_z))
                        { blocked = true; break; }
                    }
                }
                if (!blocked) {
                    result->is_separated = false;
                    Path* bp = calloc(1, sizeof(Path));
                    if (bp) {
                        bp->nodes = calloc((size_t)rl, sizeof(int));
                        if (bp->nodes) {
                            for (int k = 0; k < rl; k++)
                                bp->nodes[k] = rev[rl - 1 - k];
                            bp->length = rl;
                            result->blocking_path = bp;
                        } else free(bp);
                    }
                    free(ud); free(vis); free(par); free(q);
                    return result;
                }
                break;
            }
            for (int w = 0; w < dag->n; w++) {
                if (ud[v * dag->n + w] && !vis[w]) {
                    vis[w] = 1; par[w] = v; q[qt++] = w;
                }
            }
        }
    }

    free(ud); free(vis); free(par); free(q);
    return result;
}

/* ================================================================
 * Moral Graph: connect all parents (marry them), make undirected
 * ================================================================ */

DAG* dag_moral_graph(const DAG* dag) {
    if (!dag) return NULL;
    DAG* moral = dag_create(dag->n);
    if (!moral) return NULL;

    for (int i = 0; i < dag->n; i++)
        for (int j = 0; j < dag->n; j++)
            if (dag->adj[i * dag->n + j]) {
                moral->adj[i * moral->n + j] = 1;
                moral->adj[j * moral->n + i] = 1;
            }

    for (int v = 0; v < dag->n; v++) {
        int np = dag->n_parents[v];
        for (int i = 0; i < np; i++) {
            for (int j = i + 1; j < np; j++) {
                int p1 = dag->parents[v * dag->n + i];
                int p2 = dag->parents[v * dag->n + j];
                if (!moral->adj[p1 * moral->n + p2]) {
                    moral->adj[p1 * moral->n + p2] = 1;
                    moral->adj[p2 * moral->n + p1] = 1;
                    int c1 = moral->n_children[p1];
                    moral->children[p1 * moral->n + c1] = p2;
                    moral->n_children[p1]++;
                    int c2 = moral->n_children[p2];
                    moral->children[p2 * moral->n + c2] = p1;
                    moral->n_children[p2]++;
                }
            }
        }
    }
    return moral;
}

/* ================================================================
 * Graph Surgery ！ do(X=x): remove all edges into X
 * ================================================================ */

DAG* dag_do_intervention(const DAG* dag, int node) {
    if (!dag || node < 0 || node >= dag->n) return NULL;
    DAG* g_do = dag_create(dag->n);
    if (!g_do) return NULL;

    for (int i = 0; i < dag->n; i++)
        for (int j = 0; j < dag->n; j++)
            if (dag->adj[i * dag->n + j] && j != node)
                dag_add_edge(g_do, i, j);

    for (int i = 0; i < dag->n && i < g_do->n; i++)
        if (dag->names[i] && g_do->names[i])
            snprintf(g_do->names[i], DAG_NAME_LEN, "%s", dag->names[i]);

    return g_do;
}
