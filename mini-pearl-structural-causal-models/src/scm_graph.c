#include "scm_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

bool scm_is_collider(const SCM_Model* m, int a, int b, int c) {
    return scm_has_edge(m, a, c) && scm_has_edge(m, b, c) && !scm_has_edge(m, a, b) && !scm_has_edge(m, b, a);
}
void scm_find_colliders(const SCM_Model* m, SCM_Collider* colliders, int* n) {
    *n = 0;
    for (int i = 0; i < m->n_vars; i++) {
        SCM_VarSet pa;
        scm_get_parents(m, i, &pa);
        for (int j = 0; j < pa.n; j++)
            for (int k = j+1; k < pa.n; k++)
                if (!scm_has_edge(m, pa.nodes[j], pa.nodes[k])
                    && !scm_has_edge(m, pa.nodes[k], pa.nodes[j])) {
                    colliders[*n].parent1 = pa.nodes[j];
                    colliders[*n].parent2 = pa.nodes[k];
                    colliders[*n].child = i;
                    (*n)++;
                }
    }
}
/* Depth-first reachability test in mutilated graph.
 * Not used in current d-separation (BFS moral graph), kept for API completeness.
 *
 * Theorem: In a causal graph, X and Y are d-connected given Z iff there
 * exists an active path. This DFS tests the existence of any directed or
 * undirected path avoiding Z - used as a subroutine for path-based d-sep. */
static bool dfs_reachable(const SCM_Model* m, int cur, int target, bool* visited,
                           const SCM_VarSet* Z, bool use_directed) {
    if (cur == target) return true;
    visited[cur] = true;
    for (int v = 0; v < m->n_vars; v++) {
        if (visited[v]) continue;
        bool connected = use_directed ? m->adjacency[cur][v]
                       : (m->adjacency[cur][v] || m->adjacency[v][cur]);
        if (!connected) continue;
        if (scm_varset_contains(Z, v)) continue;
        if (dfs_reachable(m, v, target, visited, Z, use_directed)) return true;
    }
    return false;
}

bool scm_active_path_exists(const SCM_Model* m, int from, int to,
                             const SCM_VarSet* Z) {
    if (!m) return false;
    bool visited[SCM_MAX_VARS] = {false};
    return dfs_reachable(m, from, to, visited, Z, false);
}
bool scm_d_separated(const SCM_Model* m, int x, int y, const SCM_VarSet* Z) {
    SCM_VarSet an;
    scm_varset_init(&an); scm_varset_add(&an, x); scm_varset_add(&an, y);
    for (int i = 0; i < Z->n; i++) {
        SCM_VarSet an_z; scm_get_ancestors(m, Z->nodes[i], &an_z);
        scm_varset_union(&an, &an_z, &an);
    }
    bool moral[SCM_MAX_VARS][SCM_MAX_VARS] = {{false}};
    for (int i = 0; i < m->n_edges; i++) {
        int u = m->edges[i].from, v = m->edges[i].to;
        if (scm_varset_contains(&an, u) && scm_varset_contains(&an, v))
            moral[u][v] = moral[v][u] = true;
    }
    SCM_Collider colliders[64]; int nc;
    scm_find_colliders(m, colliders, &nc);
    for (int i = 0; i < nc; i++)
        if (scm_varset_contains(&an, colliders[i].parent1)
            && scm_varset_contains(&an, colliders[i].parent2)
            && scm_varset_contains(&an, colliders[i].child)
            && scm_varset_contains(Z, colliders[i].child)) {
            moral[colliders[i].parent1][colliders[i].parent2] = true;
            moral[colliders[i].parent2][colliders[i].parent1] = true;
        }
    for (int i = 0; i < Z->n; i++) {
        for (int j = 0; j < m->n_vars; j++) {
            if (moral[Z->nodes[i]][j]) {
                moral[Z->nodes[i]][j] = false;
                moral[j][Z->nodes[i]] = false;
            }
        }
    }
    bool visited[SCM_MAX_VARS] = {false};
    int q[SCM_MAX_VARS], front=0, back=0;
    q[back++] = x; visited[x] = true;
    while (front < back) {
        int u = q[front++];
        if (u == y) return false;
        for (int v = 0; v < m->n_vars; v++) {
            if (!visited[v] && moral[u][v]) {
                visited[v] = true;
                q[back++] = v;
            }
        }
    }
    return true;
}
bool scm_d_separated_set(const SCM_Model* m, const SCM_VarSet* X,
                          const SCM_VarSet* Y, const SCM_VarSet* Z) {
    for (int i = 0; i < X->n; i++)
        for (int j = 0; j < Y->n; j++)
            if (!scm_d_separated(m, X->nodes[i], Y->nodes[j], Z)) return false;
    return true;
}
void scm_moral_graph(const SCM_Model* m, bool moral[SCM_MAX_VARS][SCM_MAX_VARS]) {
    memset(moral, 0, SCM_MAX_VARS*SCM_MAX_VARS*sizeof(bool));
    for (int i = 0; i < m->n_edges; i++)
        moral[m->edges[i].from][m->edges[i].to] = moral[m->edges[i].to][m->edges[i].from] = true;
    SCM_Collider colliders[64]; int nc;
    scm_find_colliders(m, colliders, &nc);
    for (int i = 0; i < nc; i++)
        moral[colliders[i].parent1][colliders[i].parent2] = moral[colliders[i].parent2][colliders[i].parent1] = true;
}
bool scm_find_directed_path(const SCM_Model* m, int from, int to, SCM_Path* p) {
    bool visited[SCM_MAX_VARS] = {false};
    int stack[SCM_MAX_VARS], parent[SCM_MAX_VARS];
    int top = 0;
    for (int i = 0; i < m->n_vars; i++) parent[i] = -1;
    stack[top++] = from; visited[from] = true;
    while (top > 0) {
        int u = stack[--top];
        if (u == to) { break; }
        for (int v = 0; v < m->n_vars; v++) {
            if (m->adjacency[u][v] && !visited[v]) {
                visited[v] = true; parent[v] = u; stack[top++] = v;
            }
        }
    }
    if (parent[to] == -1 && from != to) return false;
    p->len = 0;
    for (int cur = to; cur != -1; cur = parent[cur]) p->nodes[p->len++] = cur;
    for (int i = 0; i < p->len/2; i++) { int t=p->nodes[i]; p->nodes[i]=p->nodes[p->len-1-i]; p->nodes[p->len-1-i]=t; }
    return true;
}
int scm_all_directed_paths(const SCM_Model* m, int from, int to, SCM_Path* paths, int max_paths) {
    int count = 0;
    int stack[SCM_MAX_VARS*SCM_MAX_PATHS];
    int path[SCM_MAX_VARS];
    bool on_path[SCM_MAX_VARS] = {false};
    int sp = 0;
    stack[sp++] = from; stack[sp++] = 0;
    on_path[from] = true; path[0] = from;
    int path_len = 1;
    while (sp > 0 && count < max_paths) {
        int idx = stack[--sp];
        int u = stack[--sp];
        if (u == to) {
            paths[count].len = path_len;
            for (int i = 0; i < path_len; i++) paths[count].nodes[i] = path[i];
            count++;
            on_path[u] = false; path_len--;
            continue;
        }
        int child_count = 0, next_child = -1;
        for (int v = 0; v < m->n_vars; v++) {
            if (m->adjacency[u][v] && !on_path[v]) {
                if (child_count == idx) { next_child = v; break; }
                child_count++;
            }
        }
        if (next_child >= 0) {
            stack[sp++] = u; stack[sp++] = idx+1;
            stack[sp++] = next_child; stack[sp++] = 0;
            on_path[next_child] = true; path[path_len++] = next_child;
        } else {
            on_path[u] = false; path_len--;
        }
    }
    return count;
}
void scm_path_print(const SCM_Path* p) {
    printf("Path[%d]: ", p->len);
    for (int i = 0; i < p->len; i++) printf("%d ", p->nodes[i]);
    printf("\n");
}
bool scm_is_blocked(const SCM_Path* p, const SCM_VarSet* Z, const SCM_Model* m) {
    for (int i = 1; i < p->len-1; i++) {
        int a = p->nodes[i-1], b = p->nodes[i], c = p->nodes[i+1];
        if (m->adjacency[a][b] && m->adjacency[c][b] && !scm_varset_contains(Z, b)) return true;
        if ((m->adjacency[a][b] || m->adjacency[b][a]) && scm_varset_contains(Z, b)) return true;
    }
    return false;
}
bool scm_path_has_collider(const SCM_Path* p, const SCM_Model* m) {
    for (int i = 1; i < p->len-1; i++)
        if (scm_is_collider(m, p->nodes[i-1], p->nodes[i+1], p->nodes[i])) return true;
    return false;
}
bool scm_is_ancestor(const SCM_Model* m, int a, int b) {
    SCM_VarSet an; scm_get_ancestors(m, b, &an);
    return scm_varset_contains(&an, a);
}
