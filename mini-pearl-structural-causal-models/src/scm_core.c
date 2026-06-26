#include "scm_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

SCM_Model* scm_create(void) {
    SCM_Model* m = calloc(1, sizeof(SCM_Model));
    if (!m) return NULL;
    m->n_vars = 0; m->n_edges = 0;
    return m;
}
void scm_free(SCM_Model* m) { free(m); }

int scm_add_variable(SCM_Model* m, const char* name, SCM_VarType type) {
    if (m->n_vars >= SCM_MAX_VARS) return -1;
    int id = m->n_vars;
    SCM_Variable* v = &m->vars[id];
    v->id = id; v->type = type; v->value = 0.0; v->noise = 0.0;
    strncpy(v->name, name, SCM_NAME_LEN-1);
    v->name[SCM_NAME_LEN-1] = '\0';
    m->n_vars++;
    return id;
}
void scm_set_equation(SCM_Model* m, int var_id, SCM_Equation eq) {
    if (var_id >= 0 && var_id < m->n_vars) m->equations[var_id] = eq;
}
bool scm_is_exogenous(const SCM_Model* m, int var_id) {
    return var_id >= 0 && var_id < m->n_vars && m->vars[var_id].type == SCM_EXOGENOUS;
}
int scm_add_edge(SCM_Model* m, int from, int to, double coef) {
    if (m->n_edges >= SCM_MAX_EDGES) return -1;
    SCM_Edge* e = &m->edges[m->n_edges];
    e->from = from; e->to = to; e->coef = coef;
    m->adjacency[from][to] = true;
    m->n_edges++;
    return m->n_edges - 1;
}
bool scm_has_edge(const SCM_Model* m, int from, int to) {
    return m->adjacency[from][to];
}
void scm_build_adjacency(SCM_Model* m) {
    memset(m->adjacency, 0, sizeof(m->adjacency));
    for (int i = 0; i < m->n_edges; i++)
        m->adjacency[m->edges[i].from][m->edges[i].to] = true;
}
void scm_topological_sort(SCM_Model* m) {
    int indeg[SCM_MAX_VARS] = {0};
    for (int i = 0; i < m->n_edges; i++) indeg[m->edges[i].to]++;
    int q[SCM_MAX_VARS], front=0, back=0;
    for (int i = 0; i < m->n_vars; i++)
        if (indeg[i] == 0) q[back++] = i;
    int idx = 0;
    while (front < back) {
        int u = q[front++];
        m->topo_order[idx++] = u;
        for (int v = 0; v < m->n_vars; v++) {
            if (m->adjacency[u][v]) {
                indeg[v]--;
                if (indeg[v] == 0) q[back++] = v;
            }
        }
    }
}
bool scm_is_dag(const SCM_Model* m) {
    int indeg[SCM_MAX_VARS] = {0};
    for (int i = 0; i < m->n_edges; i++) indeg[m->edges[i].to]++;
    int q[SCM_MAX_VARS], front=0, back=0;
    for (int i = 0; i < m->n_vars; i++)
        if (indeg[i] == 0) q[back++] = i;
    int cnt = 0;
    while (front < back) {
        int u = q[front++]; cnt++;
        for (int v = 0; v < m->n_vars; v++)
            if (m->adjacency[u][v] && --indeg[v] == 0) q[back++] = v;
    }
    return cnt == m->n_vars;
}
void scm_get_parents(const SCM_Model* m, int var_id, SCM_VarSet* pa) {
    scm_varset_init(pa);
    for (int i = 0; i < m->n_edges; i++)
        if (m->edges[i].to == var_id) scm_varset_add(pa, m->edges[i].from);
}
void scm_get_children(const SCM_Model* m, int var_id, SCM_VarSet* ch) {
    scm_varset_init(ch);
    for (int i = 0; i < m->n_edges; i++)
        if (m->edges[i].from == var_id) scm_varset_add(ch, m->edges[i].to);
}
void scm_get_ancestors(const SCM_Model* m, int var_id, SCM_VarSet* an) {
    scm_varset_init(an);
    bool visited[SCM_MAX_VARS] = {false};
    int stack[SCM_MAX_VARS], top=0;
    stack[top++] = var_id;
    while (top > 0) {
        int u = stack[--top];
        SCM_VarSet pa; scm_get_parents(m, u, &pa);
        for (int i = 0; i < pa.n; i++) {
            int p = pa.nodes[i];
            if (!visited[p]) { visited[p]=true; scm_varset_add(an,p); stack[top++]=p; }
        }
    }
}
void scm_get_descendants(const SCM_Model* m, int var_id, SCM_VarSet* de) {
    scm_varset_init(de);
    bool visited[SCM_MAX_VARS] = {false};
    int stack[SCM_MAX_VARS], top=0;
    stack[top++] = var_id;
    while (top > 0) {
        int u = stack[--top];
        for (int v = 0; v < m->n_vars; v++) {
            if (m->adjacency[u][v] && !visited[v]) {
                visited[v]=true; scm_varset_add(de,v); stack[top++]=v;
            }
        }
    }
}
void scm_markov_blanket(const SCM_Model* m, int var_id, SCM_VarSet* mb) {
    scm_varset_init(mb);
    SCM_VarSet pa, ch;
    scm_get_parents(m, var_id, &pa);
    scm_get_children(m, var_id, &ch);
    scm_varset_union(&pa, &ch, mb);
    for (int i = 0; i < ch.n; i++) {
        SCM_VarSet pa_ch; scm_get_parents(m, ch.nodes[i], &pa_ch);
        scm_varset_union(mb, &pa_ch, mb);
    }
}
void scm_sample(SCM_Model* m) {
    scm_topological_sort(m);
    for (int i = 0; i < m->n_vars; i++) {
        int id = m->topo_order[i];
        if (m->intervened[id]) {
            m->vars[id].value = m->intervention_value[id];
        } else if (m->equations[id]) {
            m->vars[id].value = m->equations[id](m, id);
        } else {
            double sum = 0.0;
            for (int j = 0; j < m->n_edges; j++)
                if (m->edges[j].to == id) sum += m->edges[j].coef * m->vars[m->edges[j].from].value;
            m->vars[id].value = sum + m->vars[id].noise;
        }
    }
}
void scm_sample_n(SCM_Model* m, double* data, int n_rows) {
    for (int r = 0; r < n_rows; r++) {
        for (int i = 0; i < m->n_vars; i++) m->vars[i].noise = ((double)rand()/RAND_MAX - 0.5)*2.0;
        scm_sample(m);
        for (int c = 0; c < m->n_vars; c++) data[r*m->n_vars + c] = m->vars[c].value;
    }
}
double scm_evaluate(SCM_Model* m, int var_id) {
    if (var_id < 0 || var_id >= m->n_vars) return 0.0;
    if (m->intervened[var_id]) return m->intervention_value[var_id];
    if (m->equations[var_id]) return m->equations[var_id](m, var_id);
    double sum = 0.0;
    for (int j = 0; j < m->n_edges; j++)
        if (m->edges[j].to == var_id) sum += m->edges[j].coef * m->vars[m->edges[j].from].value;
    return sum + m->vars[var_id].noise;
}
void scm_print(const SCM_Model* m) {
    printf("SCM: %d vars, %d edges\n", m->n_vars, m->n_edges);
    for (int i = 0; i < m->n_vars; i++)
        printf("  V[%d]=%s type=%d val=%.3f\n", i, m->vars[i].name, m->vars[i].type, m->vars[i].value);
}
double scm_linear_eq(SCM_Model* m, int var_id) {
    double sum = 0.0;
    for (int j = 0; j < m->n_edges; j++)
        if (m->edges[j].to == var_id) sum += m->edges[j].coef * m->vars[m->edges[j].from].value;
    return sum + m->vars[var_id].noise;
}
void scm_varset_init(SCM_VarSet* vs) { vs->n = 0; }
void scm_varset_add(SCM_VarSet* vs, int id) {
    if (vs->n < SCM_MAX_VARS) vs->nodes[vs->n++] = id;
}
bool scm_varset_contains(const SCM_VarSet* vs, int id) {
    for (int i = 0; i < vs->n; i++) if (vs->nodes[i]==id) return true;
    return false;
}
void scm_varset_union(SCM_VarSet* a, const SCM_VarSet* b, SCM_VarSet* out) {
    scm_varset_init(out);
    for (int i = 0; i < a->n; i++) if (!scm_varset_contains(out, a->nodes[i])) scm_varset_add(out, a->nodes[i]);
    for (int i = 0; i < b->n; i++) if (!scm_varset_contains(out, b->nodes[i])) scm_varset_add(out, b->nodes[i]);
}
void scm_varset_intersect(SCM_VarSet* a, const SCM_VarSet* b, SCM_VarSet* out) {
    scm_varset_init(out);
    for (int i = 0; i < a->n; i++) if (scm_varset_contains(b, a->nodes[i])) scm_varset_add(out, a->nodes[i]);
}
int scm_varset_size(const SCM_VarSet* vs) { return vs->n; }
