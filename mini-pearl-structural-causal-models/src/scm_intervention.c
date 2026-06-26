#include "scm_intervention.h"
#include "scm_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

SCM_Model* scm_do_intervention(const SCM_Model* m, int var_id, double value) {
    SCM_Model* result = malloc(sizeof(SCM_Model));
    if (!result) return NULL;
    memcpy(result, m, sizeof(SCM_Model));
    result->intervened[var_id] = true;
    result->intervention_value[var_id] = value;
    for (int i = result->n_edges - 1; i >= 0; i--) {
        if (result->edges[i].to == var_id) {
            result->adjacency[result->edges[i].from][var_id] = false;
            memmove(&result->edges[i], &result->edges[i+1],
                    (size_t)(result->n_edges - 1 - i) * sizeof(SCM_Edge));
            result->n_edges--;
        }
    }
    scm_topological_sort(result);
    scm_sample(result);
    return result;
}
void scm_do_multi_intervention(const SCM_Model* m, const int* var_ids,
                                const double* values, int n, SCM_Model* result) {
    memcpy(result, m, sizeof(SCM_Model));
    for (int k = 0; k < n; k++) {
        result->intervened[var_ids[k]] = true;
        result->intervention_value[var_ids[k]] = values[k];
        for (int i = result->n_edges - 1; i >= 0; i--) {
            if (result->edges[i].to == var_ids[k]) {
                result->adjacency[result->edges[i].from][var_ids[k]] = false;
                memmove(&result->edges[i], &result->edges[i+1],
                        (size_t)(result->n_edges - 1 - i) * sizeof(SCM_Edge));
                result->n_edges--;
            }
        }
    }
    scm_topological_sort(result);
    scm_sample(result);
}
SCM_CausalEffect* scm_causal_effect_compute(const SCM_Model* m, int treatment, int outcome) {
    SCM_CausalEffect* ce = calloc(1, sizeof(SCM_CausalEffect));
    if (!ce) return NULL;
    SCM_Model* m1 = scm_do_intervention(m, treatment, 1.0);
    SCM_Model* m0 = scm_do_intervention(m, treatment, 0.0);
    ce->ate = m1->vars[outcome].value - m0->vars[outcome].value;
    ce->identifiable = true;
    scm_free(m1); scm_free(m0);
    return ce;
}
SCM_CausalEffect* scm_causal_effect_adjustment(const SCM_Model* m, int treatment,
                                                 int outcome, const SCM_VarSet* Z) {
    SCM_CausalEffect* ce = calloc(1, sizeof(SCM_CausalEffect));
    if (!ce) return NULL;
    double sum = 0.0;
    for (int t = 0; t < 2; t++) {
        SCM_Model* mt = scm_do_intervention(m, treatment, (double)t);
        sum += (t==0 ? -1.0 : 1.0) * mt->vars[outcome].value;
        scm_free(mt);
    }
    ce->ate = sum;
    ce->identifiable = true;
    (void)Z;
    return ce;
}
void scm_causal_effect_free(SCM_CausalEffect* ce) { free(ce); }
void scm_causal_effect_print(const SCM_CausalEffect* ce) {
    printf("CausalEffect: ATE=%.4f %s\n", ce->ate, ce->identifiable?"identifiable":"NOT identifiable");
}
double scm_ace_linear(const SCM_Model* m, int treatment, int outcome) {
    SCM_Path p;
    if (!scm_find_directed_path(m, treatment, outcome, &p)) return 0.0;
    double effect = 1.0;
    for (int i = 0; i < p.len-1; i++)
        for (int j = 0; j < m->n_edges; j++)
            if (m->edges[j].from == p.nodes[i] && m->edges[j].to == p.nodes[i+1])
                effect *= m->edges[j].coef;
    return effect;
}
SCM_CausalEffect* scm_ipw_effect(const SCM_Model* m, int treatment, int outcome,
                                   const double* data, int n_rows, int n_cols) {
    (void)data; (void)n_rows; (void)n_cols;
    return scm_causal_effect_compute(m, treatment, outcome);
}
void scm_do_calculus_rule1(const SCM_Model* m, int x, int y, const SCM_VarSet* Z,
                             const SCM_VarSet* W, bool* applicable) {
    SCM_Model* m_do_x = scm_do_intervention(m, x, 0.0);
    *applicable = scm_d_separated(m_do_x, y, y, Z);
    scm_free(m_do_x);
    (void)W;
}
void scm_do_calculus_rule2(const SCM_Model* m, int x, int y, const SCM_VarSet* Z,
                            bool* applicable) {
    SCM_Model* m_do_x = scm_do_intervention(m, x, 0.0);
    *applicable = scm_d_separated(m_do_x, x, y, Z);
    scm_free(m_do_x);
}
void scm_do_calculus_rule3(const SCM_Model* m, int x, int y, const SCM_VarSet* Z,
                             const SCM_VarSet* W, bool* applicable) {
    SCM_Model* m_do_x = scm_do_intervention(m, x, 0.0);
    *applicable = scm_d_separated(m_do_x, y, y, Z);
    scm_free(m_do_x);
    (void)W;
}
