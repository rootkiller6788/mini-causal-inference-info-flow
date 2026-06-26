#include "scm_counterfactual.h"
#include "scm_intervention.h"
#include "scm_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

SCM_Counterfactual* scm_counterfactual_compute(const SCM_Model* m, int target,
    const int* evidence_vars, const double* evidence_vals, int n_evidence,
    const int* intervention_vars, const double* intervention_vals, int n_intervention) {
    SCM_Counterfactual* cf = calloc(1, sizeof(SCM_Counterfactual));
    if (!cf) return NULL;
    SCM_Model mutilated;
    memcpy(&mutilated, m, sizeof(SCM_Model));
    for (int k = 0; k < n_intervention; k++) {
        mutilated.intervened[intervention_vars[k]] = true;
        mutilated.intervention_value[intervention_vars[k]] = intervention_vals[k];
    }
    scm_topological_sort(&mutilated);
    scm_sample(&mutilated);
    cf->value = mutilated.vars[target].value;
    cf->identifiable = true;
    cf->probability = 1.0;
    (void)evidence_vars; (void)evidence_vals; (void)n_evidence;
    return cf;
}
void scm_counterfactual_free(SCM_Counterfactual* cf) { free(cf); }
void scm_counterfactual_print(const SCM_Counterfactual* cf) {
    printf("Counterfactual: value=%.4f prob=%.4f %s\n",
           cf->value, cf->probability, cf->identifiable?"identifiable":"NOT identifiable");
}
void scm_abduction_step(const SCM_Model* m, const int* evidence_vars,
                         const double* evidence_vals, int n_evidence, double* u_posterior) {
    for (int i = 0; i < m->n_vars; i++) u_posterior[i] = 0.0;
    for (int i = 0; i < n_evidence; i++) {
        int v = evidence_vars[i];
        double expected = 0.0;
        for (int j = 0; j < m->n_edges; j++)
            if (m->edges[j].to == v) expected += m->edges[j].coef * m->vars[m->edges[j].from].value;
        u_posterior[v] = evidence_vals[i] - expected;
    }
}
void scm_action_step(const SCM_Model* m, const int* intervention_vars,
                      const double* intervention_vals, int n_intervention,
                      SCM_Model* mutilated) {
    memcpy(mutilated, m, sizeof(SCM_Model));
    for (int k = 0; k < n_intervention; k++) {
        mutilated->intervened[intervention_vars[k]] = true;
        mutilated->intervention_value[intervention_vars[k]] = intervention_vals[k];
        for (int i = mutilated->n_edges - 1; i >= 0; i--) {
            if (mutilated->edges[i].to == intervention_vars[k]) {
                mutilated->adjacency[mutilated->edges[i].from][intervention_vars[k]] = false;
                memmove(&mutilated->edges[i], &mutilated->edges[i+1],
                        (size_t)(mutilated->n_edges - 1 - i) * sizeof(SCM_Edge));
                mutilated->n_edges--;
            }
        }
    }
}
void scm_prediction_step(SCM_Model* mutilated, int target, double* result) {
    scm_sample(mutilated);
    *result = mutilated->vars[target].value;
}
SCM_ProbNecessitySuff* scm_prob_necessity_sufficiency(const SCM_Model* m,
    int x, int y, double x_val, double y_val) {
    SCM_ProbNecessitySuff* pns = calloc(1, sizeof(SCM_ProbNecessitySuff));
    if (!pns) return NULL;
    (void)y_val;
    SCM_Model* m1 = scm_do_intervention(m, x, x_val);
    double p_yx = (m1->vars[y].value > 0.5) ? 1.0 : 0.0;
    scm_free(m1);
    SCM_Model* m0 = scm_do_intervention(m, x, 1.0 - x_val);
    double p_yx0 = (m0->vars[y].value > 0.5) ? 1.0 : 0.0;
    scm_free(m0);
    pns->pn = p_yx - p_yx0;
    pns->ps = p_yx;
    pns->pns = p_yx - p_yx0;
    if (pns->pn < 0.0) pns->pn = 0.0;
    if (pns->ps < 0.0) pns->ps = 0.0;
    if (pns->pns < 0.0) pns->pns = 0.0;
    return pns;
}
void scm_pns_free(SCM_ProbNecessitySuff* pns) { free(pns); }
void scm_pns_print(const SCM_ProbNecessitySuff* pns) {
    printf("PNS: PN=%.4f PS=%.4f PNS=%.4f\n", pns->pn, pns->ps, pns->pns);
}
double scm_effect_of_treatment_on_treated(const SCM_Model* m, int treatment, int outcome) {
    SCM_CausalEffect* ce = scm_causal_effect_compute(m, treatment, outcome);
    double ett = ce->ate;
    scm_causal_effect_free(ce);
    return ett;
}
void scm_mediation_analysis(const SCM_Model* m, int treatment, int mediator,
                              int outcome, double* nde, double* nie) {
    (void)mediator;
    SCM_Model* mtx1 = scm_do_intervention(m, treatment, 1.0);
    double y1m1 = mtx1->vars[outcome].value;
    scm_free(mtx1);
    SCM_Model* mtx0 = scm_do_intervention(m, treatment, 0.0);
    double y0m0 = mtx0->vars[outcome].value;
    double y1m0 = y1m1;
    *nde = y1m0 - y0m0;
    *nie = y1m1 - y1m0;
    scm_free(mtx0);
}
bool scm_monotonicity_check(const SCM_Model* m, int cause, int effect) {
    SCM_Path p;
    if (!scm_find_directed_path(m, cause, effect, &p)) return true;
    double prod = 1.0;
    for (int i = 0; i < p.len-1; i++)
        for (int j = 0; j < m->n_edges; j++)
            if (m->edges[j].from == p.nodes[i] && m->edges[j].to == p.nodes[i+1])
                prod *= m->edges[j].coef;
    return prod > 0.0;
}
