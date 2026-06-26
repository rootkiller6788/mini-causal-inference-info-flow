#ifndef SCM_COUNTERFACTUAL_H
#define SCM_COUNTERFACTUAL_H
#include "scm_core.h"
typedef struct { double value; double probability; bool identifiable; } SCM_Counterfactual;
typedef struct { double pn; double ps; double pns; } SCM_ProbNecessitySuff;
SCM_Counterfactual* scm_counterfactual_compute(const SCM_Model* m, int target, const int* evidence_vars, const double* evidence_vals, int n_evidence, const int* intervention_vars, const double* intervention_vals, int n_intervention);
void scm_counterfactual_free(SCM_Counterfactual* cf);
void scm_counterfactual_print(const SCM_Counterfactual* cf);
void scm_abduction_step(const SCM_Model* m, const int* evidence_vars, const double* evidence_vals, int n_evidence, double* u_posterior);
void scm_action_step(const SCM_Model* m, const int* intervention_vars, const double* intervention_vals, int n_intervention, SCM_Model* mutilated);
void scm_prediction_step(SCM_Model* mutilated, int target, double* result);
SCM_ProbNecessitySuff* scm_prob_necessity_sufficiency(const SCM_Model* m, int x, int y, double x_val, double y_val);
void scm_pns_free(SCM_ProbNecessitySuff* pns);
void scm_pns_print(const SCM_ProbNecessitySuff* pns);
double scm_effect_of_treatment_on_treated(const SCM_Model* m, int treatment, int outcome);
void scm_mediation_analysis(const SCM_Model* m, int treatment, int mediator, int outcome, double* nde, double* nie);
bool scm_monotonicity_check(const SCM_Model* m, int cause, int effect);
#endif
