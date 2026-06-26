#ifndef SCM_INTERVENTION_H
#define SCM_INTERVENTION_H
#include "scm_core.h"
typedef struct { double ate; double att; double atc; bool identifiable; } SCM_CausalEffect;
SCM_Model* scm_do_intervention(const SCM_Model* m, int var_id, double value);
void scm_do_multi_intervention(const SCM_Model* m, const int* var_ids, const double* values, int n, SCM_Model* result);
double scm_truncated_factorization(const SCM_Model* m, const double* values);
SCM_CausalEffect* scm_causal_effect_compute(const SCM_Model* m, int treatment, int outcome);
SCM_CausalEffect* scm_causal_effect_adjustment(const SCM_Model* m, int treatment, int outcome, const SCM_VarSet* adjustment);
SCM_CausalEffect* scm_ipw_effect(const SCM_Model* m, int treatment, int outcome, const double* data, int n_rows, int n_cols);
void scm_causal_effect_free(SCM_CausalEffect* ce);
void scm_causal_effect_print(const SCM_CausalEffect* ce);
double scm_ace_linear(const SCM_Model* m, int treatment, int outcome);
void scm_do_calculus_rule1(const SCM_Model* m, int x, int y, const SCM_VarSet* Z, const SCM_VarSet* W, bool* applicable);
void scm_do_calculus_rule2(const SCM_Model* m, int x, int y, const SCM_VarSet* Z, bool* applicable);
void scm_do_calculus_rule3(const SCM_Model* m, int x, int y, const SCM_VarSet* Z, const SCM_VarSet* W, bool* applicable);
#endif
