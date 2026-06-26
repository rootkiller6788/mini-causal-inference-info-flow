#ifndef CFR_COUNTERFACTUAL_H
#define CFR_COUNTERFACTUAL_H

#include "cfr_core.h"

/* ============================================================
 * cfr_counterfactual.h -- Counterfactual Computation
 *
 * Implements Pearl's 3-step counterfactual algorithm:
 *   1. Abduction:  Update P(U) given evidence E=e
 *   2. Action:     Modify SCM by do(X=x)
 *   3. Prediction: Compute Y_x in modified model
 *
 * Also covers counterfactual distributions and
 * identifiability conditions.
 * ============================================================ */

typedef struct {
    CFRSCM* original_scm;     /* original model */
    CFRSCM* modified_scm;     /* intervened model */
    double* exogenous_posterior; /* P(U | evidence) */
    int n_exogenous;
    double* counterfactual_values; /* Y_x for each sample */
    int n_samples;
    double mean;
    double variance;
    bool identifiable;
} CFRCounterfactual;

typedef struct {
    char* name;
    int var_id;
    double intervention_value;
    double observed_value;
    double counterfactual_value;
    double difference;        /* CF - observed */
} CFRCFComparison;

/* --- Lifecycle --- */
CFRCounterfactual* cfr_cf_create(CFRSCM* scm);
void               cfr_cf_free(CFRCounterfactual* cf);

/* --- 3-Step Algorithm --- */
int cfr_cf_abduction(CFRCounterfactual* cf, CFREvidence* evidence);
int cfr_cf_action(CFRCounterfactual* cf, CFRIntervention* intervention);
int cfr_cf_prediction(CFRCounterfactual* cf, int target_id);
int cfr_cf_compute_full(CFRCounterfactual* cf,
                         CFREvidence* evidence,
                         CFRIntervention* intervention,
                         int target_id);

/* --- Monte Carlo Counterfactuals --- */
int cfr_cf_monte_carlo(CFRSCM* scm, CFREvidence* evidence,
                        CFRIntervention* intervention,
                        int target_id, int n_samples,
                        double* cf_distribution);

/* --- Counterfactual Distributions --- */
double cfr_cf_expected_value(CFRCounterfactual* cf);
double cfr_cf_variance(CFRCounterfactual* cf);
double cfr_cf_probability(CFRCounterfactual* cf,
                           double threshold, bool greater_than);

/* --- Identifiability --- */
bool cfr_cf_is_identifiable(CFRSCM* scm, int target_id,
                             CFRIntervention* intv, CFREvidence* ev);

/* --- Comparison --- */
CFRCFComparison* cfr_cf_compare(CFRSCM* scm,
                                 int var_id, double observed,
                                 double intervention_value,
                                 CFREvidence* evidence);
void cfr_cf_compare_free(CFRCFComparison* cmp);
int  cfr_cf_compare_interventions(CFRSCM* scm, int target_id,
                                   CFRIntervention* intv_a, CFRIntervention* intv_b,
                                   CFREvidence* ev, double* diff);

/* --- Conditional Counterfactual --- */
int cfr_cf_conditional_counterfactual(CFRSCM* scm,
    CFREvidence* ev, CFRIntervention* intv,
    int target_id, double condition_var_id, double condition_value,
    double* cf_value);

/* --- Probability of Causation --- */
double cfr_cf_probability_of_causation(CFRSCM* scm,
    CFREvidence* ev, CFRIntervention* intv_treat,
    CFRIntervention* intv_control, int target_id);

/* --- Distribution --- */
int cfr_cf_counterfactual_distribution(CFRSCM* scm,
    CFREvidence* ev, CFRIntervention* intv,
    int target_id, int n_samples, double* distribution);

/* --- Output --- */
void cfr_cf_print(CFRCounterfactual* cf);
void cfr_cf_print_comparison(CFRCFComparison* cmp);

/* --- Path-Specific Effects --- */
int cfr_cf_path_specific_effect(CFRSCM* scm, int* path_vars, int path_len,
                                 CFRIntervention* intv, CFREvidence* ev, double* pse);

/* --- Effect of Treatment on Treated --- */
double cfr_cf_ett(CFRSCM* scm, int treatment_id, int outcome_id, CFREvidence* ev);

/* --- Twin Network --- */
CFRSCM* cfr_cf_twin_network(CFRSCM* scm, int target_id);

/* --- Counterfactual Fairness --- */
double cfr_cf_fairness_measure(CFRSCM* scm, int protected_attr_id,
                                int decision_id, int outcome_id);

/* --- Transportability --- */
double cfr_cf_transportability(CFRSCM* scm_source, CFRSCM* scm_target,
                                int treatment_id, int outcome_id);

/* --- Sequential Counterfactuals --- */
int cfr_cf_sequential(CFRSCM* scm, CFRIntervention** intv_sequence, int n_intv,
                       CFREvidence* ev, int target_id, double* seq_cf);

/* --- Controlled Direct Effect via Counterfactuals --- */
double cfr_cf_controlled_direct_effect(CFRSCM* scm, int t_id, int m_id,
                                        int y_id, double m_value);

/* --- Mediation via Counterfactuals --- */
int cfr_cf_mediation_effects(CFRSCM* scm, int t_id, int m_id, int y_id,
                              double* nde, double* nie, double* te);

#endif /* CFR_COUNTERFACTUAL_H */
