#ifndef DCI_EFFECT_H
#define DCI_EFFECT_H

#include "dci_core.h"
#include "dci_graph.h"
#include "dci_backdoor.h"

/* ---- Causal Effect Estimation ---- */

/* Average Treatment Effect (ATE):
 * ATE = E[Y|do(X=1)] - E[Y|do(X=0)]
 */
typedef struct {
    double ate;
    double att;   /* Average Treatment Effect on the Treated */
    double atu;   /* Average Treatment Effect on the Untreated */
    double causal_risk_ratio;
    double causal_odds_ratio;
    bool is_identifiable;
} DCI_CausalEffect;

/* Potential outcome */
typedef struct {
    double y1;  /* Y under do(X=1) */
    double y0;  /* Y under do(X=0) */
} DCI_PotentialOutcome;

/* ATE computation */
DCI_CausalEffect dci_compute_ate(DCI_SCM* scm, int cause, int effect,
    int n_samples);
double dci_compute_att(DCI_SCM* scm, int cause, int effect,
    int n_samples);
double dci_compute_atu(DCI_SCM* scm, int cause, int effect,
    int n_samples);

/* Conditional Average Treatment Effect (CATE) */
double dci_compute_cate(DCI_SCM* scm, int cause, int effect,
    const double* condition_values, const int* condition_vars,
    int n_conditions, int n_samples);

/* Direct and indirect effects */
double dci_natural_direct_effect(DCI_SCM* scm, int cause, int effect,
    int mediator, double cause_val, double mediator_ref, int n_samples);
double dci_natural_indirect_effect(DCI_SCM* scm, int cause, int effect,
    int mediator, double cause_val, double mediator_ref, int n_samples);
double dci_controlled_direct_effect(DCI_SCM* scm, int cause, int effect,
    int mediator, double cause_val, double mediator_val, int n_samples);

/* Covariate adjustment methods */
double dci_standardization(DCI_SCM* scm, int cause, int effect,
    double cause_val, const int* covariates, int n_cov, int n_samples);
double dci_inverse_probability_weighting(DCI_SCM* scm, int cause,
    int effect, double cause_val, const int* covariates, int n_cov,
    int n_samples);
double dci_doubly_robust(DCI_SCM* scm, int cause, int effect,
    double cause_val, const int* covariates, int n_cov, int n_samples);
double dci_stratification(DCI_SCM* scm, int cause, int effect,
    double cause_val, const int* covariates, int n_cov, int n_samples);

/* Sensitivity analysis */
double dci_e_value(DCI_CausalEffect* effect, double null_hypothesis);

/* Print effect estimate */
void dci_effect_print(const DCI_CausalEffect* effect);

typedef struct {
    double lower; double upper; double estimate; double std_error;
} DCI_ConfidenceInterval;
DCI_ConfidenceInterval dci_ate_bootstrap(DCI_SCM* scm, int cause, int effect, int n_samples, int n_bootstrap, double alpha);
double dci_matching_estimator(DCI_SCM* scm, int cause, int effect, const int* covariates, int n_cov, int n_samples);
double dci_regression_adjust(DCI_SCM* scm, int cause, int effect, const int* covariates, int n_cov, int n_samples);
double dci_marginal_structural_model(DCI_SCM* scm, int treatment, int outcome, double treat_val, int n_samples);
double dci_difference_in_differences(DCI_SCM* scm, int cause, int effect, int time_var, int n_samples);
double dci_regression_discontinuity(DCI_SCM* scm, int cause, int effect, int running_var, double cutoff, int n_samples);

#endif
