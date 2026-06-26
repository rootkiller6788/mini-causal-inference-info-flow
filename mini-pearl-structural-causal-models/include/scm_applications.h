#ifndef SCM_APPLICATIONS_H
#define SCM_APPLICATIONS_H
#include "scm_core.h"

/* ============================================================================
 * Real-world SCM Applications (L7)
 * References:
 *   Hernan & Robins (2020) "Causal Inference: What If" (CRC Press)
 *   Morgan & Winship (2015) "Counterfactuals and Causal Inference" (Cambridge)
 *   Angrist & Pischke (2009) "Mostly Harmless Econometrics" (Princeton)
 *   Imbens & Rubin (2015) "Causal Inference for Statistics, Social, and Biomedical Sciences"
 * ============================================================================ */

/* Simpson's Paradox detection result */
typedef struct {
    bool paradox_detected;      /* True if reversal found */
    double aggregate_effect;    /* Effect ignoring subgroups */
    double subgroup_effects[16];/* Effects within subgroups */
    int    n_subgroups;         /* Number of subgroups analyzed */
    double simpson_correction;  /* Causal effect after deconfounding */
} SCM_SimpsonResult;

/* Clinical trial analysis result */
typedef struct {
    double itt_effect;          /* Intention-to-treat */
    double per_protocol;        /* Per-protocol effect */
    double as_treated;          /* As-treated effect */
    double complier_ace;        /* Complier average causal effect (LATE) */
    double iv_wald_estimate;    /* Instrumental variables Wald estimator */
    bool   noncompliance_present;
} SCM_ClinicalResult;

/* Propensity score matching result */
typedef struct {
    double att;                 /* Average treatment effect on treated */
    double atc;                 /* Average treatment effect on controls */
    double ate_matched;         /* ATE after matching */
    int    n_matched_pairs;     /* Number of matched pairs */
    double balance_std_diff;    /* Max standardized difference after matching */
    bool   balance_achieved;    /* Whether balance threshold met */
} SCM_PSMResult;

/* === Simpson's Paradox (L7: Epidemiology) === */

/* Detect Simpson's paradox in data: aggregate effect reverses in subgroups.
 * treatment: column index of treatment variable
 * outcome: column index of outcome
 * subgroup: column index of categorical subgroup variable
 * data: n_rows x n_cols matrix */
SCM_SimpsonResult* scm_simpsons_paradox_detect(const double* data, int n_rows, int n_cols,
                                                int treatment, int outcome, int subgroup);

/* Resolve the paradox: compute causal effect adjusting for subgroup (back-door).
 * Returns the adjusted causal effect. */
double scm_simpsons_paradox_resolve(const SCM_Model* m, int treatment, int outcome,
                                     const SCM_VarSet* adjustment_set);

/* === Causal Paradoxes (L6: Canonical Problems) === */

/* Construct an M-bias graph and test if conditioning on M creates bias.
 * M-bias: M is a collider on a path, conditioning opens the path.
 * Returns true if conditioning on M creates spurious association. */
bool scm_m_bias_test(const SCM_Model* m, int x, int y, int m_collider);

/* Test Berkson's paradox: conditioning on a collider creates bias.
 * In hospital data, conditioning on "hospitalized" creates spurious
 * negative correlation between diseases. */
bool scm_berkson_paradox_test(const SCM_Model* m, int disease_a, int disease_b,
                               int hospitalized);

/* === Clinical Trial Analysis (L7: Medicine/Epidemiology) === */

/* Analyze a clinical trial with possible noncompliance.
 * treatment_assigned: random assignment (instrument)
 * treatment_received: actual treatment taken
 * outcome: outcome variable
 * data: observational or trial data */
SCM_ClinicalResult* scm_clinical_trial_analysis(const double* data, int n_rows, int n_cols,
                                                  int treatment_assigned, int treatment_received,
                                                  int outcome);

/* Compute Complier Average Causal Effect (CACE / LATE).
 * Under monotonicity and exclusion restriction assumptions. */
double scm_complier_ace(const double* data, int n_rows, int n_cols,
                         int instrument, int treatment, int outcome);

/* === Instrumental Variables (L7: Econometrics) === */

/* Two-Stage Least Squares (2SLS) estimation.
 * First stage: regress treatment on instrument + covariates.
 * Second stage: regress outcome on fitted treatment + covariates.
 * Returns the IV causal effect estimate. */
double scm_iv_two_stage_least_squares(const double* data, int n_rows, int n_cols,
                                       int instrument, int treatment, int outcome,
                                       const int* covariates, int n_covariates);

/* Wald estimator for binary IV and binary treatment:
 * beta_IV = (E[Y|Z=1] - E[Y|Z=0]) / (E[X|Z=1] - E[X|Z=0]) */
double scm_iv_wald_estimator(const double* data, int n_rows, int n_cols,
                              int instrument, int treatment, int outcome);

/* === Propensity Score Methods (L7) === */

/* Estimate propensity scores using logistic regression (linear approximation).
 * Returns array of propensity scores (length n_rows). */
double* scm_propensity_score_estimate(const double* data, int n_rows, int n_cols,
                                       int treatment, const int* covariates, int n_covariates);

/* Nearest-neighbor propensity score matching.
 * Matches each treated unit to the closest control unit.
 * Returns matched dataset info in the result struct. */
SCM_PSMResult* scm_propensity_score_matching(const double* data, int n_rows, int n_cols,
                                              int treatment, int outcome,
                                              const int* covariates, int n_covariates,
                                              double caliper);

/* === Difference-in-Differences (L7: Policy Evaluation) === */

/* Simple 2x2 DiD estimator.
 * pre_treatment: column of pre-period outcome
 * post_outcome: column of post-period outcome
 * treated_group: column of treatment group indicator (0/1)
 * Returns the DiD estimate (average treatment effect on treated). */
double scm_difference_in_differences(const double* data, int n_rows,
                                      int pre_outcome, int post_outcome,
                                      int treated_group);

/* === Regression Discontinuity Design (L7) === */

/* Sharp RDD: treatment = 1 if running_var >= cutoff.
 * Estimates treatment effect using local linear regression.
 * bandwidth: window around cutoff. */
double scm_regression_discontinuity(const double* data, int n_rows, int n_cols,
                                     int running_var, int outcome, int treatment,
                                     double cutoff, double bandwidth);

/* === Double/Debiased Machine Learning (L7: Advanced Econometrics) === */

/* Simplified double ML for ATE estimation.
 * Uses sample splitting + residual-on-residual regression.
 * References: Chernozhukov et al. (2018) Econometrics Journal */
double scm_double_ml_ate(const double* data, int n_rows, int n_cols,
                          int treatment, int outcome,
                          const int* covariates, int n_covariates);

/* === Mediation Analysis for Applications (L7) === */

/* Decompose total effect into direct and indirect effects
 * using the product-of-coefficients method.
 * Returns: nde = natural direct effect, nie = natural indirect effect,
 * total = nde + nie. */
void scm_mediation_decomposition(const SCM_Model* m, int treatment, int mediator,
                                  int outcome, double* nde, double* nie, double* total);

/* === Utility === */

void scm_simpson_result_free(SCM_SimpsonResult* sr);
void scm_simpson_result_print(const SCM_SimpsonResult* sr);
void scm_clinical_result_free(SCM_ClinicalResult* cr);
void scm_clinical_result_print(const SCM_ClinicalResult* cr);
void scm_psm_result_free(SCM_PSMResult* pr);
void scm_psm_result_print(const SCM_PSMResult* pr);

#endif /* SCM_APPLICATIONS_H */
