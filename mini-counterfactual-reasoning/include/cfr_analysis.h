#ifndef CFR_ANALYSIS_H
#define CFR_ANALYSIS_H

#include "cfr_core.h"

/* ============================================================
 * cfr_analysis.h -- Counterfactual Analysis Methods
 *
 * Implements Pearl's causal inference tools:
 * - Probability of necessity/sufficiency
 * - E-value sensitivity analysis
 * - IPW estimation
 * - IV estimation
 * - Front-door adjustment
 * - G-computation
 * - Manski bounds
 * - Monotonicity checks
 * ============================================================ */

/* Probability of necessity: PN = P(Y_x'=0 | X=1, Y=1) */
double cfr_probability_of_necessity(double p_y1_do_x0, double p_x1_y1,
                                     double p_y1_do_x1);

/* Probability of sufficiency: PS = P(Y_x=1 | X=0, Y=0) */
double cfr_probability_of_sufficiency(double p_y1_do_x1,
                                       double p_y1_do_x0,
                                       double p_x0_y0);

/* Individual treatment effect */
double cfr_individual_treatment_effect(double y_treated, double y_control);

/* ATT estimate */
double cfr_att_estimate(const double* y_outcomes, const int* treatments,
                         int n_samples);

/* Mediation analysis result */
typedef struct { double nde; double nie; double te; } CFRMediationResult;
CFRMediationResult cfr_mediation_analysis(double e_y1m0, double e_y0m0,
                                           double e_y1m1);

/* Propensity score matching */
int cfr_propensity_match(const double* propensity_scores, int n_treated,
                          const int* treatment, int* matched_pairs);

/* E-value for unmeasured confounding sensitivity */
double cfr_e_value(double risk_ratio);

/* IPW estimator */
double cfr_ipw_estimate(const double* outcomes, const double* propensity,
                         const int* treatment, int n_samples);

/* Manski bounds */
typedef struct { double lower; double upper; } CFRBoundsAna;
CFRBoundsAna cfr_manski_bounds(double e_y_given_t1, double e_y_given_t0,
                                 double p_t1);

/* Front-door adjustment formula */
double cfr_front_door_adjustment(const double* p_z_given_x,
                                  const double* p_y_given_xz,
                                  const double* p_x, int n_z, int n_x);

/* IV estimator: beta_IV = Cov(Z,Y) / Cov(Z,X) */
double cfr_iv_estimate(const double* z, const double* x, const double* y,
                        int n_samples);

/* G-computation formula */
double cfr_g_computation(const double* outcomes, const int* treatments,
                          const double* covariates, int n_samples,
                          int n_covariates);

/* Monotonicity check */
int cfr_check_monotonicity(const int* treatment_assigned,
                            const int* treatment_received, int n);

/* Placebo refutation test */
double cfr_placebo_refutation(const double* outcomes, const int* true_treatment,
                               const int* placebo_treatment, int n_samples);

/* Random common cause test */
double cfr_random_common_cause_test(double effect_with, double effect_without,
                                     double tol);

/* Data fusion: weighted combination */
double cfr_fuse_estimates(const double* estimates, const double* weights, int n);

#endif /* CFR_ANALYSIS_H */
