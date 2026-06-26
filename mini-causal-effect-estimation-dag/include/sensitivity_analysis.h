#ifndef SENSITIVITY_ANALYSIS_H
#define SENSITIVITY_ANALYSIS_H
#include "effect_estimators.h"

/*
 * sensitivity_analysis.h — Sensitivity Analysis for Unmeasured Confounding
 *
 * Assesses how sensitive causal effect estimates are to unmeasured
 * confounding. Key methods:
 *
 *   1. Rosenbaum Bounds (Rosenbaum, 2002): for matched observational studies
 *   2. E-value (VanderWeele & Ding, 2017): minimum strength of association
 *      an unmeasured confounder must have to explain away an effect
 *   3. Cornfield Conditions (Cornfield et al., 1959): classical sensitivity
 *   4. Bias formulas (conditional on confounder strength)
 *
 * References:
 *   Rosenbaum, "Observational Studies", 2nd ed, 2002, Springer
 *   VanderWeele & Ding, "Sensitivity Analysis in Observational Research",
 *     Annals of Internal Medicine, 2017
 *   Cornfield et al., "Smoking and Lung Cancer", JNCI, 1959
 *   Ding & VanderWeele, "Sensitivity Analysis Without Assumptions",
 *     Epidemiology, 2016
 */

typedef struct {
    double gamma;           /* Rosenbaum's Gamma: odds ratio of treatment
                               assignment due to unmeasured confounder  */
    double lower_bound;     /* lower bound of p-value under Gamma       */
    double upper_bound;     /* upper bound of p-value under Gamma       */
    double lower_ate;       /* lower bound of ATE confidence interval   */
    double upper_ate;       /* upper bound of ATE confidence interval   */
    double critical_gamma;  /* Gamma at which result becomes non-sig    */
} RosenbaumBounds;

typedef struct {
    double e_value;         /* E-value for the point estimate            */
    double e_value_ci;      /* E-value for the confidence interval limit */
    double risk_ratio;      /* observed risk ratio / effect size         */
    bool   is_robust;       /* whether effect is robust at E-value > 1.5 */
} EValue;

/* ---- Rosenbaum Bounds ---- */

/* Compute Rosenbaum bounds for a paired observational study.
 * Given matched-pair data, compute the sensitivity of the Wilcoxon
 * signed-rank test to unmeasured confounding.
 *
 * For Gamma = 1, this is the standard test (no unmeasured confounding).
 * As Gamma increases, the bounds widen until the result is non-significant.
 *
 * Parameters:
 *   y_treated[n]  — outcomes for treated units in matched pairs
 *   y_control[n]  — outcomes for control units in matched pairs
 *   n_pairs       — number of matched pairs
 *   gamma         — odds ratio of differential treatment assignment */
RosenbaumBounds* rosenbaum_bounds(const double *y_treated,
                                   const double *y_control,
                                   int n_pairs, double gamma);
void rosenbaum_free(RosenbaumBounds* rb);

/* Find the critical Gamma at which the result becomes non-significant
 * (upper bound p-value > 0.05). Returns the maximum Gamma for which
 * the result remains significant. */
double rosenbaum_critical_gamma(const double *y_treated,
                                 const double *y_control, int n_pairs);

/* ---- E-value ---- */

/* Compute the E-value for an observed causal effect.
 *
 * The E-value is the minimum strength of association, on the risk ratio
 * scale, that an unmeasured confounder would need to have with both
 * the treatment and the outcome to fully explain away the observed
 * treatment-outcome association, conditional on measured covariates.
 *
 * Formula: E-value = RR + sqrt(RR * (RR - 1))
 *   where RR = exp(ATE) for continuous outcomes transformed to RR scale.
 *
 * Interpretation:
 *   E-value > 2.0  => robust to fairly strong unmeasured confounding
 *   E-value < 1.5  => sensitive to modest unmeasured confounding
 *   E-value ~ 1.0  => very fragile, even weak confounding could explain it
 *
 * VanderWeele & Ding, Annals of Internal Medicine, 2017 */
EValue* compute_evalue(double ate, double se_ate);
void evalue_free(EValue* ev);

/* ---- Cornfield Conditions ---- */

/* Check Cornfield's conditions: an unmeasured confounder must be
 * associated with both treatment and outcome more strongly than
 * the observed treatment-outcome association to explain it away.
 *
 * Returns the minimum relative risk between confounder-treatment
 * and confounder-outcome needed to nullify the observed effect. */
typedef struct {
    double min_rr_ct;    /* minimum RR: confounder -> treatment           */
    double min_rr_co;    /* minimum RR: confounder -> outcome             */
    bool   is_plausible; /* whether these RRs are plausible in context   */
} CornfieldResult;

CornfieldResult* cornfield_check(double observed_rr,
                                  double rr_threshold_treatment,
                                  double rr_threshold_outcome);
void cornfield_free(CornfieldResult* cr);

/* ---- Bias Formulas ---- */

/* Compute bias of ATE due to unmeasured confounder U.
 *
 * Bias formula (VanderWeele & Arah, 2011):
 *   Bias = delta * gamma
 *   where delta = E[Y|U=1,X,T] - E[Y|U=0,X,T] (confounder-outcome effect)
 *         gamma = P(U=1|T=1,X) - P(U=1|T=0,X) (imbalance in U)
 *
 * Parameters:
 *   delta — effect of U on Y (outcome), on the outcome scale
 *   gamma — imbalance of U between treatment groups (difference in prevalence)
 *   scale — 1.0 for additive scale, 0.0 for ratio scale */
double confounding_bias(double delta, double gamma, double scale);

/* Compute the adjusted ATE after accounting for unmeasured confounding.
 *   ATE_adjusted = ATE_observed - Bias
 * If the unmeasured confounder is beneficial (positive confounding),
 * the true effect is smaller than observed. */
double adjusted_ate(double ate_observed, double delta, double gamma);

/* ---- Sensitivity Contour ---- */

/* Generate a sensitivity contour: for a grid of (delta, gamma) values,
 * compute the adjusted ATE and determine which combinations would
 * nullify the observed effect. */
typedef struct {
    int      n_points;    /* number of grid points per dimension         */
    double  *delta_grid;  /* delta values grid                          */
    double  *gamma_grid;  /* gamma values grid                          */
    double  *ate_matrix;  /* adjusted ATE matrix, row-major             */
    double   nullify_contour; /* threshold below which effect is null   */
} SensitivityContour;

SensitivityContour* sensitivity_contour_create(double ate_observed,
                                                double delta_max,
                                                double gamma_max,
                                                int n_points);
void sensitivity_contour_free(SensitivityContour* sc);

/* Find the (delta, gamma) combination at which ATE becomes exactly zero. */
int sensitivity_nullify_line(const SensitivityContour* sc,
                              double *delta_line, double *gamma_line,
                              int max_points);

#endif
