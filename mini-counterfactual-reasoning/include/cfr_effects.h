#ifndef CFR_EFFECTS_H
#define CFR_EFFECTS_H

#include "cfr_core.h"

/* ============================================================
 * cfr_effects.h -- Causal Effect Estimation
 *
 * Key estimands:
 *   ATE  = E[Y(1) - Y(0)]  (Average Treatment Effect)
 *   ATT  = E[Y(1) - Y(0) | T=1]  (Effect on Treated)
 *   ATU  = E[Y(1) - Y(0) | T=0]  (Effect on Untreated)
 *   CATE = E[Y(1) - Y(0) | X=x]  (Conditional ATE)
 *   ITE  = Y_i(1) - Y_i(0)       (Individual Treatment Effect)
 *
 * Refs: Imbens & Rubin (2015), Hernan & Robins (2020)
 * ============================================================ */

typedef struct {
    double ate;        /* average treatment effect */
    double att;        /* average effect on treated */
    double atu;        /* average effect on untreated */
    double se_ate;     /* standard error of ATE */
    double ci_lower;   /* 95% CI lower bound */
    double ci_upper;   /* 95% CI upper bound */
    int n_treated;
    int n_control;
} CFRCausalEffects;

typedef struct {
    double* cate;      /* conditional ATE for each covariate level */
    double* x_values;  /* covariate values */
    int n_levels;
    double se_cate;
} CFRConditionalEffects;

/* --- ATE Computation --- */
CFRCausalEffects* cfr_eff_create(void);
void              cfr_eff_free(CFRCausalEffects* eff);

/* From potential outcomes (oracle) */
void cfr_eff_ate_oracle(CFRCausalEffects* eff,
                         double* Y0, double* Y1, int n);

/* From observed data (difference in means) */
void cfr_eff_ate_diff_means(CFRCausalEffects* eff,
                             double* treatment, double* outcome, int n);

/* From observed data with regression adjustment */
void cfr_eff_ate_regression(CFRCausalEffects* eff,
                             double* treatment, double* outcome,
                             double* covariates, int n, int p);

/* Inverse Probability Weighting (IPW) */
void cfr_eff_ate_ipw(CFRCausalEffects* eff,
                      double* treatment, double* outcome,
                      double* propensity, int n);

/* Doubly Robust estimation */
void cfr_eff_ate_doubly_robust(CFRCausalEffects* eff,
                                double* treatment, double* outcome,
                                double* propensity,
                                double* Y0_pred, double* Y1_pred, int n);

/* Stratification / Subclassification */
void cfr_eff_ate_stratified(CFRCausalEffects* eff,
                             double* treatment, double* outcome,
                             int* strata, int n, int n_strata);

/* Matching */
void cfr_eff_ate_matching(CFRCausalEffects* eff,
                           double* treatment, double* outcome,
                           double* covariates, int n, int p);

/* --- Conditional Effects --- */
CFRConditionalEffects* cfr_cond_eff_create(int n_levels);
void cfr_cond_eff_free(CFRConditionalEffects* ce);

/* --- Effect Size --- */
double cfr_eff_cohens_d(double ate, double pooled_sd);
double cfr_eff_relative_risk(double* treatment, double* outcome, int n);
double cfr_eff_odds_ratio(double* treatment, double* outcome, int n);

/* --- Sensitivity --- */
void cfr_eff_sensitivity_analysis(double ate, double se, double gamma,
                                   double* lower_bound, double* upper_bound);

/* --- Multiple Testing --- */
void cfr_eff_fwer_bonferroni(double* p_values, int n_tests,
                              double alpha, int* significant);
void cfr_eff_fdr_benjamini_hochberg(double* p_values, int n_tests,
                                     double alpha, int* significant);

/* --- Forest Plot --- */
void cfr_eff_forest_plot(double* strata_ate, double* strata_se,
                          int n_strata, double overall_ate);

/* --- Matching Result --- */
typedef struct { double* weights; int* matched_indices; int n; } CFRMatchingResult;
CFRMatchingResult* cfr_eff_genetic_matching(double* treatment, double* outcome,
    double* covariates, int n, int p);

/* --- Output --- */
void cfr_eff_print(CFRCausalEffects* eff);
void cfr_cond_eff_print(CFRConditionalEffects* ce);

#endif /* CFR_EFFECTS_H */
