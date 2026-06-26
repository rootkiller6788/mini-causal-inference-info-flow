#ifndef EFFECT_ESTIMATORS_H
#define EFFECT_ESTIMATORS_H
#include "adjustment_sets.h"

/*
 * effect_estimators.h — Causal Effect Estimation Methods
 *
 * Once identification is established, we estimate the causal effect
 * from data. Key methods:
 *
 *   1. Propensity Score Matching (PSM)
 *   2. Inverse Probability Weighting (IPW)
 *   3. G-Computation (G-formula)
 *   4. Doubly Robust Estimation (DR)
 *   5. Stratification / Subclassification
 */

typedef struct {
    double ate;         /* average treatment effect              */
    double att;         /* average treatment effect on treated   */
    double se_ate;      /* standard error of ATE                 */
    double se_att;      /* standard error of ATT                 */
    int    n_treated;   /* number of treated units               */
    int    n_control;   /* number of control units               */
    char  *method;      /* estimation method used                */
} CausalEffect;

typedef struct {
    int     n;          int     d;          /* data dimensions   */
    double *X;          /* covariates, n x d, row-major          */
    int    *T;          /* treatment indicator (0 or 1)          */
    double *Y;          /* outcome                                */
} ObservationalData;

/* Lifecycle */
ObservationalData* obs_data_create(int n, int d);
void obs_data_free(ObservationalData* d);

/* Propensity score: e(x) = P(T=1 | X=x) via logistic regression */
double* propensity_score_estimate(const ObservationalData* data);
void propensity_free(double* ps);

/* IPW Estimator: ATE = (1/n) sum_i [ T_i*Y_i/e_i - (1-T_i)*Y_i/(1-e_i) ] */
CausalEffect* ipw_estimate(const ObservationalData* data, const double* ps);

/* G-Computation: fit E[Y|T,X], then ATE = (1/n) sum_i [E[Y|1,X_i] - E[Y|0,X_i]] */
CausalEffect* gcomp_estimate(const ObservationalData* data);

/* Doubly Robust: combines IPW + outcome model. Consistent if either model correct. */
CausalEffect* doubly_robust_estimate(const ObservationalData* data,
                                      const double* ps);

/* Stratification: divide into K strata by propensity score, estimate per stratum. */
CausalEffect* stratification_estimate(const ObservationalData* data,
                                       const double* ps, int K);

/* Simple difference-in-means (naive, assumes no confounding). */
CausalEffect* naive_estimate(const ObservationalData* data);

void causal_effect_free(CausalEffect* ce);

/* Matching Estimator: 1:1 Nearest Neighbor Propensity Score Matching */
CausalEffect* matching_estimate(const ObservationalData* data, const double* ps);

/* ANCOVA / Regression Adjustment (equivalent to G-Computation with linear model) */
CausalEffect* ancova_estimate(const ObservationalData* data);

/* Bootstrap Standard Error for any estimator */
typedef CausalEffect* (*estimator_fn)(const ObservationalData*, const double*);
double bootstrap_se(const ObservationalData* data, const double* ps,
                     estimator_fn est, int B, double* point_est);

/* Standardized Mean Difference for covariate balance checking */
void balance_check_smd(const ObservationalData* data, const double* ps,
                        double* smd, int d);

#define EST_DEFAULT_K  5

#endif
