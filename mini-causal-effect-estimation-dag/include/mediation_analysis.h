#ifndef MEDIATION_ANALYSIS_H
#define MEDIATION_ANALYSIS_H
#include "effect_estimators.h"

/*
 * mediation_analysis.h — Causal Mediation Analysis
 *
 * Decomposes a total causal effect into direct and indirect pathways
 * through one or more mediators. Key concepts:
 *
 *   1. Controlled Direct Effect (CDE)
 *   2. Natural Direct Effect (NDE)
 *   3. Natural Indirect Effect (NIE)
 *   4. Total Effect = NDE + NIE
 *
 * References:
 *   Baron & Kenny, JPSP, 1986
 *   Imai, Keele, Tingley, Psychological Methods, 2010
 *   VanderWeele, "Explanation in Causal Inference", 2015
 *   Pearl, "Causality", 2009, Ch.4
 */

typedef struct {
    double a0;          double a1;          double *a_cov;    int n_cov;
    double b0;          double b1;          double b2;         double *b_cov;
    double nie;         double nde;         double te;         double pm;
    double se_nie;      double se_nde;      double se_te;      double se_pm;
    double ci_lower_nie; double ci_upper_nie;
    double ci_lower_nde; double ci_upper_nde;
    double ci_lower_pm;  double ci_upper_pm;
    bool   nie_significant; bool nde_significant;
} MediationResult;

MediationResult* mediation_baron_kenny(const int *T, const double *M,
                                        const double *Y, const double *X,
                                        int n, int d);
void mediation_result_free(MediationResult* mr);

MediationResult* mediation_ikt(const int *T, const double *M,
                                const double *Y, const double *X,
                                int n, int d, int n_sim);

double controlled_direct_effect(const MediationResult* mr, double m_value);

double proportion_mediated(double nie, double nde);

double joint_indirect_effect(const double *a_coeff,
                              const double *b_coeff,
                              int n_mediators);

typedef struct {
    double *rho_grid;     double *acme_lower;   double *acme_upper;
    int     n_rho;        double  rho_at_zero;
} MediationSensitivity;

MediationSensitivity* mediation_sensitivity(const MediationResult* mr,
                                             int n_rho);
void mediation_sensitivity_free(MediationSensitivity* ms);

MediationResult* mediation_binary_mediator(const int *T, const int *M,
                                            const double *Y, const double *X,
                                            int n, int d);

typedef struct {
    double *cde_at_time;  double *nie_cumulative;  int n_timepoints;
} LongitudinalMediation;

LongitudinalMediation* longitudinal_mediation(const double *T,
                                               const double *M,
                                               const double *Y,
                                               int n, int n_timepoints);
void long_med_free(LongitudinalMediation* lm);

#endif
