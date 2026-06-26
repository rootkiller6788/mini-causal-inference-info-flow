#ifndef CFR_BOUNDS_H
#define CFR_BOUNDS_H

#include "cfr_core.h"

/* ============================================================
 * cfr_bounds.h -- Counterfactual Probability Bounds
 *
 * Probabilities of Causation (Pearl, 1999):
 *
 *   PN  = P(Y_x' = y' | X=x, Y=y)  (Probability of Necessity)
 *   PS  = P(Y_x  = y  | X=x', Y=y') (Probability of Sufficiency)
 *   PNS = P(Y_x = y, Y_x' = y')    (Probability of Necessity and Sufficiency)
 *
 * Bounds under different assumptions:
 * - No assumptions: trivial bounds [0, 1]
 * - Exogeneity: experimental bounds
 * - Monotonicity: tighter bounds
 * - Exogeneity + Monotonicity: tightest non-trivial bounds
 *
 * Refs: Pearl (1999), Tian & Pearl (2000), Mueller & Pearl (2022)
 * ============================================================ */

typedef struct {
    double pn_lower, pn_upper;    /* Probability of Necessity */
    double ps_lower, ps_upper;    /* Probability of Sufficiency */
    double pns_lower, pns_upper;  /* P(Necessity and Sufficiency) */
    double ate;                    /* ATE for context */
    bool identifiable;             /* exact identification possible */
} CFRBounds;

typedef struct {
    double* outcomes;      /* observed outcomes */
    double* treatments;    /* treatment assignments */
    int n;
    bool randomized;       /* data from RCT */
    bool monotonicity;     /* monotonicity assumption holds */
} CFRBoundsData;

/* --- Lifecycle --- */
CFRBounds*     cfr_bounds_create(void);
void           cfr_bounds_free(CFRBounds* b);
CFRBoundsData* cfr_bounds_data_create(int n);
void           cfr_bounds_data_free(CFRBoundsData* data);

/* --- Bound Computation --- */
void cfr_bounds_trivial(CFRBoundsData* data, CFRBounds* b);
void cfr_bounds_experimental(CFRBoundsData* data, CFRBounds* b);
void cfr_bounds_with_monotonicity(CFRBoundsData* data, CFRBounds* b);
void cfr_bounds_tian_pearl(CFRBoundsData* data, CFRBounds* b);

/* --- Individual-Level Bounds --- */
void cfr_bounds_individual(double y1_obs, double y0_obs,
                            double* pn_bounds, double* ps_bounds);

/* --- Proportion of Causation --- */
double cfr_bounds_excess_fraction(CFRBoundsData* data);

/* --- Output --- */
void cfr_bounds_print(CFRBounds* b);
void cfr_bounds_data_print(CFRBoundsData* data);

/* --- ATE Bounds --- */
double cfr_bounds_ate_bounds(double* treatment, double* outcome,
    int n, double y_min, double y_max, double* lower, double* upper);
void cfr_bounds_quantile_treatment_effect(
    double* treatment, double* outcome, int n,
    double q, double* qte_lower, double* qte_upper);
void cfr_bounds_sharp_bounds(double* y1_obs, double* y0_obs,
    int n1, int n0, double* lower, double* upper);
int cfr_bounds_sensitivity(double pn, double p_y_given_x,
    double sensitivity_param, double* adjusted_pn);
void cfr_bounds_plot_bounds(CFRBounds* b, double* assumptions,
    int n_assumptions, double* pn_curve, double* ps_curve);

/* --- Advanced Bound Methods --- */
int  cfr_bounds_balke_pearl_iv(double p_y1x1_z1, double p_y1x1_z0,
                                double p_y1x0_z1, double p_y1x0_z0,
                                double* lower, double* upper);
void cfr_bounds_mediation_sensitivity(double nde, double nie, double rho,
                                       double* nde_adj, double* nie_adj);
double cfr_bounds_attributable_fraction_interval(CFRBoundsData* data,
                                                  double* af_lower, double* af_upper);
void cfr_bounds_stochastic_dominance(double* y1_dist, int n1,
                                      double* y0_dist, int n0, int* dominates);
void cfr_bounds_partial_identification_map(double* y1_marginal, double* y0_marginal,
                                            int n, double* ate_bounds);
void cfr_bounds_cross_world_counterfactual(double p_y1x0, double p_y0x1,
                                            double* pn_lower, double* pn_upper);
void cfr_bounds_manski_trimming(double* outcome, double* treatment, int n,
                                 double trim_pct, double* trimmed_ate);
void cfr_bounds_horowitz_manski_adaptive(double* outcome, double* treatment,
                                          int n, double* ate, double* ci_lower, double* ci_upper);

#endif /* CFR_BOUNDS_H */
