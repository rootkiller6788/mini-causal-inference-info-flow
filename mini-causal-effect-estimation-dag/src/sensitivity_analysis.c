/*
 * sensitivity_analysis.c — Sensitivity Analysis for Unmeasured Confounding
 *
 * Implements Rosenbaum bounds for matched observational studies,
 * E-values (VanderWeele & Ding, 2017), Cornfield conditions,
 * and bias formulas for quantifying the impact of unmeasured confounders.
 *
 * References:
 *   Rosenbaum, "Observational Studies", 2nd ed, 2002, Springer
 *   VanderWeele & Ding, "Sensitivity Analysis in Observational Research",
 *     Annals of Internal Medicine, 2017
 *   Cornfield et al., "Smoking and Lung Cancer", JNCI, 1959
 *   Ding & VanderWeele, "Sensitivity Analysis Without Assumptions",
 *     Epidemiology, 2016
 *   Arah et al., "Bias Formulas for Sensitivity Analysis", Epidemiology, 2010
 */

#include "sensitivity_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * Rosenbaum Bounds for Paired Observational Studies
 *
 * For N matched pairs, let d_i = Y_{1i} - Y_{0i} be the
 * treatment-minus-control difference. The Wilcoxon signed-rank
 * statistic is:
 *   T = sum_i sgn(d_i) * rank(|d_i|)
 *
 * Under unmeasured confounding with odds ratio Gamma:
 *   P(T_i=1 | X_i) / P(T_i=0 | X_i) = Gamma for all i (hidden bias)
 *
 * The distribution of T under Gamma is bounded by two distributions:
 *   T^+ = sum_i B_i * rank(|d_i|)   where B_i ~ Bernoulli(p^+)
 *   T^- = sum_i B_i * rank(|d_i|)   where B_i ~ Bernoulli(p^-)
 * with p^+ = Gamma/(1+Gamma) and p^- = 1/(1+Gamma).
 *
 * For large N, the normal approximation:
 *   E[T^+] = p^+ * sum_i rank(|d_i|)
 *   Var[T^+] = p^+ * (1-p^+) * sum_i rank(|d_i|)^2
 *
 * The upper bound p-value = 1 - Phi((T - E[T^-]) / sqrt(Var[T^-]))
 * The lower bound p-value = 1 - Phi((T - E[T^+]) / sqrt(Var[T^+]))
 * ================================================================ */

typedef struct { double val; int idx; int sign; } RankedPair;

static int cmp_ranked(const void *a, const void *b) {
    double da = fabs(((const RankedPair*)a)->val);
    double db = fabs(((const RankedPair*)b)->val);
    return (da > db) - (da < db);
}

static double normal_cdf(double z) {
    return 0.5 * (1.0 + erf(z / sqrt(2.0)));
}

RosenbaumBounds* rosenbaum_bounds(const double *y_treated,
                                   const double *y_control,
                                   int n_pairs, double gamma) {
    if (!y_treated || !y_control || n_pairs < 2 || gamma < 1.0)
        return NULL;

    RosenbaumBounds *rb = calloc(1, sizeof(RosenbaumBounds));
    if (!rb) return NULL;
    rb->gamma = gamma;

    RankedPair *pairs = calloc((size_t)n_pairs, sizeof(RankedPair));
    if (!pairs) { free(rb); return NULL; }

    /* Compute differences and sort by absolute value */
    int n_pos = 0;
    for (int i = 0; i < n_pairs; i++) {
        pairs[i].val = y_treated[i] - y_control[i];
        pairs[i].idx = i;
        pairs[i].sign = (pairs[i].val > 0) ? 1 : ((pairs[i].val < 0) ? -1 : 0);
        if (pairs[i].sign > 0) n_pos++;
    }

    /* Sort by absolute difference (ties get average rank) */
    qsort(pairs, (size_t)n_pairs, sizeof(RankedPair), cmp_ranked);

    /* Assign ranks with tie handling */
    double *ranks = calloc((size_t)n_pairs, sizeof(double));
    if (!ranks) { free(pairs); free(rb); return NULL; }

    int i = 0;
    while (i < n_pairs) {
        int j = i;
        while (j < n_pairs &&
               fabs(fabs(pairs[j].val) - fabs(pairs[i].val)) < 1e-12)
            j++;
        double avg_rank = (double)(i + j + 1) / 2.0;
        for (int k = i; k < j; k++)
            ranks[pairs[k].idx] = avg_rank;
        i = j;
    }

    /* Compute Wilcoxon signed-rank statistic T */
    double T_stat = 0.0, sum_rank = 0.0, sum_rank_sq = 0.0;
    for (int i = 0; i < n_pairs; i++) {
        sum_rank += ranks[i];
        sum_rank_sq += ranks[i] * ranks[i];
        if (y_treated[i] > y_control[i])
            T_stat += ranks[i];
    }

    /* Bounds under Gamma */
    double p_plus  = gamma / (1.0 + gamma);   /* upper bound prob */
    double p_minus = 1.0 / (1.0 + gamma);     /* lower bound prob */

    double E_plus  = p_plus * sum_rank;
    double E_minus = p_minus * sum_rank;
    double Var_plus  = p_plus * (1.0 - p_plus) * sum_rank_sq;
    double Var_minus = p_minus * (1.0 - p_minus) * sum_rank_sq;

    double sd_plus  = sqrt(Var_plus);
    double sd_minus = sqrt(Var_minus);

    /* Standardized test statistics */
    double z_plus  = (sd_plus  > 1e-12) ? (T_stat - E_plus)  / sd_plus  : 0.0;
    double z_minus = (sd_minus > 1e-12) ? (T_stat - E_minus) / sd_minus : 0.0;

    /* One-sided p-values (upper bound is more conservative) */
    rb->upper_bound = 1.0 - normal_cdf(z_minus);
    rb->lower_bound = 1.0 - normal_cdf(z_plus);

    if (rb->upper_bound > 1.0) rb->upper_bound = 1.0;
    if (rb->upper_bound < 0.0) rb->upper_bound = 0.0;
    if (rb->lower_bound > 1.0) rb->lower_bound = 1.0;
    if (rb->lower_bound < 0.0) rb->lower_bound = 0.0;

    /* ATE bounds: compute effect estimate under Gamma */
    double ate_obs = 0.0;
    for (int i = 0; i < n_pairs; i++)
        ate_obs += (y_treated[i] - y_control[i]);
    ate_obs /= (double)n_pairs;

    /* Under unmeasured confounding of magnitude Gamma,
     * the true ATE could be shifted. Simple approximation:
     *   ATE_lower = ATE_obs - bias(Gamma)
     *   ATE_upper = ATE_obs + bias(Gamma)
     * where bias(Gamma) = (Gamma - 1) * sd_diff / sqrt(n)
     *
     * More refined: Manski bounds (nonparametric partial identification) */
    double sd_diff = 0.0;
    for (int i = 0; i < n_pairs; i++) {
        double d = (y_treated[i] - y_control[i]) - ate_obs;
        sd_diff += d * d;
    }
    sd_diff = sqrt(sd_diff / (double)(n_pairs - 1));
    double bias = (gamma - 1.0) * sd_diff / sqrt((double)n_pairs);

    rb->lower_ate = ate_obs - bias;
    rb->upper_ate = ate_obs + bias;

    free(ranks);
    free(pairs);
    return rb;
}

void rosenbaum_free(RosenbaumBounds* rb) { free(rb); }

/* Find critical Gamma via binary search.
 * The critical Gamma is the smallest Gamma >= 1 for which the
 * upper bound p-value exceeds alpha (default 0.05). */
double rosenbaum_critical_gamma(const double *y_treated,
                                 const double *y_control, int n_pairs) {
    if (!y_treated || !y_control || n_pairs < 2) return 1.0;

    double lo = 1.0, hi = 10.0;
    const double alpha = 0.05;
    const int max_iter = 30;

    /* First, check if even Gamma=1 is non-significant */
    RosenbaumBounds *rb1 = rosenbaum_bounds(y_treated, y_control, n_pairs, 1.0);
    if (!rb1) return 1.0;
    if (rb1->upper_bound > alpha) {
        rosenbaum_free(rb1);
        return 1.0;
    }
    rosenbaum_free(rb1);

    /* Expand hi until upper_bound > alpha */
    for (int i = 0; i < 10; i++) {
        RosenbaumBounds *rb = rosenbaum_bounds(y_treated, y_control, n_pairs, hi);
        if (!rb) break;
        if (rb->upper_bound > alpha) {
            rosenbaum_free(rb);
            break;
        }
        rosenbaum_free(rb);
        hi *= 2.0;
        if (hi > 1000.0) return 1000.0;
    }

    /* Binary search for critical Gamma */
    for (int iter = 0; iter < max_iter; iter++) {
        double mid = (lo + hi) / 2.0;
        RosenbaumBounds *rb = rosenbaum_bounds(y_treated, y_control, n_pairs, mid);
        if (!rb) break;
        if (rb->upper_bound > alpha) {
            hi = mid;
        } else {
            lo = mid;
        }
        rosenbaum_free(rb);
        if (hi - lo < 0.01) break;
    }

    return (lo + hi) / 2.0;
}

/* ================================================================
 * E-value (VanderWeele & Ding, 2017)
 *
 * The E-value is the minimum strength of association that an
 * unmeasured confounder must have with both treatment and outcome
 * to fully explain away an observed treatment-outcome association.
 *
 * For an observed risk ratio RR (or approximating the ATE on the RR scale):
 *
 *   E-value = RR + sqrt(RR * (RR - 1))   for RR >= 1
 *   E-value = 1/RR + sqrt(1/RR * (1/RR - 1))  for RR < 1
 *
 * For the confidence interval limit, replace RR with the CI limit
 * closer to the null (RR_ci).
 *
 * Interpretation:
 *   E-value = 2.0: an unmeasured confounder would need to be associated
 *     with both treatment and outcome by a risk ratio of at least 2
 *     each to explain away the observed effect.
 *   E-value = 1.0: even weak confounding could explain the effect.
 * ================================================================ */

/* Convert ATE to approximate risk ratio.
 * For continuous Y: RR = exp(ATE / sigma_Y) gives a rough RR scale.
 * For binary Y with rare outcome: RR ≈ OR ≈ exp(ATE).
 * Here we use: RR = exp(|ATE| / max(|ATE|, 0.1)) to map to a sensible scale. */
static double ate_to_rr(double ate, double se_ate) {
    (void)se_ate;
    double abs_ate = fabs(ate);
    if (abs_ate < 1e-10) return 1.0;
    /* Normalise: map to approximately 1.1-5.0 range */
    double scale = abs_ate / (abs_ate + 0.1);
    double rr = 1.0 + 4.0 * scale;
    if (ate < 0.0) rr = 1.0 / rr;
    return rr;
}

EValue* compute_evalue(double ate, double se_ate) {
    if (se_ate < 0.0) return NULL;

    EValue *ev = calloc(1, sizeof(EValue));
    if (!ev) return NULL;

    /* Convert ATE to risk ratio scale */
    double rr = ate_to_rr(ate, se_ate);
    ev->risk_ratio = rr;

    /* E-value for point estimate */
    if (rr >= 1.0) {
        ev->e_value = rr + sqrt(rr * (rr - 1.0));
    } else {
        double rr_inv = 1.0 / rr;
        ev->e_value = rr_inv + sqrt(rr_inv * (rr_inv - 1.0));
    }

    /* E-value for CI limit (using CI closest to null) */
    double ci_upper = ate + 1.96 * se_ate;
    double ci_lower = ate - 1.96 * se_ate;
    double ci_near_null;
    if (ate > 0.0) {
        ci_near_null = ci_lower > 0.0 ? ci_lower : 0.0;
    } else if (ate < 0.0) {
        ci_near_null = ci_upper < 0.0 ? ci_upper : 0.0;
    } else {
        ci_near_null = 0.0;
    }

    double rr_ci = ate_to_rr(ci_near_null, se_ate);
    if (rr_ci >= 1.0) {
        ev->e_value_ci = rr_ci + sqrt(rr_ci * (rr_ci - 1.0));
    } else if (rr_ci > 0.0) {
        double rr_ci_inv = 1.0 / rr_ci;
        ev->e_value_ci = rr_ci_inv + sqrt(rr_ci_inv * (rr_ci_inv - 1.0));
    } else {
        ev->e_value_ci = ev->e_value;
    }

    ev->is_robust = (ev->e_value >= 1.5);

    return ev;
}

void evalue_free(EValue* ev) { free(ev); }

/* ================================================================
 * Cornfield Conditions (Cornfield et al., 1959)
 *
 * Classic sensitivity analysis: to explain away an observed relative
 * risk RR_obs between treatment and outcome, an unmeasured confounder
 * U must have relative risks with both treatment and outcome that
 * satisfy:
 *
 *   RR_UT * RR_UY >= RR_obs
 *
 * where RR_UT is the confounder-treatment RR and RR_UY is the
 * confounder-outcome RR.
 * ================================================================ */

CornfieldResult* cornfield_check(double observed_rr,
                                  double rr_threshold_treatment,
                                  double rr_threshold_outcome) {
    if (observed_rr <= 1.0) return NULL;

    CornfieldResult *cr = calloc(1, sizeof(CornfieldResult));
    if (!cr) return NULL;

    /* Minimum required: RR_UT * RR_UY = observed_rr
     * If both are equal (worst case balance): each >= sqrt(observed_rr) */
    double min_each = sqrt(observed_rr);

    cr->min_rr_ct = min_each;
    cr->min_rr_co = min_each;

    /* Check plausibility against thresholds */
    cr->is_plausible = (min_each > rr_threshold_treatment ||
                        min_each > rr_threshold_outcome);

    /* More refined: if RR_UT is known, then RR_UY >= RR_obs / RR_UT */
    if (rr_threshold_treatment > 0.0) {
        cr->min_rr_co = observed_rr / rr_threshold_treatment;
    }
    if (rr_threshold_outcome > 0.0) {
        cr->min_rr_ct = observed_rr / rr_threshold_outcome;
    }

    return cr;
}

void cornfield_free(CornfieldResult* cr) { free(cr); }

/* ================================================================
 * Confounding Bias Formula
 *
 * Bias = delta * gamma
 *
 * where:
 *   delta = E[Y|U=1,X,T] - E[Y|U=0,X,T]  (confounder-outcome effect)
 *   gamma = P(U=1|T=1,X) - P(U=1|T=0,X)  (imbalance in confounder)
 *
 * VanderWeele & Arah, "Bias formulas for sensitivity analysis",
 *   Epidemiology, 2011
 * ================================================================ */

double confounding_bias(double delta, double gamma, double scale) {
    double bias = delta * gamma;
    if (scale <= 0.0) {
        /* Ratio scale: bias multiplicative */
        return bias;
    }
    /* Additive scale: bias additive */
    return bias;
}

double adjusted_ate(double ate_observed, double delta, double gamma) {
    return ate_observed - delta * gamma;
}

/* ================================================================
 * Sensitivity Contour
 * ================================================================ */

SensitivityContour* sensitivity_contour_create(double ate_observed,
                                                double delta_max,
                                                double gamma_max,
                                                int n_points) {
    if (n_points < 2 || delta_max <= 0.0 || gamma_max <= 0.0)
        return NULL;

    SensitivityContour *sc = calloc(1, sizeof(SensitivityContour));
    if (!sc) return NULL;

    sc->n_points = n_points;
    sc->nullify_contour = ate_observed;

    sc->delta_grid = calloc((size_t)n_points, sizeof(double));
    sc->gamma_grid = calloc((size_t)n_points, sizeof(double));
    sc->ate_matrix = calloc((size_t)n_points * n_points, sizeof(double));
    if (!sc->delta_grid || !sc->gamma_grid || !sc->ate_matrix) {
        sensitivity_contour_free(sc);
        return NULL;
    }

    for (int i = 0; i < n_points; i++) {
        sc->delta_grid[i] = delta_max * (double)i / (double)(n_points - 1);
        sc->gamma_grid[i] = gamma_max * (double)i / (double)(n_points - 1);
    }

    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n_points; j++) {
            sc->ate_matrix[i * n_points + j] =
                adjusted_ate(ate_observed, sc->delta_grid[i], sc->gamma_grid[j]);
        }
    }

    return sc;
}

void sensitivity_contour_free(SensitivityContour* sc) {
    if (!sc) return;
    free(sc->delta_grid);
    free(sc->gamma_grid);
    free(sc->ate_matrix);
    free(sc);
}

int sensitivity_nullify_line(const SensitivityContour* sc,
                              double *delta_line, double *gamma_line,
                              int max_points) {
    if (!sc || !delta_line || !gamma_line || max_points <= 0)
        return -1;

    int count = 0;

    /* Find combinations (delta, gamma) that yield ATE = 0.
     * Since ATE_adj = ATE_obs - delta*gamma, we need:
     *   delta * gamma = ATE_obs
     * So gamma = ATE_obs / delta for delta > 0. */
    for (int i = 1; i < sc->n_points && count < max_points; i++) {
        double delta = sc->delta_grid[i];
        if (delta < 1e-10) continue;
        double gamma = sc->nullify_contour / delta;
        if (gamma >= 0.0 && gamma <= sc->gamma_grid[sc->n_points - 1]) {
            delta_line[count] = delta;
            gamma_line[count] = gamma;
            count++;
        }
    }

    return count;
}
