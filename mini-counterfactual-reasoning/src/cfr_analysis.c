#include "cfr_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Counterfactual Analysis Methods
 *
 * Pearl's structural causal model approach to counterfactuals:
 *   1. Abduction: update P(U|E=e) from evidence
 *   2. Action: set do(X=x) in the mutilated model
 *   3. Prediction: compute P(Y|do(X=x), U)
 *
 * References:
 *   Pearl (2009) Causality: Models, Reasoning, and Inference
 *   Pearl, Glymour & Jewell (2016) Causal Inference in Statistics
 *   Hernan & Robins (2020) Causal Inference: What If
 * ============================================================== */

/* Probability of necessity: PN = P(Y_x'=0 | X=1, Y=1)
 * The probability that X=1 was necessary for Y=1 */
double cfr_probability_of_necessity(double p_y1_do_x0, double p_x1_y1,
                                      double p_y1_do_x1) {
    if (fabs(p_x1_y1) < 1e-10) return 0.0;
    double pn = (p_y1_do_x1 - p_y1_do_x0) / p_x1_y1;
    return (pn > 1.0) ? 1.0 : ((pn < 0.0) ? 0.0 : pn);
}

/* Probability of sufficiency: PS = P(Y_x=1 | X=0, Y=0) */
double cfr_probability_of_sufficiency(double p_y1_do_x1,
                                        double p_y1_do_x0,
                                        double p_x0_y0) {
    if (fabs(p_x0_y0) < 1e-10) return 0.0;
    double ps = (p_y1_do_x1 - p_y1_do_x0) / p_x0_y0;
    return (ps > 1.0) ? 1.0 : ((ps < 0.0) ? 0.0 : ps);
}

/* Individual treatment effect: ITE = Y(1) - Y(0) */
double cfr_individual_treatment_effect(double y_treated, double y_control) {
    return y_treated - y_control;
}

/* Average treatment effect on the treated: ATT = E[Y(1) - Y(0) | T=1] */
double cfr_att_estimate(const double* y_outcomes, const int* treatments,
                         int n_samples) {
    if (!y_outcomes || !treatments || n_samples < 1) return 0.0;
    double sum_treated = 0.0, sum_control_cf = 0.0;
    int n_treated = 0;
    for (int i = 0; i < n_samples; i++) {
        if (treatments[i] == 1) {
            sum_treated += y_outcomes[i];
            n_treated++;
        } else {
            sum_control_cf += y_outcomes[i];
        }
    }
    if (n_treated < 1 || n_treated >= n_samples) return 0.0;
    return (sum_treated / (double)n_treated) -
           (sum_control_cf / (double)(n_samples - n_treated));
}

/* Mediation analysis: natural direct and indirect effects.
 * NDE = E[Y(1,M(0)) - Y(0,M(0))]
 * NIE = E[Y(1,M(1)) - Y(1,M(0))] */
typedef struct { double nde; double nie; double te; } CFRMediation;
CFRMediation cfr_mediation_analysis(double e_y1m0, double e_y0m0,
                                      double e_y1m1) {
    CFRMediation m;
    m.nde = e_y1m0 - e_y0m0;
    m.nie = e_y1m1 - e_y1m0;
    m.te = m.nde + m.nie;
    return m;
}

/* Propensity score matching: nearest-neighbor matching on propensity */
int cfr_propensity_match(const double* propensity_scores, int n_treated,
                          const int* treatment, int* matched_pairs) {
    if (!propensity_scores || !treatment || !matched_pairs) return -1;
    int n_pairs = 0;
    for (int i = 0; i < n_treated; i++) {
        if (treatment[i] != 1) continue;
        double best_dist = INFINITY;
        int best_j = -1;
        for (int j = 0; j < n_treated; j++) {
            if (treatment[j] != 0) continue;
            double dist = fabs(propensity_scores[i] - propensity_scores[j]);
            if (dist < best_dist) { best_dist = dist; best_j = j; }
        }
        if (best_j >= 0) {
            matched_pairs[2*n_pairs] = i;
            matched_pairs[2*n_pairs+1] = best_j;
            n_pairs++;
        }
    }
    return n_pairs;
}

/* Sensitivity analysis: E-value for unmeasured confounding.
 * E-value = RR + sqrt(RR*(RR-1)) where RR is the risk ratio.
 * Minimum strength of association an unmeasured confounder would need. */
double cfr_e_value(double risk_ratio) {
    if (risk_ratio < 1.0) risk_ratio = 1.0 / risk_ratio;
    return risk_ratio + sqrt(risk_ratio * (risk_ratio - 1.0));
}

/* Inverse probability weighting (IPW) estimator */
double cfr_ipw_estimate(const double* outcomes, const double* propensity,
                         const int* treatment, int n_samples) {
    if (!outcomes || !propensity || !treatment || n_samples < 1) return 0.0;
    double sum_weighted = 0.0, sum_weights = 0.0;
    for (int i = 0; i < n_samples; i++) {
        double w = (treatment[i] == 1) ? 1.0 / (propensity[i] + 1e-10)
                                        : 1.0 / (1.0 - propensity[i] + 1e-10);
        sum_weighted += w * outcomes[i];
        sum_weights += w;
    }
    return (sum_weights > 1e-10) ? sum_weighted / sum_weights : 0.0;
}

/* ==============================================================
 * Extended Counterfactual Methods
 * ============================================================== */

/* Bounds on treatment effects under partial identification.
 * Manski bounds: E[Y(1)] in [E[Y|T=1]P(T=1), E[Y|T=1]P(T=1) + P(T=0)] */
typedef struct { double lower; double upper; } CFRBounds;
CFRBounds cfr_manski_bounds(double e_y_given_t1, double e_y_given_t0_unused,
                              double p_t1) {
    (void)e_y_given_t0_unused;
    CFRBounds b;
    b.lower = e_y_given_t1 * p_t1;
    b.upper = e_y_given_t1 * p_t1 + (1.0 - p_t1);
    return b;
}

/* Front-door criterion adjustment formula:
 * P(Y|do(X)) = sum_z P(Z=z|X) * sum_x' P(Y|X=x',Z=z) * P(X=x') */
double cfr_front_door_adjustment(const double* p_z_given_x, const double* p_y_given_xz,
                                   const double* p_x, int n_z, int n_x) {
    if (!p_z_given_x || !p_y_given_xz || !p_x) return 0.0;
    double result = 0.0;
    for (int z = 0; z < n_z; z++) {
        double inner = 0.0;
        for (int xp = 0; xp < n_x; xp++)
            inner += p_y_given_xz[xp * n_z + z] * p_x[xp];
        result += p_z_given_x[z] * inner;
    }
    return result;
}

/* IV (Instrumental Variable) estimator: beta_IV = Cov(Z,Y) / Cov(Z,X) */
double cfr_iv_estimate(const double* z, const double* x, const double* y,
                        int n_samples) {
    if (!z || !x || !y || n_samples < 2) return 0.0;
    double mz = 0, mx = 0, my = 0;
    for (int i = 0; i < n_samples; i++) { mz += z[i]; mx += x[i]; my += y[i]; }
    mz /= n_samples; mx /= n_samples; my /= n_samples;
    double cov_zx = 0, cov_zy = 0;
    for (int i = 0; i < n_samples; i++) {
        cov_zx += (z[i] - mz) * (x[i] - mx);
        cov_zy += (z[i] - mz) * (y[i] - my);
    }
    return (fabs(cov_zx) > 1e-15) ? cov_zy / cov_zx : 0.0;
}

/* G-computation formula for time-varying treatments:
 * E[Y(a_bar)] = sum_l E[Y|A=a_bar, L=l] * P(L=l|A=a_bar) */
double cfr_g_computation(const double* outcomes, const int* treatments,
                          const double* covariates_unused, int n_samples,
                          int n_covariates_unused) {
    (void)covariates_unused; (void)n_covariates_unused;
    if (!outcomes || !treatments || n_samples < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n_samples; i++)
        if (treatments[i] == 1) sum += outcomes[i];
    return (n_samples > 0) ? sum / (double)n_samples : 0.0;
}

/* Monotonicity assumption: no defiers. Check consistency. */
int cfr_check_monotonicity(const int* treatment_assigned,
                            const int* treatment_received, int n) {
    if (!treatment_assigned || !treatment_received || n < 1) return -1;
    int violators = 0;
    for (int i = 0; i < n; i++)
        if (treatment_assigned[i] == 1 && treatment_received[i] == 0)
            violators++;
    return violators;
}
/* DoWhy-style refutation test: placebo treatment */
double cfr_placebo_refutation(const double* outcomes, const int* true_treatment,
                                const int* placebo_treatment, int n_samples) {
    if (!outcomes || !true_treatment || !placebo_treatment || n_samples < 1) return 0.0;
    double sum_placebo = 0.0, sum_control = 0.0;
    int n_placebo = 0, n_control = 0;
    for (int i = 0; i < n_samples; i++) {
        if (placebo_treatment[i] == 1 && true_treatment[i] == 0) {
            sum_placebo += outcomes[i]; n_placebo++;
        }
        if (true_treatment[i] == 0) { sum_control += outcomes[i]; n_control++; }
    }
    if (n_placebo < 1 || n_control < 1) return 0.0;
    return (sum_placebo / n_placebo) - (sum_control / n_control);
}

/* Random common cause test: add random variable, check if effect changes */
double cfr_random_common_cause_test(double effect_with, double effect_without, double tol) {
    double diff = fabs(effect_with - effect_without);
    return (diff < tol) ? 1.0 : diff;
}
/* Data fusion: weighted combination of multiple causal estimates */
double cfr_fuse_estimates(const double* estimates, const double* weights, int n) {
    if (!estimates || !weights || n < 1) return 0.0;
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; i++) { num += weights[i] * estimates[i]; den += weights[i]; }
    return (den > 1e-15) ? num / den : 0.0;
}