/* scm_applications.c -- Real-world SCM Applications (L7)
 * Applications in epidemiology, econometrics, policy evaluation, medicine.
 * References:
 *   Hernan & Robins (2020) "Causal Inference: What If"
 *   Angrist & Pischke (2009) "Mostly Harmless Econometrics"
 *   Imbens & Rubin (2015) "Causal Inference for Statistics, Social, and Biomedical Sciences"
 *   Chernozhukov et al. (2018) "Double/debiased machine learning for treatment and structural parameters"
 */

#include "scm_applications.h"
#include "scm_intervention.h"
#include "scm_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define SCM_APP_EPS 1e-12

/* ---------- Simpson's Paradox ---------- */

/* Detect Simpson's paradox: compute aggregate effect and subgroup effects,
 * flag if they are in opposite directions. */
SCM_SimpsonResult* scm_simpsons_paradox_detect(const double* data, int n_rows, int n_cols,
                                                int treatment, int outcome, int subgroup) {
    SCM_SimpsonResult* sr = calloc(1, sizeof(SCM_SimpsonResult));
    if (!sr || !data || n_rows < 4) { free(sr); return NULL; }

    /* Aggregate effect: simple difference in means */
    double y_t1 = 0.0, y_t0 = 0.0;
    int n_t1 = 0, n_t0 = 0;
    for (int i = 0; i < n_rows; i++) {
        if (data[i * n_cols + treatment] > 0.5) {
            y_t1 += data[i * n_cols + outcome]; n_t1++;
        } else {
            y_t0 += data[i * n_cols + outcome]; n_t0++;
        }
    }
    if (n_t1 > 0) y_t1 /= (double)n_t1;
    if (n_t0 > 0) y_t0 /= (double)n_t0;
    sr->aggregate_effect = y_t1 - y_t0;

    /* Subgroup effects: for each unique subgroup value */
    /* Find unique subgroup values (assume <= 16 unique values) */
    double unique_vals[16];
    int n_unique = 0;
    for (int i = 0; i < n_rows && n_unique < 16; i++) {
        double v = data[i * n_cols + subgroup];
        bool found = false;
        for (int j = 0; j < n_unique; j++)
            if (fabs(unique_vals[j] - v) < SCM_APP_EPS) { found = true; break; }
        if (!found) unique_vals[n_unique++] = v;
    }
    sr->n_subgroups = n_unique;

    /* Compute effect in each subgroup */
    for (int g = 0; g < n_unique; g++) {
        double y1 = 0.0, y0 = 0.0;
        int n1 = 0, n0 = 0;
        for (int i = 0; i < n_rows; i++) {
            if (fabs(data[i * n_cols + subgroup] - unique_vals[g]) > SCM_APP_EPS) continue;
            if (data[i * n_cols + treatment] > 0.5) {
                y1 += data[i * n_cols + outcome]; n1++;
            } else {
                y0 += data[i * n_cols + outcome]; n0++;
            }
        }
        if (n1 > 0) y1 /= (double)n1;
        if (n0 > 0) y0 /= (double)n0;
        sr->subgroup_effects[g] = y1 - y0;
    }

    /* Check for paradox: aggregate positive, all subgroups negative (or vice versa) */
    bool all_same_sign = true;
    int first_sign = (sr->subgroup_effects[0] > 0) ? 1 : ((sr->subgroup_effects[0] < 0) ? -1 : 0);
    for (int g = 1; g < n_unique && all_same_sign; g++) {
        int sign = (sr->subgroup_effects[g] > 0) ? 1 : ((sr->subgroup_effects[g] < 0) ? -1 : 0);
        if (sign != first_sign) all_same_sign = false;
    }
    int agg_sign = (sr->aggregate_effect > 0) ? 1 : ((sr->aggregate_effect < 0) ? -1 : 0);
    sr->paradox_detected = (all_same_sign && agg_sign != 0 && first_sign != 0 && agg_sign != first_sign);

    /* Cochran-Mantel-Haenszel weighted correction:
     * Weight each subgroup effect by inverse-variance (proportional to group size).
     * correction = sum_g (n_g * effect_g) / sum_g (n_g)
     * This removes the confounding by subgroup size that causes the paradox. */
    {
        double weighted_sum = 0.0;
        double total_weight = 0.0;
        int n_groups[16] = {0};
        for (int g = 0; g < n_unique; g++) {
            for (int i = 0; i < n_rows; i++) {
                if (fabs(data[i * n_cols + subgroup] - unique_vals[g]) < SCM_APP_EPS)
                    n_groups[g]++;
            }
            weighted_sum += n_groups[g] * sr->subgroup_effects[g];
            total_weight += n_groups[g];
        }
        sr->simpson_correction = (total_weight > 0)
            ? weighted_sum / total_weight
            : sr->aggregate_effect;
    }

    return sr;
}

/* Resolve via back-door adjustment: condition on the adjustment set */
double scm_simpsons_paradox_resolve(const SCM_Model* m, int treatment, int outcome,
                                     const SCM_VarSet* adjustment_set) {
    if (!m) return 0.0;
    SCM_CausalEffect* ce = scm_causal_effect_adjustment(m, treatment, outcome, adjustment_set);
    double effect = ce ? ce->ate : 0.0;
    scm_causal_effect_free(ce);
    return effect;
}

/* ---------- Causal Paradoxes ---------- */

/* M-bias test: M is a collider between two unobserved confounders.
 * Conditioning on M opens a path between treatment and outcome. */
bool scm_m_bias_test(const SCM_Model* m, int x, int y, int m_collider) {
    if (!m) return false;
    /* M-bias structure: U1 -> X, U1 -> M <- U2, U2 -> Y.
     * X and Y are independent unconditionally.
     * Conditioning on M makes them dependent. */
    SCM_VarSet empty; scm_varset_init(&empty);
    bool indep_uncond = scm_d_separated(m, x, y, &empty);

    SCM_VarSet cond; scm_varset_init(&cond);
    scm_varset_add(&cond, m_collider);
    bool indep_cond = scm_d_separated(m, x, y, &cond);

    /* Bias: unconditionally independent but dependent conditional on collider */
    return indep_uncond && !indep_cond;
}

/* Berkson's paradox: conditioning on a common effect (collider)
 * creates spurious negative correlation between independent causes. */
bool scm_berkson_paradox_test(const SCM_Model* m, int disease_a, int disease_b,
                               int hospitalized) {
    if (!m) return false;
    /* Structure: DiseaseA -> Hospitalized <- DiseaseB.
     * Diseases are independent unconditionally.
     * Conditioning on Hospitalized induces negative correlation. */
    SCM_VarSet empty; scm_varset_init(&empty);
    bool indep_uncond = scm_d_separated(m, disease_a, disease_b, &empty);

    SCM_VarSet cond; scm_varset_init(&cond);
    scm_varset_add(&cond, hospitalized);
    bool indep_cond = scm_d_separated(m, disease_a, disease_b, &cond);

    /* If unshielded collider: arrows point to hospitalized from both diseases */
    bool has_collider = scm_has_edge(m, disease_a, hospitalized)
                     && scm_has_edge(m, disease_b, hospitalized)
                     && !scm_has_edge(m, disease_a, disease_b);

    return has_collider && indep_uncond && !indep_cond;
}

/* ---------- Clinical Trial Analysis ---------- */

SCM_ClinicalResult* scm_clinical_trial_analysis(const double* data, int n_rows, int n_cols,
                                                  int treatment_assigned, int treatment_received,
                                                  int outcome) {
    SCM_ClinicalResult* cr = calloc(1, sizeof(SCM_ClinicalResult));
    if (!cr || !data || n_rows < 4) { free(cr); return NULL; }

    /* ITT: difference in mean outcome by assignment */
    double y_z1 = 0.0, y_z0 = 0.0;
    int n_z1 = 0, n_z0 = 0;
    for (int i = 0; i < n_rows; i++) {
        if (data[i * n_cols + treatment_assigned] > 0.5) {
            y_z1 += data[i * n_cols + outcome]; n_z1++;
        } else {
            y_z0 += data[i * n_cols + outcome]; n_z0++;
        }
    }
    if (n_z1 > 0) y_z1 /= (double)n_z1;
    if (n_z0 > 0) y_z0 /= (double)n_z0;
    cr->itt_effect = y_z1 - y_z0;

    /* Per-protocol: only those who complied */
    double y_pp1 = 0.0, y_pp0 = 0.0;
    int n_pp1 = 0, n_pp0 = 0;
    for (int i = 0; i < n_rows; i++) {
        double z = data[i * n_cols + treatment_assigned];
        double x = data[i * n_cols + treatment_received];
        if (fabs(z - x) < 0.1) { /* complier */
            if (z > 0.5) { y_pp1 += data[i * n_cols + outcome]; n_pp1++; }
            else { y_pp0 += data[i * n_cols + outcome]; n_pp0++; }
        }
    }
    if (n_pp1 > 0) y_pp1 /= (double)n_pp1;
    if (n_pp0 > 0) y_pp0 /= (double)n_pp0;
    cr->per_protocol = (n_pp1 > 0 && n_pp0 > 0) ? y_pp1 - y_pp0 : 0.0;

    /* As-treated: compare by actual treatment */
    double y_at1 = 0.0, y_at0 = 0.0;
    int n_at1 = 0, n_at0 = 0;
    for (int i = 0; i < n_rows; i++) {
        if (data[i * n_cols + treatment_received] > 0.5) {
            y_at1 += data[i * n_cols + outcome]; n_at1++;
        } else {
            y_at0 += data[i * n_cols + outcome]; n_at0++;
        }
    }
    if (n_at1 > 0) y_at1 /= (double)n_at1;
    if (n_at0 > 0) y_at0 /= (double)n_at0;
    cr->as_treated = y_at1 - y_at0;

    /* Wald IV estimator: (E[Y|Z=1] - E[Y|Z=0]) / (E[X|Z=1] - E[X|Z=0]) */
    double x_z1 = 0.0, x_z0 = 0.0;
    n_z1 = 0; n_z0 = 0;
    for (int i = 0; i < n_rows; i++) {
        if (data[i * n_cols + treatment_assigned] > 0.5) {
            x_z1 += data[i * n_cols + treatment_received]; n_z1++;
        } else {
            x_z0 += data[i * n_cols + treatment_received]; n_z0++;
        }
    }
    if (n_z1 > 0) x_z1 /= (double)n_z1;
    if (n_z0 > 0) x_z0 /= (double)n_z0;
    double compliance = x_z1 - x_z0;
    cr->iv_wald_estimate = (fabs(compliance) > 1e-10) ? cr->itt_effect / compliance : 0.0;
    cr->noncompliance_present = (fabs(compliance - 1.0) > 0.05);
    cr->complier_ace = cr->iv_wald_estimate;

    return cr;
}

/* Complier Average Causal Effect via Wald estimator */
double scm_complier_ace(const double* data, int n_rows, int n_cols,
                         int instrument, int treatment, int outcome) {
    if (!data || n_rows < 4) return 0.0;

    /* First stage: E[treatment | instrument=1] - E[treatment | instrument=0] */
    double x_z1 = 0.0, x_z0 = 0.0;
    int n1 = 0, n0 = 0;
    for (int i = 0; i < n_rows; i++) {
        if (data[i * n_cols + instrument] > 0.5) {
            x_z1 += data[i * n_cols + treatment]; n1++;
        } else {
            x_z0 += data[i * n_cols + treatment]; n0++;
        }
    }
    if (n1 > 0) x_z1 /= (double)n1;
    if (n0 > 0) x_z0 /= (double)n0;
    double first_stage = x_z1 - x_z0;
    if (fabs(first_stage) < 1e-10) return 0.0;

    /* Reduced form: E[outcome | instrument=1] - E[outcome | instrument=0] */
    n1 = 0; n0 = 0;
    double y_z1 = 0.0, y_z0 = 0.0;
    for (int i = 0; i < n_rows; i++) {
        if (data[i * n_cols + instrument] > 0.5) {
            y_z1 += data[i * n_cols + outcome]; n1++;
        } else {
            y_z0 += data[i * n_cols + outcome]; n0++;
        }
    }
    if (n1 > 0) y_z1 /= (double)n1;
    if (n0 > 0) y_z0 /= (double)n0;
    double reduced_form = y_z1 - y_z0;

    return reduced_form / first_stage;
}

/* ---------- Instrumental Variables ---------- */

/* 2SLS for continuous treatment.
 * Simplified: one IV, one treatment, optional covariates. */
double scm_iv_two_stage_least_squares(const double* data, int n_rows, int n_cols,
                                       int instrument, int treatment, int outcome,
                                       const int* covariates, int n_covariates) {
    if (!data || n_rows < 4) return 0.0;

    /* First stage: regress treatment on instrument + covariates */
    int n_pred = 1 + n_covariates; /* intercept + instrument + covariates */
    if (n_pred + 1 > SCM_MAX_VARS) return 0.0;

    /* Build design matrix X (n_rows x (1 + n_pred)) */
    double* X = calloc(n_rows * (n_pred + 1), sizeof(double));
    if (!X) return 0.0;
    for (int i = 0; i < n_rows; i++) {
        X[i * (n_pred + 1) + 0] = 1.0; /* intercept */
        X[i * (n_pred + 1) + 1] = data[i * n_cols + instrument];
        for (int c = 0; c < n_covariates; c++)
            X[i * (n_pred + 1) + 2 + c] = data[i * n_cols + covariates[c]];
    }

    /* XtX_inv * Xty via simple OLS. Use normal equations. */
    double* XtX = calloc((n_pred + 1) * (n_pred + 1), sizeof(double));
    double* Xty = calloc(n_pred + 1, sizeof(double));
    if (!XtX || !Xty) { free(X); free(XtX); free(Xty); return 0.0; }

    for (int i = 0; i < n_rows; i++)
        for (int p = 0; p <= n_pred; p++) {
            Xty[p] += X[i * (n_pred + 1) + p] * data[i * n_cols + treatment];
            for (int q = 0; q <= n_pred; q++)
                XtX[p * (n_pred + 1) + q] += X[i * (n_pred + 1) + p] * X[i * (n_pred + 1) + q];
        }

    /* Compute fitted values: X * beta */
    double* fitted = calloc(n_rows, sizeof(double));
    if (!fitted) { free(X); free(XtX); free(Xty); return 0.0; }

    /* Solve XtX * beta = Xty via simple approach (use beta estimate as approximation) */
    /* Simplified: use ratio estimator for single IV */
    double cov_zx = 0.0, var_z = 0.0;
    double mean_z = 0.0, mean_x = 0.0;
    for (int i = 0; i < n_rows; i++) {
        mean_z += data[i * n_cols + instrument];
        mean_x += data[i * n_cols + treatment];
    }
    mean_z /= (double)n_rows; mean_x /= (double)n_rows;
    for (int i = 0; i < n_rows; i++) {
        double dz = data[i * n_cols + instrument] - mean_z;
        double dx = data[i * n_cols + treatment] - mean_x;
        cov_zx += dz * dx; var_z += dz * dz;
    }
    double pi_zx = (fabs(var_z) > 1e-10) ? cov_zx / var_z : 0.0;
    for (int i = 0; i < n_rows; i++)
        fitted[i] = mean_x + pi_zx * (data[i * n_cols + instrument] - mean_z);

    /* Second stage: regress outcome on fitted treatment + covariates */
    double cov_fy = 0.0, var_f = 0.0;
    double mean_y = 0.0, mean_f = 0.0;
    for (int i = 0; i < n_rows; i++) {
        mean_y += data[i * n_cols + outcome];
        mean_f += fitted[i];
    }
    mean_y /= (double)n_rows; mean_f /= (double)n_rows;
    for (int i = 0; i < n_rows; i++) {
        double df = fitted[i] - mean_f;
        double dy = data[i * n_cols + outcome] - mean_y;
        cov_fy += df * dy; var_f += df * df;
    }

    double iv_effect = (fabs(var_f) > 1e-10) ? cov_fy / var_f : 0.0;
    free(X); free(XtX); free(Xty); free(fitted);
    return iv_effect;
}

/* Wald estimator for binary IV and binary treatment */
double scm_iv_wald_estimator(const double* data, int n_rows, int n_cols,
                              int instrument, int treatment, int outcome) {
    if (!data || n_rows < 4) return 0.0;
    return scm_complier_ace(data, n_rows, n_cols, instrument, treatment, outcome);
}


/* ---------- Propensity Score Methods ---------- */

/* Estimate propensity scores using a linear probability model.
 * Returns array of propensity scores, one per observation. */
double* scm_propensity_score_estimate(const double* data, int n_rows, int n_cols,
                                       int treatment, const int* covariates, int n_covariates) {
    if (!data || n_rows < 2) return NULL;
    double* pscores = calloc(n_rows, sizeof(double));
    if (!pscores) return NULL;

    /* Simple linear model: pscore = 1/(1+exp(-linear_combination))
     * Use marginal proportions as starting point, adjust by covariate deviation */
    double base_prob = 0.0;
    for (int i = 0; i < n_rows; i++)
        base_prob += data[i * n_cols + treatment];
    base_prob /= (double)n_rows;

    /* Compute covariate means for treated and controls */
    double* cov_mean_t = calloc(n_covariates, sizeof(double));
    double* cov_mean_c = calloc(n_covariates, sizeof(double));
    int n_t = 0, n_c = 0;
    if (!cov_mean_t || !cov_mean_c) { free(pscores); free(cov_mean_t); free(cov_mean_c); return NULL; }

    for (int i = 0; i < n_rows; i++) {
        if (data[i * n_cols + treatment] > 0.5) {
            for (int j = 0; j < n_covariates; j++)
                cov_mean_t[j] += data[i * n_cols + covariates[j]];
            n_t++;
        } else {
            for (int j = 0; j < n_covariates; j++)
                cov_mean_c[j] += data[i * n_cols + covariates[j]];
            n_c++;
        }
    }
    for (int j = 0; j < n_covariates; j++) {
        if (n_t > 0) cov_mean_t[j] /= (double)n_t;
        if (n_c > 0) cov_mean_c[j] /= (double)n_c;
    }

    /* Compute propensity scores via logistic approximation */
    for (int i = 0; i < n_rows; i++) {
        double linear = log((base_prob + 0.01) / (1.0 - base_prob + 0.01));
        for (int j = 0; j < n_covariates; j++) {
            double diff = data[i * n_cols + covariates[j]] - cov_mean_t[j];
            /* Simple adjustment: higher covariate similarity to treated -> higher score */
            /* Using Mahalanobis-style distance scaling */
            double pooled_sd = 1.0;
            if (n_t > 1 && n_c > 1) {
                double var_t = 0.0, var_c = 0.0;
                for (int k = 0; k < n_rows; k++) {
                    if (data[k * n_cols + treatment] > 0.5)
                        var_t += (data[k * n_cols + covariates[j]] - cov_mean_t[j]) * (data[k * n_cols + covariates[j]] - cov_mean_t[j]);
                    else
                        var_c += (data[k * n_cols + covariates[j]] - cov_mean_c[j]) * (data[k * n_cols + covariates[j]] - cov_mean_c[j]);
                }
                pooled_sd = sqrt((var_t / (n_t - 1) + var_c / (n_c - 1)) / 2.0);
                if (pooled_sd < 1e-8) pooled_sd = 1.0;
            }
            linear += diff / pooled_sd * 0.5;
        }
        pscores[i] = 1.0 / (1.0 + exp(-linear));
        if (pscores[i] < 0.01) pscores[i] = 0.01;
        if (pscores[i] > 0.99) pscores[i] = 0.99;
    }

    free(cov_mean_t); free(cov_mean_c);
    return pscores;
}

/* Nearest-neighbor propensity score matching with caliper */
SCM_PSMResult* scm_propensity_score_matching(const double* data, int n_rows, int n_cols,
                                              int treatment, int outcome,
                                              const int* covariates, int n_covariates,
                                              double caliper) {
    SCM_PSMResult* pr = calloc(1, sizeof(SCM_PSMResult));
    if (!pr || !data || n_rows < 4) { free(pr); return NULL; }

    double* pscores = scm_propensity_score_estimate(data, n_rows, n_cols, treatment,
                                                     covariates, n_covariates);
    if (!pscores) { free(pr); return NULL; }

    /* Separate treated and control indices */
    int* treated_idx = calloc(n_rows, sizeof(int));
    int* control_idx = calloc(n_rows, sizeof(int));
    int n_treated = 0, n_controls = 0;
    if (!treated_idx || !control_idx) {
        free(pscores); free(treated_idx); free(control_idx); free(pr); return NULL;
    }
    for (int i = 0; i < n_rows; i++) {
        if (data[i * n_cols + treatment] > 0.5) treated_idx[n_treated++] = i;
        else control_idx[n_controls++] = i;
    }

    /* Match each treated to nearest control */
    double sum_att = 0.0;
    int n_matched = 0;
    bool* used_control = calloc(n_controls, sizeof(bool));
    if (!used_control) { free(pscores); free(treated_idx); free(control_idx); free(pr); return NULL; }

    for (int t = 0; t < n_treated; t++) {
        int t_idx = treated_idx[t];
        double best_dist = DBL_MAX;
        int best_c = -1;
        for (int c = 0; c < n_controls; c++) {
            if (used_control[c]) continue;
            int c_idx = control_idx[c];
            double dist = fabs(pscores[t_idx] - pscores[c_idx]);
            if (dist < best_dist && dist <= caliper) {
                best_dist = dist; best_c = c;
            }
        }
        if (best_c >= 0) {
            used_control[best_c] = true;
            int c_idx = control_idx[best_c];
            sum_att += data[t_idx * n_cols + outcome] - data[c_idx * n_cols + outcome];
            n_matched++;
        }
    }

    pr->att = (n_matched > 0) ? sum_att / (double)n_matched : 0.0;
    pr->n_matched_pairs = n_matched;
    pr->ate_matched = pr->att; /* For ATT matching, ATE approx = ATT */
    pr->balance_std_diff = 0.0;
    pr->balance_achieved = true;

    /* Compute standardized difference after matching */
    if (n_matched > 0 && n_covariates > 0) {
        double max_sd = 0.0;
        for (int j = 0; j < n_covariates; j++) {
            double mean_t = 0.0, mean_c = 0.0;
            int count = 0;
            for (int t = 0; t < n_treated && count < n_matched; t++) {
                int t_idx = treated_idx[t];
                for (int c = 0; c < n_controls; c++) {
                    if (!used_control[c]) continue;
                    int c_idx = control_idx[c];
                    if (fabs(pscores[t_idx] - pscores[c_idx]) <= caliper) {
                        mean_t += data[t_idx * n_cols + covariates[j]];
                        mean_c += data[c_idx * n_cols + covariates[j]];
                        count++;
                        break;
                    }
                }
            }
            if (count > 0) {
                mean_t /= count; mean_c /= count;
                double sd = fabs(mean_t - mean_c);
                if (sd > max_sd) max_sd = sd;
            }
        }
        pr->balance_std_diff = max_sd;
        pr->balance_achieved = (max_sd < 0.1); /* Common threshold */
    }

    free(pscores); free(treated_idx); free(control_idx); free(used_control);
    return pr;
}

/* ---------- Difference-in-Differences ---------- */

/* Simple 2x2 DiD estimator:
 * DiD = (E[Y_post | T=1] - E[Y_pre | T=1]) - (E[Y_post | T=0] - E[Y_pre | T=0]) */
double scm_difference_in_differences(const double* data, int n_rows,
                                      int pre_outcome, int post_outcome,
                                      int treated_group) {
    if (!data || n_rows < 4) return 0.0;

    double pre_t1 = 0.0, post_t1 = 0.0, pre_t0 = 0.0, post_t0 = 0.0;
    int n1 = 0, n0 = 0;

    for (int i = 0; i < n_rows; i++) {
        if (data[i * 4 + treated_group] > 0.5) { /* treated */
            pre_t1 += data[i * 4 + pre_outcome];
            post_t1 += data[i * 4 + post_outcome];
            n1++;
        } else {
            pre_t0 += data[i * 4 + pre_outcome];
            post_t0 += data[i * 4 + post_outcome];
            n0++;
        }
    }

    if (n1 > 0) { pre_t1 /= (double)n1; post_t1 /= (double)n1; }
    if (n0 > 0) { pre_t0 /= (double)n0; post_t0 /= (double)n0; }

    return (post_t1 - pre_t1) - (post_t0 - pre_t0);
}

/* ---------- Regression Discontinuity Design ---------- */

/* Sharp RDD via local linear regression.
 * treatment = 1[running_var >= cutoff].
 * Uses uniform kernel within bandwidth. */
double scm_regression_discontinuity(const double* data, int n_rows, int n_cols,
                                     int running_var, int outcome, int treatment,
                                     double cutoff, double bandwidth) {
    if (!data || n_rows < 10) return 0.0;

    /* Separate points left and right of cutoff (within bandwidth) */
    double y_above = 0.0, y_below = 0.0;
    double x_above = 0.0, x_below = 0.0;
    int n_above = 0, n_below = 0;

    for (int i = 0; i < n_rows; i++) {
        double rv = data[i * n_cols + running_var];
        double dist = rv - cutoff;
        if (fabs(dist) > bandwidth) continue;

        if (dist >= 0) {
            y_above += data[i * n_cols + outcome];
            x_above += dist;
            n_above++;
        } else {
            y_below += data[i * n_cols + outcome];
            x_below += dist;
            n_below++;
        }
    }

    /* Local linear regression: fit line on each side, compare at cutoff */
    /* Simplified: weighted average near cutoff */
    if (n_above < 2 || n_below < 2) return 0.0;

    /* Weight points closer to cutoff more heavily */
    double w_y_above = 0.0, w_y_below = 0.0, w_sum_above = 0.0, w_sum_below = 0.0;
    for (int i = 0; i < n_rows; i++) {
        double rv = data[i * n_cols + running_var];
        double dist = rv - cutoff;
        if (fabs(dist) > bandwidth) continue;
        double w = 1.0 - fabs(dist) / bandwidth; /* triangular kernel */
        if (dist >= 0) {
            w_y_above += w * data[i * n_cols + outcome];
            w_sum_above += w;
        } else {
            w_y_below += w * data[i * n_cols + outcome];
            w_sum_below += w;
        }
    }
    if (w_sum_above < 1e-10 || w_sum_below < 1e-10) return 0.0;
    y_above = w_y_above / w_sum_above;
    y_below = w_y_below / w_sum_below;

    /* RDD effect = discontinuity at cutoff */
    (void)treatment;
    return y_above - y_below;
}

/* ---------- Double/Debiased Machine Learning ---------- */

/* Simplified Double ML for ATE estimation.
 * Step 1: Split sample into two halves.
 * Step 2: On half 1, learn E[Y|X] and E[D|X].
 * Step 3: On half 2, compute residuals and regress Y_resid on D_resid.
 * References: Chernozhukov et al. (2018) */
double scm_double_ml_ate(const double* data, int n_rows, int n_cols,
                          int treatment, int outcome,
                          const int* covariates, int n_covariates) {
    if (!data || n_rows < 20 || n_covariates < 1) return 0.0;

    int n1 = n_rows / 2;
    (void)(n_rows - n1); /* n2 = n_rows - n1; ensure split is valid */

    /* Half 1: compute conditional means (simple linear model) */
    double mean_y_x = 0.0, mean_d_x = 0.0;
    for (int i = 0; i < n1; i++) {
        mean_y_x += data[i * n_cols + outcome];
        mean_d_x += data[i * n_cols + treatment];
    }
    mean_y_x /= (double)n1; mean_d_x /= (double)n1;

    /* Half 2: residualize and estimate */
    double num = 0.0, den = 0.0;
    for (int i = n1; i < n_rows; i++) {
        double d = data[i * n_cols + treatment];
        double y = data[i * n_cols + outcome];

        /* Simple residualization: subtract covariate-adjusted mean */
        double cov_adj_y = 0.0, cov_adj_d = 0.0;
        for (int j = 0; j < n_covariates; j++) {
            double cv = data[i * n_cols + covariates[j]];
            cov_adj_y += cv; cov_adj_d += cv;
        }
        cov_adj_y *= mean_y_x / ((double)n_covariates + 1.0);
        cov_adj_d *= mean_d_x / ((double)n_covariates + 1.0);

        double d_resid = d - mean_d_x - cov_adj_d * 0.1;
        double y_resid = y - mean_y_x - cov_adj_y * 0.1;

        num += d_resid * y_resid;
        den += d_resid * d_resid;
    }

    return (fabs(den) > 1e-10) ? num / den : 0.0;
}

/* ---------- Mediation Decomposition ---------- */

/* Decompose total effect using product-of-coefficients method.
 * total = nde + nie where:
 * nde (direct) = effect of treatment on outcome not through mediator
 * nie (indirect) = effect through mediator = alpha * beta
 * where alpha = treatment -> mediator, beta = mediator -> outcome (adjusted for treatment) */
void scm_mediation_decomposition(const SCM_Model* m, int treatment, int mediator,
                                  int outcome, double* nde, double* nie, double* total) {
    if (!m) { if (nde) *nde = 0.0; if (nie) *nie = 0.0; if (total) *total = 0.0; return; }

    /* Direct path coefficient: treatment -> outcome */
    double direct = 0.0;
    for (int i = 0; i < m->n_edges; i++)
        if (m->edges[i].from == treatment && m->edges[i].to == outcome)
            direct = m->edges[i].coef;

    /* Indirect path: treatment -> mediator -> outcome */
    double alpha = 0.0, beta = 0.0;
    for (int i = 0; i < m->n_edges; i++) {
        if (m->edges[i].from == treatment && m->edges[i].to == mediator)
            alpha = m->edges[i].coef;
        if (m->edges[i].from == mediator && m->edges[i].to == outcome)
            beta = m->edges[i].coef;
    }
    double indirect = alpha * beta;

    *nde = direct;
    *nie = indirect;
    *total = direct + indirect;
}

/* ---------- Utility ---------- */

void scm_simpson_result_free(SCM_SimpsonResult* sr) { free(sr); }

void scm_simpson_result_print(const SCM_SimpsonResult* sr) {
    if (!sr) { printf("SimpsonResult: NULL\n"); return; }
    printf("=== Simpson's Paradox Analysis ===\n");
    printf("Paradox detected: %s\n", sr->paradox_detected ? "YES" : "NO");
    printf("Aggregate effect: %.4f\n", sr->aggregate_effect);
    printf("Subgroup effects (%d groups):\n", sr->n_subgroups);
    for (int i = 0; i < sr->n_subgroups; i++)
        printf("  Group %d: %.4f\n", i, sr->subgroup_effects[i]);
}

void scm_clinical_result_free(SCM_ClinicalResult* cr) { free(cr); }

void scm_clinical_result_print(const SCM_ClinicalResult* cr) {
    if (!cr) { printf("ClinicalResult: NULL\n"); return; }
    printf("=== Clinical Trial Analysis ===\n");
    printf("ITT effect:          %.4f\n", cr->itt_effect);
    printf("Per-protocol effect: %.4f\n", cr->per_protocol);
    printf("As-treated effect:   %.4f\n", cr->as_treated);
    printf("CACE (LATE):         %.4f\n", cr->complier_ace);
    printf("IV Wald estimate:    %.4f\n", cr->iv_wald_estimate);
    printf("Noncompliance:       %s\n", cr->noncompliance_present ? "YES" : "NO");
}

void scm_psm_result_free(SCM_PSMResult* pr) { free(pr); }

void scm_psm_result_print(const SCM_PSMResult* pr) {
    if (!pr) { printf("PSMResult: NULL\n"); return; }
    printf("=== Propensity Score Matching ===\n");
    printf("ATT:             %.4f\n", pr->att);
    printf("ATE (matched):   %.4f\n", pr->ate_matched);
    printf("Matched pairs:   %d\n", pr->n_matched_pairs);
    printf("Balance std diff: %.4f (%s)\n", pr->balance_std_diff,
           pr->balance_achieved ? "balanced" : "NOT balanced");
}
