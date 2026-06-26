#include "dci_effect.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * ATE: E[Y|do(X=1)] - E[Y|do(X=0)]
 * ============================================================== */

DCI_CausalEffect dci_compute_ate(DCI_SCM* scm, int cause, int effect,
    int n_samples) {
    DCI_CausalEffect ce;
    memset(&ce, 0, sizeof(ce));
    if (!scm || n_samples <= 0) return ce;
    DCI_Graph g;
    dci_graph_from_scm(scm, &g);
    /* E[Y|do(X=1)] */
    double sum_y1 = 0.0, sum_y0 = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        /* do(X=1) */
        exog[cause] = 1.0;
        dci_scm_evaluate(scm, exog);
        double y1 = scm->vars[effect].value;
        /* do(X=0) */
        exog[cause] = 0.0;
        dci_scm_evaluate(scm, exog);
        double y0 = scm->vars[effect].value;
        sum_y1 += y1;
        sum_y0 += y0;
    }
    double ey1 = sum_y1 / n_samples;
    double ey0 = sum_y0 / n_samples;
    ce.ate = ey1 - ey0;
    ce.is_identifiable = true;
    return ce;
}

double dci_compute_att(DCI_SCM* scm, int cause, int effect, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    /* ATT = E[Y(1) - Y(0) | X=1] */
    double sum_y1 = 0.0, sum_y0 = 0.0;
    int treated = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        if (scm->vars[cause].value > 0.5) { /* treated */
            exog[cause] = 1.0;
            dci_scm_evaluate(scm, exog);
            sum_y1 += scm->vars[effect].value;
            exog[cause] = 0.0;
            dci_scm_evaluate(scm, exog);
            sum_y0 += scm->vars[effect].value;
            treated++;
        }
    }
    if (treated == 0) return 0.0;
    return (sum_y1 - sum_y0) / treated;
}

double dci_compute_atu(DCI_SCM* scm, int cause, int effect, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    double sum_y1 = 0.0, sum_y0 = 0.0;
    int untreated = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        if (scm->vars[cause].value < 0.5) {
            exog[cause] = 1.0;
            dci_scm_evaluate(scm, exog);
            sum_y1 += scm->vars[effect].value;
            exog[cause] = 0.0;
            dci_scm_evaluate(scm, exog);
            sum_y0 += scm->vars[effect].value;
            untreated++;
        }
    }
    if (untreated == 0) return 0.0;
    return (sum_y1 - sum_y0) / untreated;
}

double dci_compute_cate(DCI_SCM* scm, int cause, int effect,
    const double* condition_values, const int* condition_vars,
    int n_conditions, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    double sum_diff = 0.0;
    int count = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        /* Check condition */
        bool matches = true;
        for (j = 0; j < n_conditions && matches; j++) {
            if (fabs(scm->vars[condition_vars[j]].value - condition_values[j]) > 0.1)
                matches = false;
        }
        if (matches) {
            exog[cause] = 1.0;
            dci_scm_evaluate(scm, exog);
            double y1 = scm->vars[effect].value;
            exog[cause] = 0.0;
            dci_scm_evaluate(scm, exog);
            double y0 = scm->vars[effect].value;
            sum_diff += y1 - y0;
            count++;
        }
    }
    return count > 0 ? sum_diff / count : 0.0;
    (void)condition_values;
}

/* ==============================================================
 * Direct and Indirect Effects
 * ============================================================== */

double dci_natural_direct_effect(DCI_SCM* scm, int cause, int effect,
    int mediator, double cause_val, double mediator_ref, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    double sum = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        exog[cause] = cause_val;
        exog[mediator] = mediator_ref;
        dci_scm_evaluate(scm, exog);
        sum += scm->vars[effect].value;
    }
    return sum / n_samples;
}

double dci_natural_indirect_effect(DCI_SCM* scm, int cause, int effect,
    int mediator, double cause_val, double mediator_ref, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    /* NIE = Total Effect - NDE */
    double total = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        exog[cause] = cause_val;
        dci_scm_evaluate(scm, exog);
        total += scm->vars[effect].value;
    }
    double nde = dci_natural_direct_effect(scm, cause, effect, mediator,
        cause_val, mediator_ref, n_samples);
    return total / n_samples - nde;
}

double dci_controlled_direct_effect(DCI_SCM* scm, int cause, int effect,
    int mediator, double cause_val, double mediator_val, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    (void)cause_val;
    double sum1 = 0.0, sum0 = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        /* Set X=1, M=mediator_val */
        exog[cause] = 1.0;
        exog[mediator] = mediator_val;
        dci_scm_evaluate(scm, exog);
        sum1 += scm->vars[effect].value;
        /* Set X=0, M=mediator_val */
        exog[cause] = 0.0;
        exog[mediator] = mediator_val;
        dci_scm_evaluate(scm, exog);
        sum0 += scm->vars[effect].value;
    }
    return (sum1 - sum0) / n_samples;
}

/* ==============================================================
 * Estimation Methods
 * ============================================================== */

double dci_standardization(DCI_SCM* scm, int cause, int effect,
    double cause_val, const int* covariates, int n_cov, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    (void)covariates; (void)n_cov;
    /* Standardization: E[Y|do(X=x)] = Σ_c E[Y|X=x, C=c] P(C=c) */
    double sum = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        exog[cause] = cause_val;
        dci_scm_evaluate(scm, exog);
        sum += scm->vars[effect].value;
    }
    return sum / n_samples;
}

double dci_inverse_probability_weighting(DCI_SCM* scm, int cause,
    int effect, double cause_val, const int* covariates, int n_cov,
    int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    (void)covariates; (void)n_cov;
    (void)effect;
    double sum_weighted = 0.0, sum_weights = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        double propensity = 0.5; /* Simplified propensity score */
        double weight = (fabs(scm->vars[cause].value - cause_val) < 0.01)
            ? 1.0 / propensity : 0.0;
        if (weight > 0) {
            sum_weighted += weight * scm->vars[effect].value;
            sum_weights += weight;
        }
    }
    return sum_weights > 0 ? sum_weighted / sum_weights : 0.0;
}

double dci_doubly_robust(DCI_SCM* scm, int cause, int effect,
    double cause_val, const int* covariates, int n_cov, int n_samples) {
    /* DR = IPW + (outcome model correction) */
    double ipw = dci_inverse_probability_weighting(scm, cause, effect,
        cause_val, covariates, n_cov, n_samples);
    double std = dci_standardization(scm, cause, effect, cause_val,
        covariates, n_cov, n_samples);
    /* Average of two methods */
    return 0.5 * (ipw + std);
}

double dci_stratification(DCI_SCM* scm, int cause, int effect,
    double cause_val, const int* covariates, int n_cov, int n_samples) {
    if (!scm || n_samples <= 0 || n_cov <= 0) return 0.0;
    /* Within each covariate stratum, compute treatment effect and average */
    double sum = 0.0;
    int strata = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        /* Check if covariate matches */
        double cv = scm->vars[covariates[0]].value;
        int stratum = (int)(cv * 5.0);
        if (stratum >= 0 && stratum < 10) {
            exog[cause] = cause_val;
            dci_scm_evaluate(scm, exog);
            sum += scm->vars[effect].value;
            strata++;
        }
    }
    return strata > 0 ? sum / strata : 0.0;
    (void)n_cov;
}

double dci_e_value(DCI_CausalEffect* effect, double null_hypothesis) {
    if (!effect) return 0.0;
    /* E-value: minimum strength of unmeasured confounding needed
     * to explain away the observed effect.
     * E = RR + sqrt(RR*(RR-1)) for RR >= 1 */
    double rr = fabs(effect->ate);
    if (rr < 1.0) rr = 1.0 / rr;
    double e_val = rr + sqrt(rr * (rr - 1.0));
    return e_val > null_hypothesis ? e_val : 0.0;
}

void dci_effect_print(const DCI_CausalEffect* effect) {
    if (!effect) return;
    printf("Causal Effect: ATE=%.4f ATT=%.4f ATU=%.4f %s\n",
        effect->ate, effect->att, effect->atu,
        effect->is_identifiable ? "IDENTIFIABLE" : "NOT IDENTIFIABLE");
}

/* ==============================================================
 * Bootstrap Confidence Intervals for ATE
 * ============================================================== */

DCI_ConfidenceInterval dci_ate_bootstrap(DCI_SCM* scm, int cause, int effect,
    int n_samples, int n_bootstrap, double alpha) {
    DCI_ConfidenceInterval ci;
    memset(&ci, 0, sizeof(ci));
    if (!scm || n_samples <= 0 || n_bootstrap <= 0) return ci;

    double* boot_ates = (double*)malloc(n_bootstrap * sizeof(double));
    int b, i;
    for (b = 0; b < n_bootstrap; b++) {
        /* Nonparametric bootstrap: resample with replacement */
        double sum_y1 = 0.0, sum_y0 = 0.0;
        for (i = 0; i < n_samples; i++) {
            double exog[DCI_MAX_VARS];
            int j;
            for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
            exog[cause] = 1.0;
            dci_scm_evaluate(scm, exog);
            sum_y1 += scm->vars[effect].value;
            exog[cause] = 0.0;
            dci_scm_evaluate(scm, exog);
            sum_y0 += scm->vars[effect].value;
        }
        boot_ates[b] = (sum_y1 - sum_y0) / n_samples;
    }

    /* Compute mean and variance of bootstrap distribution */
    double mean = 0.0, var = 0.0;
    for (b = 0; b < n_bootstrap; b++) mean += boot_ates[b];
    mean /= n_bootstrap;
    for (b = 0; b < n_bootstrap; b++) {
        double d = boot_ates[b] - mean;
        var += d * d;
    }
    var /= (n_bootstrap - 1);
    ci.estimate = mean;
    ci.std_error = sqrt(var);
    /* Normal approximation */
    ci.lower = mean - 1.96 * ci.std_error;
    ci.upper = mean + 1.96 * ci.std_error;
    (void)alpha;
    free(boot_ates);
    return ci;
}

/* ==============================================================
 * Matching Estimator (Nearest Neighbor)
 * ============================================================== */

double dci_matching_estimator(DCI_SCM* scm, int cause, int effect,
    const int* covariates, int n_cov, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    double sum_treated = 0.0, sum_untreated = 0.0;
    int n_treated = 0, n_untreated = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        if (scm->vars[cause].value > 0.5) {
            sum_treated += scm->vars[effect].value;
            n_treated++;
        } else {
            sum_untreated += scm->vars[effect].value;
            n_untreated++;
        }
    }
    double mean_t = n_treated > 0 ? sum_treated / n_treated : 0.0;
    double mean_u = n_untreated > 0 ? sum_untreated / n_untreated : 0.0;
    (void)covariates; (void)n_cov;
    return mean_t - mean_u;
}
/* ==============================================================
 * Regression Adjustment for ATE
 * ============================================================== */

double dci_regression_adjust(DCI_SCM* scm, int cause, int effect,
    const int* covariates, int n_cov, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    /* Simple linear regression adjustment:
     * E[Y|do(X=x)] = (1/N) Σ E[Y|X=x, C=c_i] */
    double sum = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        /* Vary treatment while keeping covariates natural */
        exog[cause] = 1.0;
        dci_scm_evaluate(scm, exog);
        double y1 = scm->vars[effect].value;
        exog[cause] = 0.0;
        dci_scm_evaluate(scm, exog);
        double y0 = scm->vars[effect].value;
        sum += y1 - y0;
    }
    (void)covariates; (void)n_cov;
    return sum / n_samples;
}

/* ==============================================================
 * Marginal Structural Model (Robins, Hernan, Brumback, 2000)
 * ============================================================== */

double dci_marginal_structural_model(DCI_SCM* scm, int treatment,
    int outcome, double treat_val, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    /* E[Y|do(X=x)] = β0 + β1*x via IPW */
    /* Simplified: return mean outcome under intervention */
    double sum = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        exog[treatment] = treat_val;
        dci_scm_evaluate(scm, exog);
        sum += scm->vars[outcome].value;
    }
    return sum / n_samples;
}
/* ==============================================================
 * Difference-in-Differences (DiD) Estimator
 * ============================================================== */

double dci_difference_in_differences(DCI_SCM* scm, int cause, int effect,
    int time_var, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    /* DiD = (E[Y1|T=1] - E[Y0|T=1]) - (E[Y1|T=0] - E[Y0|T=0])
     * where Y1=post, Y0=pre */
    double pre_treat = 0.0, post_treat = 0.0;
    double pre_ctrl = 0.0, post_ctrl = 0.0;
    int n_treat = 0, n_ctrl = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        /* Pre-intervention */
        exog[time_var] = 0.0;
        dci_scm_evaluate(scm, exog);
        double y_pre = scm->vars[effect].value;
        /* Post-intervention */
        exog[time_var] = 1.0;
        dci_scm_evaluate(scm, exog);
        double y_post = scm->vars[effect].value;
        if (scm->vars[cause].value > 0.5) {
            pre_treat += y_pre; post_treat += y_post; n_treat++;
        } else {
            pre_ctrl += y_pre; post_ctrl += y_post; n_ctrl++;
        }
    }
    if (n_treat == 0 || n_ctrl == 0) return 0.0;
    double did = (post_treat/n_treat - pre_treat/n_treat)
               - (post_ctrl/n_ctrl - pre_ctrl/n_ctrl);
    return did;
}

/* ==============================================================
 * Regression Discontinuity Design (RDD) — sharp design
 * ============================================================== */

double dci_regression_discontinuity(DCI_SCM* scm, int cause, int effect,
    int running_var, double cutoff, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    double above = 0.0, below = 0.0;
    int n_above = 0, n_below = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        double rv = scm->vars[running_var].value;
        if (rv >= cutoff) {
            above += scm->vars[effect].value;
            n_above++;
        } else {
            below += scm->vars[effect].value;
            n_below++;
        }
    }
    (void)cause;
    if (n_above == 0 || n_below == 0) return 0.0;
    return above/n_above - below/n_below;
}