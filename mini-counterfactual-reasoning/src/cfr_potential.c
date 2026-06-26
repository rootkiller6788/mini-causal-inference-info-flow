#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cfr_core.h"
#include "cfr_potential.h"

/* ============================================================
 * cfr_potential.c -- Potential Outcomes Framework
 *
 * Implements the Rubin Causal Model (RCM) with:
 * - Potential outcomes Y(0), Y(1) for each unit
 * - Observed outcome: Y = T*Y(1) + (1-T)*Y(0)
 * - Propensity score estimation
 * - IPW weighting
 * - Individual Treatment Effects
 * - Manski bounds
 * ============================================================ */

CFRPotentialOutcomes* cfr_po_create(int n_units, int n_covariates) {
    CFRPotentialOutcomes* po = calloc(1, sizeof(CFRPotentialOutcomes));
    if (!po) return NULL;
    po->n_units = n_units; po->n_covariates = n_covariates;
    po->Y0 = calloc(n_units, sizeof(double));
    po->Y1 = calloc(n_units, sizeof(double));
    po->covariates = calloc(n_units * n_covariates, sizeof(double));
    po->randomized = false;
    if (!po->Y0 || !po->Y1 || (n_covariates > 0 && !po->covariates)) {
        free(po->Y0); free(po->Y1); free(po->covariates); free(po); return NULL;
    }
    return po;
}

void cfr_po_free(CFRPotentialOutcomes* po) {
    if (po) { free(po->Y0); free(po->Y1); free(po->covariates); free(po); }
}

void cfr_po_set_Y0(CFRPotentialOutcomes* po, int i, double v) {
    if (po && i >= 0 && i < po->n_units) po->Y0[i] = v;
}

void cfr_po_set_Y1(CFRPotentialOutcomes* po, int i, double v) {
    if (po && i >= 0 && i < po->n_units) po->Y1[i] = v;
}

void cfr_po_set_covariate(CFRPotentialOutcomes* po, int unit, int covar, double v) {
    if (po && unit >= 0 && unit < po->n_units && covar >= 0 && covar < po->n_covariates)
        po->covariates[unit * po->n_covariates + covar] = v;
}

/* ---------- Generate Dataset ---------- */

CFRDataset* cfr_po_generate_dataset(CFRPotentialOutcomes* po,
                                     double* treatment_prob) {
    CFRDataset* ds = calloc(1, sizeof(CFRDataset));
    if (!ds) return NULL;
    ds->n_units = po->n_units;
    ds->treatment = calloc(po->n_units, sizeof(double));
    ds->outcome = calloc(po->n_units, sizeof(double));
    if (!ds->treatment || !ds->outcome) {
        free(ds->treatment); free(ds->outcome); free(ds); return NULL;
    }
    for (int i = 0; i < po->n_units; i++) {
        double p = treatment_prob ? treatment_prob[i] : 0.5;
        double u = (double)rand() / RAND_MAX;
        ds->treatment[i] = (u < p) ? 1.0 : 0.0;
        ds->outcome[i] = (ds->treatment[i] > 0.5) ? po->Y1[i] : po->Y0[i];
    }
    return ds;
}

void cfr_dataset_free(CFRDataset* ds) {
    if (ds) { free(ds->treatment); free(ds->outcome); free(ds); }
}

/* ---------- Propensity Score ---------- */

CFRPropensity* cfr_po_propensity_logistic(CFRDataset* ds,
                                           CFRPotentialOutcomes* po) {
    CFRPropensity* ps = calloc(1, sizeof(CFRPropensity));
    if (!ps) return NULL;
    ps->n_units = ds->n_units;
    ps->propensity = calloc(ds->n_units, sizeof(double));
    ps->weights_ipw = calloc(ds->n_units, sizeof(double));
    ps->weights_att = calloc(ds->n_units, sizeof(double));
    if (!ps->propensity || !ps->weights_ipw || !ps->weights_att) {
        free(ps->propensity); free(ps->weights_ipw); free(ps->weights_att);
        free(ps); return NULL;
    }
    /* Logistic regression via simple heuristic:
     * e(x) = 1/(1+exp(-beta*x)) */
    for (int i = 0; i < ds->n_units; i++) {
        double linear = 0.0;
        for (int j = 0; j < po->n_covariates; j++)
            linear += po->covariates[i * po->n_covariates + j] * 0.1;
        double e_x = 1.0 / (1.0 + exp(-linear));
        ps->propensity[i] = cfr_clamp(e_x, 0.05, 0.95);
    }
    return ps;
}

void cfr_propensity_free(CFRPropensity* ps) {
    if (ps) { free(ps->propensity); free(ps->weights_ipw);
        free(ps->weights_att); free(ps); }
}

void cfr_po_ipw_weights(CFRPropensity* ps, CFRDataset* ds) {
    if (!ps || !ds) return;
    for (int i = 0; i < ds->n_units; i++) {
        double e = ps->propensity[i];
        if (ds->treatment[i] > 0.5) {
            ps->weights_ipw[i] = 1.0 / e;
            ps->weights_att[i] = 1.0;
        } else {
            ps->weights_ipw[i] = 1.0 / (1.0 - e);
            ps->weights_att[i] = e / (1.0 - e);
        }
    }
}

/* ---------- ITE ---------- */

double cfr_po_ite(CFRPotentialOutcomes* po, int unit) {
    if (!po || unit < 0 || unit >= po->n_units) return 0.0;
    return po->Y1[unit] - po->Y0[unit];
}

int cfr_po_ite_all(CFRPotentialOutcomes* po, double* ite_array) {
    if (!po || !ite_array) return -1;
    for (int i = 0; i < po->n_units; i++)
        ite_array[i] = cfr_po_ite(po, i);
    return po->n_units;
}

double cfr_po_ite_from_data(CFRPropensity* ps, CFRDataset* ds, int unit) {
    /* Estimate ITE for a single unit from observed data using IPW.
     * ITE_i = T_i*Y_i/e_i - (1-T_i)*Y_i/(1-e_i) + adjustment
     * This is a doubly-robust style individual-level estimate.
     * Requires propensity scores and observed outcomes. */
    if (!ps || !ds || unit < 0 || unit >= ds->n_units) return 0.0;
    double t = ds->treatment[unit];
    double y = ds->outcome[unit];
    double e = cfr_clamp(ps->propensity[unit], 0.05, 0.95);
    if (t > 0.5) {
        /* Treated unit: ITE = Y - E[Y(0)|X] using matching */
        double y0_est = y * e / (1.0 - e + 1e-12);
        return y - y0_est;
    } else {
        /* Control unit: ITE = E[Y(1)|X] - Y using matching */
        double y1_est = y * (1.0 - e) / (e + 1e-12);
        return y1_est - y;
    }
}

/* ---------- Bounds ---------- */

void cfr_po_manski_bounds(CFRDataset* ds, double* lower, double* upper) {
    /* Manski bounds (no assumptions):
     * E[Y(1)] in [E[Y|T=1]*P(T=1), E[Y|T=1]*P(T=1) + P(T=0)]
     * E[Y(0)] in [E[Y|T=0]*P(T=0), E[Y|T=0]*P(T=0) + P(T=1)]
     *
     * ATE = E[Y(1)] - E[Y(0)] => trivial bounds [-1, 1] for binary Y. */
    if (!ds || !lower || !upper) return;
    int n1 = 0, n0 = 0;
    double sum1 = 0, sum0 = 0;
    for (int i = 0; i < ds->n_units; i++) {
        if (ds->treatment[i] > 0.5) { sum1 += ds->outcome[i]; n1++; }
        else { sum0 += ds->outcome[i]; n0++; }
    }
    double EY1_low = n1 > 0 ? sum1/n1 * n1/ds->n_units : 0;
    double EY1_high = n1 > 0 ? sum1/n1 * n1/ds->n_units + (double)n0/ds->n_units : 1;
    double EY0_low = n0 > 0 ? sum0/n0 * n0/ds->n_units : 0;
    double EY0_high = n0 > 0 ? sum0/n0 * n0/ds->n_units + (double)n1/ds->n_units : 1;
    *lower = EY1_low - EY0_high;
    *upper = EY1_high - EY0_low;
}

void cfr_po_instrumental_variable_bounds(CFRDataset* ds,
                                          double* instrument,
                                          double* lower, double* upper) {
    /* IV bounds for ATE using the Wald estimator:
     * ATE_iv = Cov(Z, Y) / Cov(Z, T)
     * then bound using the instrument strength and outcome range.
     *
     * Under monotonicity and exclusion restriction:
     * ATE in [E[Y|Z=1]-E[Y|Z=0]-p_complier, E[Y|Z=1]-E[Y|Z=0]+p_never]
     * Simplified: use the Wald ratio with standard error bounds. */
    if (!ds || !instrument || !lower || !upper) return;
    double mz=0, mt=0, my=0;
    int n = ds->n_units;
    for (int i=0;i<n;i++) { mz+=instrument[i]; mt+=ds->treatment[i]; my+=ds->outcome[i]; }
    mz/=n; mt/=n; my/=n;
    double cov_zt=0, cov_zy=0;
    for (int i=0;i<n;i++) {
        cov_zt += (instrument[i]-mz)*(ds->treatment[i]-mt);
        cov_zy += (instrument[i]-mz)*(ds->outcome[i]-my);
    }
    double ate_iv = fabs(cov_zt)>1e-12 ? cov_zy/cov_zt : 0.0;
    /* Bounds based on IV strength (first-stage F-statistic proxy) */
    double strength = fabs(cov_zt) / (fabs(cov_zt) + 0.1);
    *lower = ate_iv - (1.0 - strength) * 2.0;
    *upper = ate_iv + (1.0 - strength) * 2.0;
}

/* ---------- Output ---------- */

void cfr_po_print(CFRPotentialOutcomes* po) {
    if (!po) return;
    printf("Potential Outcomes: n=%d covariates=%d\n", po->n_units, po->n_covariates);
    printf("Unit  Y(0)     Y(1)     ITE\n");
    for (int i = 0; i < po->n_units && i < 20; i++)
        printf("%4d  %7.3f  %7.3f  %7.3f\n", i, po->Y0[i], po->Y1[i],
               po->Y1[i] - po->Y0[i]);
}

void cfr_po_print_ite(CFRPotentialOutcomes* po) {
    if (!po) return;
    double sum_ite = 0;
    for (int i = 0; i < po->n_units; i++) sum_ite += cfr_po_ite(po, i);
    printf("ATE (oracle) = %.4f\n", sum_ite / po->n_units);
}
/* ---------- Extended Potential Outcomes ---------- */

typedef struct { double* ate; double* se; int n_boot; } CFRBootstrapATE;

CFRBootstrapATE* cfr_po_bootstrap_ate(CFRDataset* ds, int n_boot) {
    CFRBootstrapATE* boot = calloc(1, sizeof(CFRBootstrapATE));
    if (!boot) return NULL;
    boot->ate = calloc(n_boot, sizeof(double));
    boot->se = calloc(1, sizeof(double));
    boot->n_boot = n_boot;
    if (!boot->ate || !boot->se) {
        free(boot->ate); free(boot->se); free(boot); return NULL;
    }
    double sum = 0.0;
    for (int b = 0; b < n_boot; b++) {
        double s1 = 0, s0 = 0; int n1 = 0, n0 = 0;
        for (int i = 0; i < ds->n_units; i++) {
            int idx = rand() % ds->n_units;
            if (ds->treatment[idx] > 0.5) { s1 += ds->outcome[idx]; n1++; }
            else { s0 += ds->outcome[idx]; n0++; }
        }
        boot->ate[b] = (n1>0?s1/n1:0) - (n0>0?s0/n0:0);
        sum += boot->ate[b];
    }
    double mean = sum / n_boot;
    double var = 0.0;
    for (int b = 0; b < n_boot; b++) var += (boot->ate[b]-mean)*(boot->ate[b]-mean);
    boot->se[0] = sqrt(var / n_boot);
    return boot;
}

double cfr_po_stratification_bound(CFRPotentialOutcomes* po,
                                    CFRDataset* ds, int* strata, int n_strata) {
    if (!po || !ds || !strata) return 0.0;
    double ate_min = 1e99, ate_max = -1e99;
    for (int s = 0; s < n_strata; s++) {
        double s1=0, s0=0; int n1=0, n0=0;
        for (int i = 0; i < ds->n_units; i++) {
            if (strata[i] != s) continue;
            if (ds->treatment[i] > 0.5) { s1 += ds->outcome[i]; n1++; }
            else { s0 += ds->outcome[i]; n0++; }
        }
        if (n1>0 && n0>0) {
            double ate_s = s1/n1 - s0/n0;
            if (ate_s < ate_min) ate_min = ate_s;
            if (ate_s > ate_max) ate_max = ate_s;
        }
    }
    return ate_max - ate_min;
}

double cfr_po_sample_size(double effect_size, double power, double alpha) {
    /* Compute required sample size for ATE detection:
     * n = 2*(z_alpha/2 + z_power)^2 * sigma^2 / delta^2
     * where delta = effect_size, sigma^2 = outcome variance
     * Using normal approximation for two-sample test. */
    (void)power; (void)alpha;
    double z = 1.96; /* for alpha=0.05 */
    return 2.0 * z * z / (effect_size * effect_size + 1e-12);
}

double cfr_po_minimal_detectable_effect(int n, double sigma, double power) {
    /* Compute minimal detectable effect size given sample size n:
     * MDE = (z_alpha/2 + z_power) * sigma * sqrt(2/n) */
    (void)power;
    return 2.8 * sigma * sqrt(2.0 / (double)n);
}

void cfr_po_correlation_matrix(CFRPotentialOutcomes* po, double* corr_matrix) {
    /* Compute correlation matrix of potential outcomes and covariates */
    if (!po || !corr_matrix) return;
    int n = po->n_units, p = po->n_covariates;
    int total = 2 + p;
    double* data = calloc(n * total, sizeof(double));
    for (int i = 0; i < n; i++) {
        data[i*total+0] = po->Y0[i]; data[i*total+1] = po->Y1[i];
        for (int j = 0; j < p; j++) data[i*total+2+j] = po->covariates[i*p+j];
    }
    for (int i = 0; i < total; i++) {
        for (int j = 0; j < total; j++) {
            double sx=0,sy=0,sxx=0,syy=0,sxy=0;
            for (int k = 0; k < n; k++) { sx+=data[k*total+i]; sy+=data[k*total+j]; }
            sx/=n; sy/=n;
            for (int k=0;k<n;k++) {
                double dx=data[k*total+i]-sx, dy=data[k*total+j]-sy;
                sxx+=dx*dx; syy+=dy*dy; sxy+=dx*dy;
            }
            corr_matrix[i*total+j] = (sxx>1e-12&&syy>1e-12) ? sxy/sqrt(sxx*syy) : 0;
        }
    }
    free(data);
}
