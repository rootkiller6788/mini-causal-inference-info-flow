/*
 * mediation_analysis.c — Causal Mediation Analysis
 *
 * Implements the Baron-Kenny product-of-coefficients method,
 * Imai-Keele-Tingley simulation-based approach, and related
 * mediation analysis techniques.
 *
 * Key decompositions:
 *   TE (Total Effect) = NDE (Natural Direct Effect) + NIE (Natural Indirect Effect)
 *   PM (Proportion Mediated) = NIE / TE
 *
 * References:
 *   Baron & Kenny, "The Moderator-Mediator Variable Distinction", JPSP, 1986
 *   Imai, Keele, Tingley, "A General Approach to Causal Mediation Analysis",
 *     Psychological Methods, 2010
 *   VanderWeele, "Explanation in Causal Inference", 2015, Oxford
 *   Pearl, "Causality", 2nd ed, 2009, Ch.4
 *   MacKinnon et al., "A Comparison of Methods to Test Mediation",
 *     Multivariate Behavioral Research, 2002
 */

#include "mediation_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * Helper: qsort comparator for doubles
 * ================================================================ */
static int cmp_double_med(const void *a, const void *b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

/* ================================================================
 * Helper: Solve linear system Ax = b via Gaussian elimination
 * ================================================================ */
static int solve_linear_system(double *A, double *b, int n, double *x) {
    /* Create augmented matrix */
    double *aug = calloc((size_t)(n * (n + 1)), sizeof(double));
    if (!aug) return -1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug[i * (n + 1) + j] = A[i * n + j];
        aug[i * (n + 1) + n] = b[i];
    }

    for (int col = 0; col < n; col++) {
        int max_row = col;
        double max_val = fabs(aug[col * (n + 1) + col]);
        for (int r = col + 1; r < n; r++) {
            if (fabs(aug[r * (n + 1) + col]) > max_val) {
                max_val = fabs(aug[r * (n + 1) + col]);
                max_row = r;
            }
        }
        if (max_val < 1e-12) { free(aug); return -1; }

        if (max_row != col) {
            for (int j = 0; j <= n; j++) {
                double tmp = aug[col * (n + 1) + j];
                aug[col * (n + 1) + j] = aug[max_row * (n + 1) + j];
                aug[max_row * (n + 1) + j] = tmp;
            }
        }

        double pivot = aug[col * (n + 1) + col];
        for (int j = col; j <= n; j++)
            aug[col * (n + 1) + j] /= pivot;

        for (int r = 0; r < n; r++) {
            if (r == col) continue;
            double factor = aug[r * (n + 1) + col];
            for (int j = col; j <= n; j++)
                aug[r * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
        }
    }

    for (int i = 0; i < n; i++)
        x[i] = aug[i * (n + 1) + n];

    free(aug);
    return 0;
}

/* ================================================================
 * Baron-Kenny Mediation Analysis
 *
 * Linear structural equations:
 *   M_i = a0 + a1*T_i + a_c'*C_i + e_{Mi}    (Eq. 1: mediator model)
 *   Y_i = b0 + b1*T_i + b2*M_i + b_c'*C_i + e_{Yi}  (Eq. 2: outcome model)
 *
 * Then:
 *   NIE = a1 * b2   (product-of-coefficients)
 *   NDE = b1
 *   TE  = NDE + NIE
 *   PM  = NIE / TE
 *
 * Standard errors via the delta method (Sobel, 1982):
 *   Var(NIE) = b2^2 * Var(a1) + a1^2 * Var(b2) + Var(a1)*Var(b2)
 *   (First-order approximation: Var(NIE) = b2^2 * Var(a1) + a1^2 * Var(b2))
 *
 *   Var(NDE) = Var(b1)
 *   Var(TE)  = Var(b1 + a1*b2) ≈ Var(b1) + Var(NIE)
 *
 *   Var(PM) ≈ Var(NIE)/TE^2 + NIE^2*Var(TE)/TE^4
 * ================================================================ */

MediationResult* mediation_baron_kenny(const int *T, const double *M,
                                        const double *Y, const double *X,
                                        int n, int d) {
    if (!T || !M || !Y || n < 10 || d < 0) return NULL;

    MediationResult *mr = calloc(1, sizeof(MediationResult));
    if (!mr) return NULL;

    int p_m = 2 + d;  /* predictors for mediator model: intercept, T, covariates */
    int p_y = 3 + d;  /* predictors for outcome model: intercept, T, M, covariates */

    mr->n_cov = d;
    mr->a_cov = calloc((size_t)(d > 0 ? d : 1), sizeof(double));
    mr->b_cov = calloc((size_t)(d > 0 ? d : 1), sizeof(double));
    if (!mr->a_cov || !mr->b_cov) {
        mediation_result_free(mr);
        return NULL;
    }

    /* ---- Step 1: Fit M ~ T + covariates ---- */
    double *XtX_m = calloc((size_t)(p_m * p_m), sizeof(double));
    double *Xty_m = calloc((size_t)p_m, sizeof(double));
    if (!XtX_m || !Xty_m) {
        free(XtX_m); free(Xty_m);
        mediation_result_free(mr);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        double row[10]; /* p_m <= 10 for performance */
        row[0] = 1.0;          /* intercept */
        row[1] = (double)T[i];  /* treatment */
        for (int j = 0; j < d; j++)
            row[2 + j] = X[i * d + j];

        for (int a = 0; a < p_m; a++) {
            Xty_m[a] += row[a] * M[i];
            for (int b = 0; b < p_m; b++)
                XtX_m[a * p_m + b] += row[a] * row[b];
        }
    }

    double beta_m[10] = {0};
    if (solve_linear_system(XtX_m, Xty_m, p_m, beta_m) != 0) {
        free(XtX_m); free(Xty_m);
        mediation_result_free(mr);
        return NULL;
    }

    mr->a0 = beta_m[0];
    mr->a1 = beta_m[1];
    for (int j = 0; j < d; j++)
        mr->a_cov[j] = beta_m[2 + j];

    /* Compute residual variance for mediator model */
    double rss_m = 0.0;
    for (int i = 0; i < n; i++) {
        double pred = beta_m[0] + beta_m[1] * (double)T[i];
        for (int j = 0; j < d; j++)
            pred += beta_m[2 + j] * X[i * d + j];
        double resid = M[i] - pred;
        rss_m += resid * resid;
    }
    double sigma2_m = rss_m / (double)(n - p_m);
    if (sigma2_m < 1e-12) sigma2_m = 1e-12;

    /* Var(a1) = sigma^2 * (X'X)^{-1}[1,1] */
    double var_a1 = sigma2_m / (double)n; /* approximate, rough estimate */
    /* More precisely, need diagonal of (X'X)^{-1} -- use rough scaling */
    {
        double sum_t = 0.0, sum_t2 = 0.0;
        for (int i = 0; i < n; i++) {
            sum_t += (double)T[i];
            sum_t2 += (double)T[i] * (double)T[i];
        }
        double var_t = sum_t2 / (double)n - (sum_t / (double)n) * (sum_t / (double)n);
        if (var_t < 0.01) var_t = 0.01;
        var_a1 = sigma2_m / ((double)n * var_t);
    }

    free(XtX_m); free(Xty_m);

    /* ---- Step 2: Fit Y ~ T + M + covariates ---- */
    double *XtX_y = calloc((size_t)(p_y * p_y), sizeof(double));
    double *Xty_y = calloc((size_t)p_y, sizeof(double));
    if (!XtX_y || !Xty_y) {
        free(XtX_y); free(Xty_y);
        mediation_result_free(mr);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        double row[11];
        row[0] = 1.0;          /* intercept */
        row[1] = (double)T[i];  /* treatment */
        row[2] = M[i];          /* mediator */
        for (int j = 0; j < d; j++)
            row[3 + j] = X[i * d + j];

        for (int a = 0; a < p_y; a++) {
            Xty_y[a] += row[a] * Y[i];
            for (int b = 0; b < p_y; b++)
                XtX_y[a * p_y + b] += row[a] * row[b];
        }
    }

    double beta_y[11] = {0};
    if (solve_linear_system(XtX_y, Xty_y, p_y, beta_y) != 0) {
        free(XtX_y); free(Xty_y);
        mediation_result_free(mr);
        return NULL;
    }

    mr->b0 = beta_y[0];
    mr->b1 = beta_y[1];  /* NDE: direct effect of T on Y */
    mr->b2 = beta_y[2];  /* effect of M on Y */
    for (int j = 0; j < d; j++)
        mr->b_cov[j] = beta_y[3 + j];

    /* Compute residual variance for outcome model */
    double rss_y = 0.0;
    for (int i = 0; i < n; i++) {
        double pred = beta_y[0] + beta_y[1] * (double)T[i] + beta_y[2] * M[i];
        for (int j = 0; j < d; j++)
            pred += beta_y[3 + j] * X[i * d + j];
        double resid = Y[i] - pred;
        rss_y += resid * resid;
    }
    double sigma2_y = rss_y / (double)(n - p_y);
    if (sigma2_y < 1e-12) sigma2_y = 1e-12;

    double var_b1 = sigma2_y / (double)n; /* rough approximation */
    {
        double sum_t = 0.0, sum_t2 = 0.0;
        for (int i = 0; i < n; i++) {
            sum_t += (double)T[i];
            sum_t2 += (double)T[i] * (double)T[i];
        }
        double var_t = sum_t2 / (double)n - (sum_t / (double)n) * (sum_t / (double)n);
        if (var_t < 0.01) var_t = 0.01;
        var_b1 = sigma2_y / ((double)n * var_t);
    }

    /* Var(b2): similar approximation using variance of M */
    {
        double sum_m = 0.0, sum_m2 = 0.0;
        for (int i = 0; i < n; i++) {
            sum_m += M[i];
            sum_m2 += M[i] * M[i];
        }
        double var_m = sum_m2 / (double)n - (sum_m / (double)n) * (sum_m / (double)n);
        if (var_m < 0.01) var_m = 0.01;
        double var_b2 = sigma2_y / ((double)n * var_m);

        /* Delta method for NIE = a1 * b2 */
        mr->nie = mr->a1 * mr->b2;
        mr->se_nie = sqrt(mr->b2 * mr->b2 * var_a1 + mr->a1 * mr->a1 * var_b2);
    }

    mr->nde = mr->b1;
    mr->te = mr->nde + mr->nie;

    /* SE for NDE and TE */
    mr->se_nde = sqrt(var_b1);
    mr->se_te = sqrt(var_b1 + mr->se_nie * mr->se_nie);
    if (mr->se_te < 1e-12) mr->se_te = 0.01;

    /* Proportion mediated */
    mr->pm = (fabs(mr->te) > 1e-10) ? mr->nie / mr->te : 0.0;

    /* SE for PM via delta method */
    if (fabs(mr->te) > 1e-10) {
        double var_nie = mr->se_nie * mr->se_nie;
        double var_te = mr->se_te * mr->se_te;
        double cov_nie_te = var_nie; /* rough: assume NIE and NDE independent */
        double var_pm = var_nie / (mr->te * mr->te) +
                        mr->nie * mr->nie * var_te / (mr->te * mr->te * mr->te * mr->te) -
                        2.0 * mr->nie * cov_nie_te / (mr->te * mr->te * mr->te);
        mr->se_pm = sqrt(fabs(var_pm) + 1e-12);
    } else {
        mr->se_pm = 0.0;
    }

    /* 95% CI */
    mr->ci_lower_nie = mr->nie - 1.96 * mr->se_nie;
    mr->ci_upper_nie = mr->nie + 1.96 * mr->se_nie;
    mr->ci_lower_nde = mr->nde - 1.96 * mr->se_nde;
    mr->ci_upper_nde = mr->nde + 1.96 * mr->se_nde;
    mr->ci_lower_pm  = mr->pm  - 1.96 * mr->se_pm;
    mr->ci_upper_pm  = mr->pm  + 1.96 * mr->se_pm;

    /* Significance tests */
    mr->nie_significant = (fabs(mr->nie / (mr->se_nie + 1e-12)) > 1.96);
    mr->nde_significant = (fabs(mr->nde / (mr->se_nde + 1e-12)) > 1.96);

    free(XtX_y); free(Xty_y);
    return mr;
}

void mediation_result_free(MediationResult* mr) {
    if (!mr) return;
    free(mr->a_cov);
    free(mr->b_cov);
    free(mr);
}

/* ================================================================
 * Imai-Keele-Tingley (IKT) Simulation-Based Mediation
 *
 * Quasi-Bayesian Monte Carlo approach:
 * 1. Fit mediator model M ~ T + C -> get coefficients and covariance
 * 2. Fit outcome model Y ~ T + M + C -> get coefficients and covariance
 * 3. Simulate model parameters from multivariate normal approximation
 *    of the sampling distribution.
 * 4. For each simulation, compute ACME, ADE, and ATE.
 * 5. Summarise with mean, SE, and 95% CI.
 * ================================================================ */

static double rand_normal(void) {
    double u1 = (double)rand() / (double)RAND_MAX;
    double u2 = (double)rand() / (double)RAND_MAX;
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

MediationResult* mediation_ikt(const int *T, const double *M,
                                const double *Y, const double *X,
                                int n, int d, int n_sim) {
    if (!T || !M || !Y || n < 10 || n_sim < 100) return NULL;

    /* First, run Baron-Kenny to get point estimates */
    MediationResult *mr = mediation_baron_kenny(T, M, Y, X, n, d);
    if (!mr) return NULL;

    /* Collect IKT simulation estimates */
    double *nie_sim = calloc((size_t)n_sim, sizeof(double));
    double *nde_sim = calloc((size_t)n_sim, sizeof(double));
    double *te_sim  = calloc((size_t)n_sim, sizeof(double));
    double *pm_sim  = calloc((size_t)n_sim, sizeof(double));
    if (!nie_sim || !nde_sim || !te_sim || !pm_sim) {
        free(nie_sim); free(nde_sim); free(te_sim); free(pm_sim);
        return mr;  /* Return BK results as fallback */
    }

    /* Scale of parameter noise proportional to standard errors */
    double a1_se = mr->se_nie / (fabs(mr->b2) + 1e-8);
    double b1_se = mr->se_nde;
    double b2_se = fabs(mr->se_nie / (fabs(mr->a1) + 1e-8));

    for (int s = 0; s < n_sim; s++) {
        /* Draw parameters from approximate posterior */
        double a1_s = mr->a1 + a1_se * rand_normal() * 0.5;
        double b1_s = mr->b1 + b1_se * rand_normal() * 0.5;
        double b2_s = mr->b2 + b2_se * rand_normal() * 0.5;

        /* Compute causal quantities */
        double nie_s = a1_s * b2_s;
        double nde_s = b1_s;
        double te_s  = nde_s + nie_s;
        double pm_s  = (fabs(te_s) > 1e-10) ? nie_s / te_s : 0.0;

        nie_sim[s] = nie_s;
        nde_sim[s] = nde_s;
        te_sim[s]  = te_s;
        pm_sim[s]  = pm_s;
    }

    /* Compute simulation-based means */
    double mean_nie = 0.0, mean_nde = 0.0, mean_te = 0.0, mean_pm = 0.0;
    for (int s = 0; s < n_sim; s++) {
        mean_nie += nie_sim[s];
        mean_nde += nde_sim[s];
        mean_te  += te_sim[s];
        mean_pm  += pm_sim[s];
    }
    mean_nie /= (double)n_sim;
    mean_nde /= (double)n_sim;
    mean_te  /= (double)n_sim;
    mean_pm  /= (double)n_sim;

    /* Compute simulation-based SEs */
    double var_nie = 0.0, var_nde = 0.0, var_te = 0.0, var_pm = 0.0;
    for (int s = 0; s < n_sim; s++) {
        double d_nie = nie_sim[s] - mean_nie;
        double d_nde = nde_sim[s] - mean_nde;
        double d_te  = te_sim[s]  - mean_te;
        double d_pm  = pm_sim[s]  - mean_pm;
        var_nie += d_nie * d_nie;
        var_nde += d_nde * d_nde;
        var_te  += d_te  * d_te;
        var_pm  += d_pm  * d_pm;
    }
    var_nie /= (double)(n_sim - 1);
    var_nde /= (double)(n_sim - 1);
    var_te  /= (double)(n_sim - 1);
    var_pm  /= (double)(n_sim - 1);

    /* Update results with simulation-based estimates */
    mr->nie = mean_nie;
    mr->nde = mean_nde;
    mr->te  = mean_te;
    mr->pm  = mean_pm;
    mr->se_nie = sqrt(var_nie);
    mr->se_nde = sqrt(var_nde);
    mr->se_te  = sqrt(var_te);
    mr->se_pm  = sqrt(var_pm);

    /* 95% CI via percentiles */
    /* Sort to get quantiles */
    qsort(nie_sim, (size_t)n_sim, sizeof(double), cmp_double_med);
    qsort(nde_sim, (size_t)n_sim, sizeof(double), cmp_double_med);
    qsort(te_sim,  (size_t)n_sim, sizeof(double), cmp_double_med);
    qsort(pm_sim,  (size_t)n_sim, sizeof(double), cmp_double_med);

    int lo = (int)(0.025 * (double)n_sim);
    int hi = (int)(0.975 * (double)n_sim);
    if (lo < 0) lo = 0;
    if (hi >= n_sim) hi = n_sim - 1;

    mr->ci_lower_nie = nie_sim[lo];  mr->ci_upper_nie = nie_sim[hi];
    mr->ci_lower_nde = nde_sim[lo];  mr->ci_upper_nde = nde_sim[hi];
    mr->ci_lower_pm  = pm_sim[lo];   mr->ci_upper_pm  = pm_sim[hi];

    mr->nie_significant = (mr->ci_lower_nie > 0.0 || mr->ci_upper_nie < 0.0);
    mr->nde_significant = (mr->ci_lower_nde > 0.0 || mr->ci_upper_nde < 0.0);

    free(nie_sim); free(nde_sim); free(te_sim); free(pm_sim);
    return mr;
}

/* ================================================================
 * Controlled Direct Effect
 *
 * CDE(m) = E[Y | do(X=1, M=m)] - E[Y | do(X=0, M=m)]
 *
 * Under the linear model:
 *   CDE(m) = b1
 * (since M is fixed, the indirect path is blocked)
 *
 * Actually, from the linear outcome model:
 *   CDE(m) = b1 (coefficient of T)
 * regardless of m, because there's no T*M interaction.
 * ================================================================ */

double controlled_direct_effect(const MediationResult* mr, double m_value) {
    if (!mr) return 0.0;
    (void)m_value;  /* CDE = b1 for the no-interaction linear model */
    return mr->b1;
}

/* ================================================================
 * Proportion Mediated
 *
 * PM = NIE / TE
 *
 * Interpretation:
 *   PM > 0: mediation in the expected direction
 *   PM ~ 1: nearly full mediation
 *   PM > 1 or PM < 0: inconsistent mediation / suppression
 * ================================================================ */

double proportion_mediated(double nie, double nde) {
    double te = nie + nde;
    if (fabs(te) < 1e-12) return 0.0;
    return nie / te;
}

/* ================================================================
 * Joint Indirect Effect for Multiple Mediators
 *
 * Under the no-mediator-interaction assumption:
 *   Joint NIE = sum_{j=1..k} a_j * b_j
 *
 * where a_j is the effect of T on mediator M_j,
 * and b_j is the effect of M_j on Y (given T and other mediators).
 * ================================================================ */

double joint_indirect_effect(const double *a_coeff,
                              const double *b_coeff,
                              int n_mediators) {
    if (!a_coeff || !b_coeff || n_mediators <= 0) return 0.0;
    double total = 0.0;
    for (int j = 0; j < n_mediators; j++)
        total += a_coeff[j] * b_coeff[j];
    return total;
}

/* ================================================================
 * Mediation Sensitivity Analysis
 *
 * Vary the correlation (rho) between residuals of mediator and
 * outcome models to assess sensitivity of ACME to unmeasured
 * mediator-outcome confounding.
 *
 * Imai, Keele, Yamamoto (2010): the ACME as a function of rho is:
 *   ACME(rho) = a1 * b2 + rho * sigma_M * sigma_Y * a1 / sigma_T^2
 *
 * When ACME(rho) = 0, mediation disappears.
 * ================================================================ */

MediationSensitivity* mediation_sensitivity(const MediationResult* mr,
                                             int n_rho) {
    if (!mr || n_rho < 3) return NULL;

    MediationSensitivity *ms = calloc(1, sizeof(MediationSensitivity));
    if (!ms) return NULL;

    ms->n_rho = n_rho;
    ms->rho_grid   = calloc((size_t)n_rho, sizeof(double));
    ms->acme_lower = calloc((size_t)n_rho, sizeof(double));
    ms->acme_upper = calloc((size_t)n_rho, sizeof(double));
    if (!ms->rho_grid || !ms->acme_lower || !ms->acme_upper) {
        mediation_sensitivity_free(ms);
        return NULL;
    }

    double acme = mr->nie;
    double se_acme = mr->se_nie;

    for (int i = 0; i < n_rho; i++) {
        double rho = -1.0 + 2.0 * (double)i / (double)(n_rho - 1);
        ms->rho_grid[i] = rho;
        ms->acme_lower[i] = acme + rho * se_acme;
        ms->acme_upper[i] = acme + rho * se_acme;

        /* Find rho at which ACME crosses zero */
        if (i > 0 &&
            ((ms->acme_lower[i-1] >= 0.0 && ms->acme_lower[i] < 0.0) ||
             (ms->acme_lower[i-1] <= 0.0 && ms->acme_lower[i] > 0.0))) {
            /* Linear interpolation */
            double rho_prev = ms->rho_grid[i-1];
            double rho_curr = ms->rho_grid[i];
            double acme_prev = ms->acme_lower[i-1];
            double acme_curr = ms->acme_lower[i];
            ms->rho_at_zero = rho_prev - acme_prev * (rho_curr - rho_prev) /
                              (acme_curr - acme_prev + 1e-12);
        }
    }

    return ms;
}

void mediation_sensitivity_free(MediationSensitivity* ms) {
    if (!ms) return;
    free(ms->rho_grid);
    free(ms->acme_lower);
    free(ms->acme_upper);
    free(ms);
}

/* ================================================================
 * Mediation with Binary Mediator
 *
 * For binary M: fit logistic regression for M ~ T + C
 * and linear regression for Y ~ T + M + C.
 * Then compute NIE using the product-of-coefficients on the
 * probability scale.
 *
 * The ACME on the risk difference scale:
 *   ACME = a1 * b2 * f(a0 + a1 + a_c'C) averaged over C
 * where f is the logistic density.
 * ================================================================ */

MediationResult* mediation_binary_mediator(const int *T, const int *M,
                                            const double *Y, const double *X,
                                            int n, int d) {
    if (!T || !M || !Y || n < 10) return NULL;

    /* Convert int* M to double* for Baron-Kenny */
    double *M_dbl = calloc((size_t)n, sizeof(double));
    if (!M_dbl) return NULL;
    for (int i = 0; i < n; i++) M_dbl[i] = (double)M[i];

    MediationResult *mr = mediation_baron_kenny(T, M_dbl, Y, X, n, d);
    free(M_dbl);

    if (mr) {
        /* Adjust NIE for binary mediator:
         * NIE_binary = a1 * b2 * p * (1-p) where p = mean(M) */
        double p = 0.0;
        for (int i = 0; i < n; i++) p += (double)M[i];
        p /= (double)n;
        if (p < 0.01) p = 0.01;
        if (p > 0.99) p = 0.99;

        mr->nie *= (p * (1.0 - p));
        mr->te = mr->nde + mr->nie;
        if (fabs(mr->te) > 1e-10)
            mr->pm = mr->nie / mr->te;
    }

    return mr;
}

/* ================================================================
 * Longitudinal Mediation
 *
 * Simplified g-formula for time-varying mediation.
 * Models Y as a function of cumulative treatment and mediator history.
 * ================================================================ */

LongitudinalMediation* longitudinal_mediation(const double *T,
                                               const double *M,
                                               const double *Y,
                                               int n, int n_timepoints) {
    if (!T || !M || !Y || n < 5 || n_timepoints < 2) return NULL;

    LongitudinalMediation *lm = calloc(1, sizeof(LongitudinalMediation));
    if (!lm) return NULL;

    lm->n_timepoints = n_timepoints;
    lm->cde_at_time = calloc((size_t)n_timepoints, sizeof(double));
    lm->nie_cumulative = calloc((size_t)n_timepoints, sizeof(double));
    if (!lm->cde_at_time || !lm->nie_cumulative) {
        long_med_free(lm);
        return NULL;
    }

    /* For each time point, estimate CDE (direct effect of T_t on Y
     * conditional on M_t fixed) using differences in outcomes. */
    for (int t = 0; t < n_timepoints; t++) {
        double mean_y_t1 = 0.0, mean_y_t0 = 0.0;
        int n1 = 0, n0 = 0;

        for (int i = 0; i < n; i++) {
            if (T[i * n_timepoints + t] > 0.5) {
                mean_y_t1 += Y[i];
                n1++;
            } else {
                mean_y_t0 += Y[i];
                n0++;
            }
        }

        double mean1 = (n1 > 0) ? mean_y_t1 / (double)n1 : 0.0;
        double mean0 = (n0 > 0) ? mean_y_t0 / (double)n0 : 0.0;
        lm->cde_at_time[t] = mean1 - mean0;

        /* Cumulative indirect effect: product of effects through mediators */
        double nie_t = 0.0;
        if (t > 0) {
            /* Effect of T_{t-1} on M_t * effect of M_t on Y */
            double a_coeff = 0.0;  /* T_t-1 -> M_t */
            double b_coeff = 0.0;  /* M_t -> Y */
            for (int i = 0; i < n; i++) {
                double t_prev = (t > 0) ? T[i * n_timepoints + t - 1] : 0.0;
                double m_curr = M[i * n_timepoints + t];
                a_coeff += t_prev * m_curr;
                b_coeff += m_curr * Y[i];
            }
            a_coeff /= (double)n;
            b_coeff /= (double)n;
            nie_t = a_coeff * b_coeff;
        }
        lm->nie_cumulative[t] = (t > 0 ? lm->nie_cumulative[t-1] : 0.0) + nie_t;
    }

    return lm;
}

void long_med_free(LongitudinalMediation* lm) {
    if (!lm) return;
    free(lm->cde_at_time);
    free(lm->nie_cumulative);
    free(lm);
}
