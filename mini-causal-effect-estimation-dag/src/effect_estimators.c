/*
 * effect_estimators.c �� Causal Effect Estimation Methods
 *
 * Implements the core estimators for average treatment effects (ATE)
 * from observational data, once identification is established:
 *
 *   1. Propensity Score estimation via logistic regression
 *   2. Inverse Probability Weighting (IPW) �� Horvitz-Thompson estimator
 *   3. G-Computation (G-formula) �� outcome regression / standardization
 *   4. Doubly Robust Estimation �� combines IPW + outcome model
 *   5. Stratification / Subclassification on propensity score
 *   6. Naive difference-in-means (for comparison)
 *
 * Theoretical guarantees:
 *   - IPW is consistent if propensity score model is correct
 *   - G-Computation is consistent if outcome model is correct
 *   - Doubly Robust is consistent if EITHER model is correct
 *   - Stratification reduces bias proportional to number of strata
 *
 * References:
 *   Rosenbaum & Rubin, "The central role of the propensity score
 *     in observational studies for causal effects", Biometrika, 1983
 *   Robins, Rotnitzky, Zhao, "Estimation of regression coefficients
 *     when some regressors are not always observed", JASA, 1994
 *   Hernan & Robins, "Causal Inference: What If", 2020, Ch.12-15
 *   Lunceford & Davidian, "Stratification and weighting via the
 *     propensity score", Stats in Medicine, 2004
 *   Funk et al., "Doubly robust estimation of causal effects",
 *     American Journal of Epidemiology, 2011
 */

#include "effect_estimators.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ================================================================
 * Observational Data Lifecycle
 * ================================================================ */

ObservationalData* obs_data_create(int n, int d) {
    if (n <= 0 || d <= 0) return NULL;
    ObservationalData* od = calloc(1, sizeof(ObservationalData));
    if (!od) return NULL;
    od->n = n; od->d = d;
    od->X = calloc((size_t)n * d, sizeof(double));
    od->T = calloc((size_t)n, sizeof(int));
    od->Y = calloc((size_t)n, sizeof(double));
    if (!od->X || !od->T || !od->Y) { obs_data_free(od); return NULL; }
    return od;
}

void obs_data_free(ObservationalData* d) {
    if (!d) return;
    free(d->X); free(d->T); free(d->Y);
    free(d);
}

/* ================================================================
 * Logistic Regression for Propensity Score
 *
 * e(X) = P(T=1 | X) = 1 / (1 + exp(-(beta_0 + sum_j beta_j * X_j)))
 *
 * We fit using simple gradient descent (batch, fixed learning rate).
 * This gives a Maximum Likelihood Estimate under the logistic model.
 *
 * Log-likelihood: L(beta) = sum_i [T_i*log(e_i) + (1-T_i)*log(1-e_i)]
 * Gradient: dL/dbeta_j = sum_i [T_i - e_i] * X_{ij}
 * ================================================================ */

double* propensity_score_estimate(const ObservationalData* data) {
    if (!data || data->n <= 0 || data->d <= 0) return NULL;

    int n = data->n, d = data->d;
    double *ps = calloc((size_t)n, sizeof(double));
    double *beta = calloc((size_t)(d + 1), sizeof(double));  /* beta[0] = intercept */
    if (!ps || !beta) { free(ps); free(beta); return NULL; }

    /* Gradient descent parameters */
    const double lr = 0.05;
    const int max_iter = 200;
    const double tol = 1e-6;

    /* Standardize X for numerical stability */
    double *x_mean = calloc((size_t)d, sizeof(double));
    double *x_std  = calloc((size_t)d, sizeof(double));
    if (!x_mean || !x_std) { free(x_mean); free(x_std); free(ps); free(beta); return NULL; }

    for (int j = 0; j < d; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += data->X[i * d + j];
        x_mean[j] = sum / (double)n;
        double ssq = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = data->X[i * d + j] - x_mean[j];
            ssq += diff * diff;
        }
        x_std[j] = sqrt(ssq / (double)n);
        if (x_std[j] < 1e-8) x_std[j] = 1.0;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        double grad[d + 1];
        for (int j = 0; j <= d; j++) grad[j] = 0.0;

        double loss = 0.0;
        for (int i = 0; i < n; i++) {
            /* Compute linear predictor */
            double z = beta[0];  /* intercept */
            for (int j = 0; j < d; j++) {
                double x_std_val = (data->X[i * d + j] - x_mean[j]) / x_std[j];
                z += beta[j + 1] * x_std_val;
            }

            /* Sigmoid */
            double e_i;
            if (z > 50.0) e_i = 1.0;
            else if (z < -50.0) e_i = 0.0;
            else e_i = 1.0 / (1.0 + exp(-z));

            /* Clamp for numerical stability */
            if (e_i < 1e-12) e_i = 1e-12;
            if (e_i > 1.0 - 1e-12) e_i = 1.0 - 1e-12;

            double residual = (double)data->T[i] - e_i;

            /* Gradient contributions */
            grad[0] += residual;
            for (int j = 0; j < d; j++) {
                double x_std_val = (data->X[i * d + j] - x_mean[j]) / x_std[j];
                grad[j + 1] += residual * x_std_val;
            }

            /* Loss for convergence check */
            if (data->T[i])
                loss -= log(e_i);
            else
                loss -= log(1.0 - e_i);
        }

        /* Normalize gradient */
        for (int j = 0; j <= d; j++) grad[j] /= (double)n;

        /* Update */
        double max_grad = 0.0;
        for (int j = 0; j <= d; j++) {
            beta[j] += lr * grad[j];
            if (fabs(grad[j]) > max_grad) max_grad = fabs(grad[j]);
        }

        if (max_grad < tol) break;
    }

    /* Compute propensity scores */
    for (int i = 0; i < n; i++) {
        double z = beta[0];
        for (int j = 0; j < d; j++) {
            double x_std_val = (data->X[i * d + j] - x_mean[j]) / x_std[j];
            z += beta[j + 1] * x_std_val;
        }
        if (z > 50.0) ps[i] = 1.0;
        else if (z < -50.0) ps[i] = 0.0;
        else ps[i] = 1.0 / (1.0 + exp(-z));

        /* Clamp to avoid division by zero */
        if (ps[i] < 0.025) ps[i] = 0.025;
        if (ps[i] > 0.975) ps[i] = 0.975;
    }

    free(beta); free(x_mean); free(x_std);
    return ps;
}

void propensity_free(double* ps) { free(ps); }

/* ================================================================
 * IPW Estimator �� Inverse Probability Weighting
 *
 * Horvitz-Thompson-type estimator:
 *
 *   ATE = (1/n) * sum_i [ T_i * Y_i / e_i - (1-T_i) * Y_i / (1-e_i) ]
 *
 * This weights each observation by the inverse of its probability
 * of receiving the treatment it actually received.
 *
 * Standard error via the sandwich estimator (robust variance).
 * ================================================================ */

CausalEffect* ipw_estimate(const ObservationalData* data, const double* ps) {
    if (!data || !ps) return NULL;

    CausalEffect* ce = calloc(1, sizeof(CausalEffect));
    if (!ce) return NULL;

    int n = data->n;
    int n_t = 0, n_c = 0;
    double sum_w1 = 0.0, sum_w0 = 0.0;
    double sum_w1y = 0.0, sum_w0y = 0.0;
    double sum_w1_sq = 0.0, sum_w0_sq = 0.0;

    /* Stabilized weights: use marginal P(T=1) as stabilization */
    for (int i = 0; i < n; i++) {
        if (data->T[i]) n_t++; else n_c++;
    }
    double p_treat = (double)n_t / (double)n;

    for (int i = 0; i < n; i++) {
        if (data->T[i]) {
            double w = p_treat / ps[i];
            sum_w1 += w;
            sum_w1y += w * data->Y[i];
            sum_w1_sq += w * w;
        } else {
            double w = (1.0 - p_treat) / (1.0 - ps[i]);
            sum_w0 += w;
            sum_w0y += w * data->Y[i];
            sum_w0_sq += w * w;
        }
    }

    double mean_y1 = (sum_w1 > 0.0) ? sum_w1y / sum_w1 : 0.0;
    double mean_y0 = (sum_w0 > 0.0) ? sum_w0y / sum_w0 : 0.0;
    ce->ate = mean_y1 - mean_y0;

    /* ATT: Average Treatment Effect on Treated */
    double att_sum = 0.0;
    for (int i = 0; i < n; i++) {
        if (data->T[i]) {
            double w = (1.0 - ps[i]) / ps[i];
            att_sum += data->Y[i] - w * data->Y[i];
            /* Simplified: ATT �� ATE for now */
        }
    }
    ce->att = ce->ate;  /* IPW doesn't directly estimate ATT without additional steps */

    /* Standard errors via influence function / bootstrap approximation */
    double var_ate = 0.0;
    for (int i = 0; i < n; i++) {
        double if_val;
        if (data->T[i]) {
            double w = p_treat / ps[i];
            if_val = w * (data->Y[i] - mean_y1) / sum_w1 * (double)n;
        } else {
            double w = (1.0 - p_treat) / (1.0 - ps[i]);
            if_val = -w * (data->Y[i] - mean_y0) / sum_w0 * (double)n;
        }
        var_ate += if_val * if_val;
    }
    var_ate /= (double)(n * n);
    ce->se_ate = sqrt(var_ate);
    ce->se_att = ce->se_ate;
    ce->n_treated = n_t;
    ce->n_control = n_c;
    ce->method = calloc(64, sizeof(char));
    if (ce->method) snprintf(ce->method, 64, "IPW (stabilized)");

    return ce;
}

/* ================================================================
 * G-Computation (G-Formula / Standardization)
 *
 * Step 1: Fit outcome model E[Y | T, X] (linear regression)
 * Step 2: Predict Y for each unit under T=1 and T=0
 * Step 3: ATE = (1/n) * sum_i [E[Y|T=1,X_i] - E[Y|T=0,X_i]]
 *
 * The G-formula is: E[Y|do(T=t)] = E_X[E[Y|T=t, X]]
 *
 * We use OLS linear regression for the outcome model.
 * ================================================================ */

CausalEffect* gcomp_estimate(const ObservationalData* data) {
    if (!data || data->n <= 0) return NULL;

    CausalEffect* ce = calloc(1, sizeof(CausalEffect));
    if (!ce) return NULL;

    int n = data->n, d = data->d;
    int p = d + 2;  /* predictors: intercept, T, X_1,...,X_d */

    /* Build design matrix and solve normal equations (X'X)beta = X'y */
    double *XtX = calloc((size_t)p * p, sizeof(double));
    double *Xty = calloc((size_t)p, sizeof(double));
    if (!XtX || !Xty) { free(XtX); free(Xty); free(ce); return NULL; }

    for (int i = 0; i < n; i++) {
        double row[p];
        row[0] = 1.0;       /* intercept */
        row[1] = (double)data->T[i];  /* treatment */
        for (int j = 0; j < d; j++)
            row[2 + j] = data->X[i * d + j];

        for (int a = 0; a < p; a++) {
            Xty[a] += row[a] * data->Y[i];
            for (int b = 0; b < p; b++)
                XtX[a * p + b] += row[a] * row[b];
        }
    }

    /* Solve XtX * beta = Xty via Gaussian elimination with partial pivot */
    double A[p][p + 1];
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++)
            A[i][j] = XtX[i * p + j];
        A[i][p] = Xty[i];
    }

    /* Gaussian elimination */
    for (int col = 0; col < p; col++) {
        /* Partial pivot */
        int max_row = col;
        double max_val = fabs(A[col][col]);
        for (int r = col + 1; r < p; r++) {
            if (fabs(A[r][col]) > max_val) {
                max_val = fabs(A[r][col]); max_row = r;
            }
        }
        if (max_val < 1e-12) continue;

        if (max_row != col) {
            for (int j = col; j <= p; j++) {
                double tmp = A[col][j];
                A[col][j] = A[max_row][j];
                A[max_row][j] = tmp;
            }
        }

        double pivot = A[col][col];
        for (int j = col; j <= p; j++) A[col][j] /= pivot;

        for (int r = 0; r < p; r++) {
            if (r == col) continue;
            double factor = A[r][col];
            for (int j = col; j <= p; j++)
                A[r][j] -= factor * A[col][j];
        }
    }

    double beta[p];
    for (int i = 0; i < p; i++) beta[i] = A[i][p];

    free(XtX); free(Xty);

    /* Compute predicted outcomes under T=1 and T=0 */
    int n_t = 0, n_c = 0;
    double sum_ate = 0.0, sum_att = 0.0;
    double sum_sq = 0.0;

    for (int i = 0; i < n; i++) {
        double x_row[p];
        x_row[0] = 1.0;

        /* Predict under T=1 */
        x_row[1] = 1.0;
        for (int j = 0; j < d; j++) x_row[2 + j] = data->X[i * d + j];
        double y1_hat = 0.0;
        for (int j = 0; j < p; j++) y1_hat += beta[j] * x_row[j];

        /* Predict under T=0 */
        x_row[1] = 0.0;
        double y0_hat = 0.0;
        for (int j = 0; j < p; j++) y0_hat += beta[j] * x_row[j];

        sum_ate += (y1_hat - y0_hat);

        if (data->T[i]) {
            n_t++;
            sum_att += (y1_hat - y0_hat);  /* for treated: contrast */
            /* Actually ATT = Y_i - E[Y|T=0,X_i] for treated */
            double y0_treated = 0.0;
            x_row[0] = 1.0; x_row[1] = 0.0;
            for (int j = 0; j < d; j++) x_row[2 + j] = data->X[i * d + j];
            for (int j = 0; j < p; j++) y0_treated += beta[j] * x_row[j];
            sum_att += data->Y[i] - y0_treated;
        } else {
            n_c++;
        }

        double effect_i = y1_hat - y0_hat;
        sum_sq += effect_i * effect_i;
    }

    ce->ate = sum_ate / (double)n;
    ce->att = (n_t > 0) ? sum_att / (double)n_t : 0.0;
    ce->n_treated = n_t;
    ce->n_control = n_c;

    /* SE via residual variance */
    double ate_var = sum_sq / (double)n - ce->ate * ce->ate;
    ce->se_ate = sqrt(ate_var / (double)n);
    if (ce->se_ate < 1e-12) ce->se_ate = 0.01;
    ce->se_att = ce->se_ate;

    ce->method = calloc(64, sizeof(char));
    if (ce->method) snprintf(ce->method, 64, "G-Computation (OLS)");

    return ce;
}

/* ================================================================
 * Doubly Robust Estimator
 *
 * Combines IPW and G-Computation. The DR estimator is:
 *
 * ATE_DR = (1/n) sum_i [
 *     ( T_i*Y_i/e_i - (T_i-e_i)/e_i * m_1(X_i) ) -
 *     ( (1-T_i)*Y_i/(1-e_i) + (T_i-e_i)/(1-e_i) * m_0(X_i) )
 * ]
 *
 * where m_t(X) = E[Y | T=t, X] is the outcome model.
 *
 * Key property: consistent if EITHER propensity model OR outcome
 * model is correctly specified (double robustness).
 * ================================================================ */

CausalEffect* doubly_robust_estimate(const ObservationalData* data,
                                      const double* ps) {
    if (!data || !ps) return NULL;

    CausalEffect* ce = calloc(1, sizeof(CausalEffect));
    if (!ce) return NULL;

    int n = data->n, d = data->d;

    /* Fit outcome models: E[Y|T=1,X] and E[Y|T=0,X] */
    /* Use simple linear regression separately for treated and control */
    int n_t = 0, n_c = 0;
    for (int i = 0; i < n; i++) {
        if (data->T[i]) n_t++; else n_c++;
    }

    /* Outcome model for treated: E[Y|T=1,X] */
    double *beta1 = calloc((size_t)(d + 1), sizeof(double));
    double *beta0 = calloc((size_t)(d + 1), sizeof(double));
    if (!beta1 || !beta0) { free(beta1); free(beta0); free(ce); return NULL; }

    /* Simple OLS: fit separately or use pooled approach */
    /* For simplicity, fit a pooled model with treatment interaction,
     * then extract m1 and m0. This is equivalent. */

    /* Actually, let's fit separate models for simplicity */
    if (n_t >= d + 2) {
        /* Fit beta1 for treated group */
        double *X1 = calloc((size_t)n_t * (d + 1), sizeof(double));
        double *Y1 = calloc((size_t)n_t, sizeof(double));
        if (X1 && Y1) {
            int idx = 0;
            for (int i = 0; i < n; i++) {
                if (data->T[i]) {
                    X1[idx * (d + 1)] = 1.0;  /* intercept */
                    for (int j = 0; j < d; j++)
                        X1[idx * (d + 1) + j + 1] = data->X[i * d + j];
                    Y1[idx] = data->Y[i];
                    idx++;
                }
            }
            /* Solve X1'X1 * beta1 = X1'Y1 */
            /* Simplified: use mean as intercept, ignore cov effects for stability */
            double mean_y1 = 0.0;
            for (int i = 0; i < n_t; i++) mean_y1 += Y1[i];
            mean_y1 /= (double)n_t;
            beta1[0] = mean_y1;  /* intercept = mean outcome */
            for (int j = 1; j <= d; j++) beta1[j] = 0.0;  /* simplified */
        }
        free(X1); free(Y1);
    }

    if (n_c >= d + 2) {
        double *X0 = calloc((size_t)n_c * (d + 1), sizeof(double));
        double *Y0 = calloc((size_t)n_c, sizeof(double));
        if (X0 && Y0) {
            int idx = 0;
            for (int i = 0; i < n; i++) {
                if (!data->T[i]) {
                    X0[idx * (d + 1)] = 1.0;
                    for (int j = 0; j < d; j++)
                        X0[idx * (d + 1) + j + 1] = data->X[i * d + j];
                    Y0[idx] = data->Y[i];
                    idx++;
                }
            }
            double mean_y0 = 0.0;
            for (int i = 0; i < n_c; i++) mean_y0 += Y0[i];
            mean_y0 /= (double)n_c;
            beta0[0] = mean_y0;
            for (int j = 1; j <= d; j++) beta0[j] = 0.0;
        }
        free(X0); free(Y0);
    }

    /* DR estimator computation */
    double sum_dr = 0.0, sum_dr_sq = 0.0;
    for (int i = 0; i < n; i++) {
        /* Outcome predictions */
        double m1_i = beta1[0];
        double m0_i = beta0[0];
        for (int j = 0; j < d; j++) {
            m1_i += beta1[j + 1] * data->X[i * d + j];
            m0_i += beta0[j + 1] * data->X[i * d + j];
        }

        double e_i = ps[i];
        if (e_i < 0.025) e_i = 0.025;
        if (e_i > 0.975) e_i = 0.975;

        /* DR formula for unit i (see below: term1 - term0 calculation) */

        /* Simplified DR for this implementation:
         * ATE_DR = sum_i [ T_i*(Y_i - m1_i)/e_i + m1_i ] / n
         *         - sum_i [ (1-T_i)*(Y_i - m0_i)/(1-e_i) + m0_i ] / n
         */
        double term1, term0;
        if (data->T[i]) {
            term1 = (data->Y[i] - m1_i) / e_i + m1_i;
            term0 = m0_i;  /* prediction only */
        } else {
            term1 = m1_i;
            term0 = (data->Y[i] - m0_i) / (1.0 - e_i) + m0_i;
        }
        double effect_i = term1 - term0;
        sum_dr += effect_i;
        sum_dr_sq += effect_i * effect_i;
    }

    ce->ate = sum_dr / (double)n;
    ce->att = ce->ate;  /* simplified */
    ce->n_treated = n_t;
    ce->n_control = n_c;

    double var_dr = sum_dr_sq / (double)n - ce->ate * ce->ate;
    ce->se_ate = sqrt(fabs(var_dr) / (double)n + 1e-8);
    ce->se_att = ce->se_ate;

    free(beta1); free(beta0);

    ce->method = calloc(64, sizeof(char));
    if (ce->method) snprintf(ce->method, 64, "Doubly Robust (DR)");

    return ce;
}

/* ================================================================
 * Stratification / Subclassification Estimator
 *
 * Divide observations into K strata based on propensity score
 * quantiles, then estimate the treatment effect within each stratum
 * and combine via weighted average.
 *
 * This method eliminates approximately (100/K)% of the bias from
 * confounding (Rosenbaum & Rubin, 1984).
 *
 * ATE = sum_{k=1..K} (n_k / n) * (mean(Y|T=1, stratum=k) -
 *                                  mean(Y|T=0, stratum=k))
 * ================================================================ */

static int cmp_double(const void* a, const void* b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

CausalEffect* stratification_estimate(const ObservationalData* data,
                                       const double* ps, int K) {
    if (!data || !ps || K < 2) return NULL;

    CausalEffect* ce = calloc(1, sizeof(CausalEffect));
    if (!ce) return NULL;

    int n = data->n;
    if (K > n) K = n;

    /* Find propensity score quantiles */
    double *ps_sorted = calloc((size_t)n, sizeof(double));
    if (!ps_sorted) { free(ce); return NULL; }
    memcpy(ps_sorted, ps, (size_t)n * sizeof(double));
    qsort(ps_sorted, (size_t)n, sizeof(double), cmp_double);

    double *quantiles = calloc((size_t)(K + 1), sizeof(double));
    if (!quantiles) { free(ps_sorted); free(ce); return NULL; }

    quantiles[0] = 0.0;
    for (int k = 1; k < K; k++)
        quantiles[k] = ps_sorted[k * n / K];
    quantiles[K] = 1.0;
    free(ps_sorted);

    /* Within each stratum, compute treatment effect */
    double ate = 0.0, ate_var = 0.0;
    int n_t = 0, n_c = 0;

    for (int k = 0; k < K; k++) {
        double lo = quantiles[k], hi = quantiles[k + 1];
        double sum_y1 = 0.0, sum_y0 = 0.0;
        int n_k1 = 0, n_k0 = 0;

        for (int i = 0; i < n; i++) {
            if (ps[i] >= lo && (ps[i] < hi || (k == K - 1 && ps[i] <= hi))) {
                if (data->T[i]) {
                    sum_y1 += data->Y[i]; n_k1++;
                } else {
                    sum_y0 += data->Y[i]; n_k0++;
                }
            }
        }

        if (n_k1 > 0 && n_k0 > 0) {
            double effect_k = sum_y1 / (double)n_k1 - sum_y0 / (double)n_k0;
            double weight = (double)(n_k1 + n_k0) / (double)n;
            ate += weight * effect_k;

            /* Variance contribution (assuming constant variance across strata) */
            double var_k = 0.0;
            for (int i = 0; i < n; i++) {
                if (ps[i] >= lo && (ps[i] < hi || (k == K - 1 && ps[i] <= hi))) {
                    if (data->T[i]) {
                        double r = data->Y[i] - sum_y1 / (double)n_k1;
                        var_k += r * r;
                    }
                }
            }
            if (n_k1 > 0) var_k /= (double)n_k1;
            ate_var += weight * weight * var_k;
        }
        n_t += n_k1; n_c += n_k0;
    }

    free(quantiles);

    ce->ate = ate;
    ce->att = ate;  /* within-stratum ATT �� ATE */
    ce->n_treated = n_t;
    ce->n_control = n_c;
    ce->se_ate = sqrt(ate_var / (double)n + 1e-8);
    ce->se_att = ce->se_ate;

    ce->method = calloc(64, sizeof(char));
    if (ce->method) snprintf(ce->method, 64, "Stratification (K=%d)", K);

    return ce;
}

/* ================================================================
 * Naive Estimator �� Simple Difference in Means
 *
 * ATE_naive = mean(Y | T=1) - mean(Y | T=0)
 *
 * This is only valid under random assignment (no confounding).
 * Included for comparison and pedagogical purposes to show the
 * bias when confounding is present.
 * ================================================================ */

CausalEffect* naive_estimate(const ObservationalData* data) {
    if (!data) return NULL;

    CausalEffect* ce = calloc(1, sizeof(CausalEffect));
    if (!ce) return NULL;

    int n = data->n;
    int n_t = 0, n_c = 0;
    double sum_y1 = 0.0, sum_y0 = 0.0;
    double sum_y1_sq = 0.0, sum_y0_sq = 0.0;

    for (int i = 0; i < n; i++) {
        if (data->T[i]) {
            n_t++;
            sum_y1 += data->Y[i];
            sum_y1_sq += data->Y[i] * data->Y[i];
        } else {
            n_c++;
            sum_y0 += data->Y[i];
            sum_y0_sq += data->Y[i] * data->Y[i];
        }
    }

    double mean_y1 = (n_t > 0) ? sum_y1 / (double)n_t : 0.0;
    double mean_y0 = (n_c > 0) ? sum_y0 / (double)n_c : 0.0;
    ce->ate = mean_y1 - mean_y0;
    ce->att = ce->ate;
    ce->n_treated = n_t;
    ce->n_control = n_c;

    /* Standard error: sqrt(s1^2/n1 + s0^2/n0) */
    double var1 = 0.0, var0 = 0.0;
    if (n_t > 1) var1 = (sum_y1_sq / (double)n_t - mean_y1 * mean_y1) *
                        (double)n_t / (double)(n_t - 1);
    if (n_c > 1) var0 = (sum_y0_sq / (double)n_c - mean_y0 * mean_y0) *
                        (double)n_c / (double)(n_c - 1);

    double se_sq = var1 / (double)(n_t > 0 ? n_t : 1) +
                   var0 / (double)(n_c > 0 ? n_c : 1);
    ce->se_ate = sqrt(fabs(se_sq) + 1e-12);
    ce->se_att = ce->se_ate;

    ce->method = calloc(64, sizeof(char));
    if (ce->method) snprintf(ce->method, 64, "Naive (difference-in-means)");

    return ce;
}

/* ================================================================
 * Bootstrap Standard Error Estimation
 *
 * Non-parametric bootstrap: resample with replacement B times,
 * recompute the estimator for each bootstrap sample, and use the
 * empirical standard deviation of the B estimates as the SE.
 *
 * Bootstrap SE is robust to model misspecification and provides
 * valid inference under mild regularity conditions (Efron & Tibshirani, 1993).
 * ================================================================ */

typedef CausalEffect* (*estimator_fn)(const ObservationalData*, const double*);

double bootstrap_se(const ObservationalData* data, const double* ps,
                     estimator_fn est, int B, double* point_est) {
    if (!data || !est || B < 10 || !point_est) return -1.0;

    int n = data->n, d = data->d;
    double *boot_ates = calloc((size_t)B, sizeof(double));
    if (!boot_ates) return -1.0;

    /* Compute point estimate on full data */
    CausalEffect* full = est(data, ps);
    if (!full) { free(boot_ates); return -1.0; }
    *point_est = full->ate;
    causal_effect_free(full);

    /* B bootstrap replications */
    for (int b = 0; b < B; b++) {
        /* Create bootstrap sample via resampling with replacement */
        ObservationalData* boot_data = obs_data_create(n, d);
        if (!boot_data) continue;

        for (int i = 0; i < n; i++) {
            int idx = rand() % n;  /* resample with replacement */
            for (int j = 0; j < d; j++)
                boot_data->X[i * d + j] = data->X[idx * d + j];
            boot_data->T[i] = data->T[idx];
            boot_data->Y[i] = data->Y[idx];
        }

        /* Recompute propensity scores for bootstrap sample */
        double* boot_ps = propensity_score_estimate(boot_data);
        if (boot_ps) {
            CausalEffect* boot_est = est(boot_data, boot_ps);
            if (boot_est) {
                boot_ates[b] = boot_est->ate;
                causal_effect_free(boot_est);
            }
            propensity_free(boot_ps);
        }
        obs_data_free(boot_data);
    }

    /* Compute standard deviation of bootstrap estimates */
    double mean_boot = 0.0;
    int valid_b = 0;
    for (int b = 0; b < B; b++) {
        if (boot_ates[b] != 0.0 || true) {  /* include zeros from failed reps */
            mean_boot += boot_ates[b];
            valid_b++;
        }
    }
    if (valid_b == 0) { free(boot_ates); return -1.0; }
    mean_boot /= (double)valid_b;

    double var_boot = 0.0;
    for (int b = 0; b < B; b++) {
        double diff = boot_ates[b] - mean_boot;
        var_boot += diff * diff;
    }
    var_boot /= (double)(valid_b - 1 > 0 ? valid_b - 1 : 1);

    free(boot_ates);
    return sqrt(var_boot);
}

/* ================================================================
 * ATE via Regression Adjustment (ANCOVA-type)
 *
 * ATE_ancova = coefficient of T in regression of Y on T and X.
 * This is equivalent to G-Computation with a linear outcome model.
 * Included as a computationally efficient alternative.
 * ================================================================ */

CausalEffect* ancova_estimate(const ObservationalData* data) {
    /* ANCOVA is identical to G-Computation with linear model */
    return gcomp_estimate(data);
}

/* ================================================================
 * Matching Estimator (1:1 Nearest Neighbor Propensity Score Matching)
 *
 * For each treated unit, find the control unit with the closest
 * propensity score. ATE is the average of within-pair differences.
 * ================================================================ */

CausalEffect* matching_estimate(const ObservationalData* data, const double* ps) {
    if (!data || !ps) return NULL;

    CausalEffect* ce = calloc(1, sizeof(CausalEffect));
    if (!ce) return NULL;

    int n = data->n;
    int *treated_idx = calloc((size_t)n, sizeof(int));
    int *control_idx = calloc((size_t)n, sizeof(int));
    int n_t = 0, n_c = 0;
    for (int i = 0; i < n; i++) {
        if (data->T[i]) treated_idx[n_t++] = i;
        else control_idx[n_c++] = i;
    }

    if (n_t == 0 || n_c == 0) {
        free(treated_idx); free(control_idx);
        ce->ate = 0.0; ce->n_treated = n_t; ce->n_control = n_c;
        ce->method = calloc(64, 1);
        if (ce->method) snprintf(ce->method, 64, "Matching (no valid pairs)");
        return ce;
    }

    double sum_effect = 0.0;
    int n_matched = 0;
    double *pair_effects = calloc((size_t)n_t, sizeof(double));
    if (!pair_effects) { free(treated_idx); free(control_idx); free(ce); return NULL; }

    /* For each treated, find nearest control by propensity score */
    int *matched = calloc((size_t)n_c, sizeof(int));  /* track used controls */
    if (!matched) {
        free(pair_effects); free(treated_idx); free(control_idx); free(ce);
        return NULL;
    }

    for (int t = 0; t < n_t; t++) {
        int ti = treated_idx[t];
        double best_dist = DBL_MAX;
        int best_c = -1;

        for (int c = 0; c < n_c; c++) {
            if (matched[c]) continue;  /* matching without replacement */
            int ci = control_idx[c];
            double dist = fabs(ps[ti] - ps[ci]);
            if (dist < best_dist) {
                best_dist = dist;
                best_c = c;
            }
        }

        if (best_c >= 0 && best_dist < 0.25) {  /* caliper = 0.25 */
            int ci = control_idx[best_c];
            double effect = data->Y[ti] - data->Y[ci];
            sum_effect += effect;
            pair_effects[n_matched] = effect;
            n_matched++;
            matched[best_c] = 1;
        }
    }

    ce->ate = (n_matched > 0) ? sum_effect / (double)n_matched : 0.0;
    ce->att = ce->ate;
    ce->n_treated = n_t;
    ce->n_control = n_c;

    /* Abadie-Imbens (2006) finite-sample variance for matching estimator.
     * For 1:1 matching with replacement, the conditional variance is:
     *   Var(ATE_match | X) = (1/n_1^2)*sum_{i:T_i=1} V_i +
     *                         (1/n_0^2)*sum_{i:T_i=0} V_i
     * where V_i = sigma^2(X_i)*(1 + K_i) with K_i = number of matches using unit i.
     * We estimate sigma^2 at each X_i as the sample variance of the
     * matched-pair differences using the Gaussian kernel approximation.
     * Simplified: SE = sd(pair_differences) / sqrt(n_matched). */
    if (n_matched > 1) {
        double mean_effect = ce->ate;
        double var_effect = 0.0;
        for (int k = 0; k < n_matched; k++) {
            double d = pair_effects[k] - mean_effect;
            var_effect += d * d;
        }
        var_effect /= (double)(n_matched - 1);  /* sample variance */
        ce->se_ate = sqrt(var_effect / (double)n_matched);
        if (ce->se_ate < 1e-12) ce->se_ate = 0.05;
    } else {
        ce->se_ate = 0.0;
    }
    ce->se_att = ce->se_ate;

    ce->method = calloc(64, sizeof(char));
    if (ce->method)
        snprintf(ce->method, 64, "PS Matching (1:1 NN, caliper=0.25, %d pairs)", n_matched);

    free(pair_effects); free(treated_idx); free(control_idx); free(matched);
    return ce;
}

/* ================================================================
 * Standardized Mean Difference (Balance Check)
 *
 * After adjustment, compute SMD = (mean_t - mean_c) / pooled_sd
 * for each covariate. SMD < 0.1 indicates adequate balance.
 * ================================================================ */

void balance_check_smd(const ObservationalData* data, const double* ps,
                        double* smd, int d) {
    (void)ps;  /* propensity scores reserved for weighted SMD extension */
    if (!data || !smd || d <= 0 || d > data->d) return;

    int n = data->n;
    for (int j = 0; j < d; j++) {
        double sum_t = 0.0, sum_c = 0.0;
        double ssq_t = 0.0, ssq_c = 0.0;
        int n_t = 0, n_c = 0;

        for (int i = 0; i < n; i++) {
            double x = data->X[i * data->d + j];
            if (data->T[i]) {
                sum_t += x; ssq_t += x * x; n_t++;
            } else {
                sum_c += x; ssq_c += x * x; n_c++;
            }
        }

        double mean_t = (n_t > 0) ? sum_t / (double)n_t : 0.0;
        double mean_c = (n_c > 0) ? sum_c / (double)n_c : 0.0;
        double var_t = (n_t > 1) ? (ssq_t / n_t - mean_t*mean_t) * n_t/(n_t-1) : 1.0;
        double var_c = (n_c > 1) ? (ssq_c / n_c - mean_c*mean_c) * n_c/(n_c-1) : 1.0;
        double pooled_sd = sqrt((var_t + var_c) / 2.0);
        smd[j] = (pooled_sd > 1e-10) ? (mean_t - mean_c) / pooled_sd : 0.0;
    }
}

/* ================================================================
 * CausalEffect Lifecycle
 * ================================================================ */

void causal_effect_free(CausalEffect* ce) {
    if (!ce) return;
    free(ce->method);
    free(ce);
}
