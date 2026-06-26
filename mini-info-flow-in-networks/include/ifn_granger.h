#ifndef IFN_GRANGER_H
#define IFN_GRANGER_H
#include "ifn_core.h"

/* Granger causality: X Granger-causes Y if past values of X improve
 * prediction of Y beyond past values of Y alone.
 * Tested via F-test on VAR model residuals. */

IFN_GrangerResult ifn_granger_test(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, int max_lag);

IFN_GrangerResult ifn_granger_multivariate(const IFN_TimeSeries* multi_ts,
    int n_vars, int target_idx, int source_idx, int max_lag);

void ifn_granger_fit_var(const IFN_TimeSeries* y,
    const IFN_TimeSeries** predictors, int n_predictors, int lag,
    double* coefficients, double* residuals);

double ifn_granger_predict(const IFN_TimeSeries* y,
    const IFN_TimeSeries** predictors, int n_pred, int lag,
    const double* coeffs, int t_index);

void ifn_granger_aggregate_matrix(const IFN_TimeSeries* multi_ts,
    int n_vars, int max_lag, double* gc_matrix);

double ifn_granger_causality_index(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, int max_lag);

bool ifn_granger_is_significant(const IFN_GrangerResult* result, double alpha);

/* VAR OLS fitting: solve normal equations for VAR coefficients.
 * Ref: Hamilton (1994) Time Series Analysis, Chapter 11. */
void ifn_var_ols_fit(const double* y, int T, const double** pred,
    int n_pred, int lag, double* coeffs);

void ifn_granger_free(IFN_GrangerResult* gr);
#endif
