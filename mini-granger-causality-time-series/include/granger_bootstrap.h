#ifndef GRANGER_BOOTSTRAP_H
#define GRANGER_BOOTSTRAP_H

/* Bootstrap inference for Granger causality.
 * Provides non-asymptotic p-values and confidence intervals
 * through parametric and non-parametric bootstrap resampling.
 *
 * Reference: Efron & Tibshirani (1993). An Introduction to the Bootstrap.
 *            Hacker & Hatemi-J (2006), Applied Economics, 38(13):1489-1500.
 */
#include "ts_core.h"

typedef struct {
    double observed_F;
    double bootstrap_p_value;
    double* surrogate_F;
    int n_bootstrap;
    double ci_lower;
    double ci_upper;
    double bias_correction;
} BootstrapGrangerResult;

BootstrapGrangerResult* boot_create(int n_bootstrap);
void boot_free(BootstrapGrangerResult* b);
int boot_parametric(BootstrapGrangerResult* b, const double* x, const double* y,
    int length, int max_lag, double alpha);
int boot_block(BootstrapGrangerResult* b, const double* x, const double* y,
    int length, int max_lag, double alpha, int block_len);
void boot_print(const BootstrapGrangerResult* b);

#endif
