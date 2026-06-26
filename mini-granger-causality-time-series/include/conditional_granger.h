#ifndef CONDITIONAL_GRANGER_H
#define CONDITIONAL_GRANGER_H
#include "ts_core.h"

typedef struct {
    double f_statistic, p_value, rss_without_x, rss_with_x;
    int lag_order, n_obs, n_conditioning;
    bool is_significant;
    double significance_level;
} ConditionalGranger;

ConditionalGranger* cg_create(void);
void cg_free(ConditionalGranger* cg);
int cg_test(ConditionalGranger* cg, const double* x, const double* y,
    const double* z, int length, int max_lag, double alpha);
int cg_test_multivariate(ConditionalGranger* cg, const double* x,
    const double* y, const TimeSeries* Z, int length, int max_lag, double alpha);
void cg_print(const ConditionalGranger* cg);
#endif
