#ifndef GRANGER_TEST_H
#define GRANGER_TEST_H
#include "ts_core.h"

typedef struct {
    double f_statistic;
    double p_value;
    double rss_restricted;
    double rss_unrestricted;
    int lag_order;
    int n_obs;
    bool is_significant;
    double significance_level;
} GrangerTestResult;

GrangerTestResult* gt_create(void);
void gt_free(GrangerTestResult* gt);
int gt_test(GrangerTestResult* gt, const double* x, const double* y,
    int length, int max_lag, double alpha);
int gt_test_bidirectional(GrangerTestResult* gt_xy, GrangerTestResult* gt_yx,
    const double* x, const double* y, int length, int max_lag, double alpha);
double gt_f_distribution_pvalue(double f_stat, int df1, int df2);
int gt_optimal_lag(const double* x, const double* y, int length, int max_lag);
void gt_print(const GrangerTestResult* gt);
#endif
