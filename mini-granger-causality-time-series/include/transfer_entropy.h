#ifndef TRANSFER_ENTROPY_H
#define TRANSFER_ENTROPY_H
#include "ts_core.h"

typedef struct {
    double te_xy, te_yx, hx, hy, mi_xy, p_value;
    int n_bins, n_samples;
    bool is_significant;
} TransferEntropy;

TransferEntropy* te_create(void);
void te_free(TransferEntropy* te);
int te_compute(TransferEntropy* te, const double* x, const double* y,
    int length, int n_bins, int k_history);
int te_significance(TransferEntropy* te, const double* x, const double* y,
    int length, int n_bins, int k_history, int n_surrogates, double alpha);
double entropy_histogram(const double* data, int n, int n_bins);
double mutual_information_histogram(const double* x, const double* y, int n, int n_bins);
void te_print(const TransferEntropy* te);
#endif
