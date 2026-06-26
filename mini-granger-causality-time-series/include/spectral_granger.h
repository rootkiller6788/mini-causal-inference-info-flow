#ifndef SPECTRAL_GRANGER_H
#define SPECTRAL_GRANGER_H
#include "ts_core.h"

typedef struct {
    double* frequencies;
    double* gc_x_to_y;
    double* gc_y_to_x;
    int n_freqs;
    double bandwidth;
} SpectralGranger;

SpectralGranger* sg_create(int n_freqs);
void sg_free(SpectralGranger* sg);
int sg_compute(SpectralGranger* sg, const double* x, const double* y,
    int length, int n_freqs, int max_lag);
int sg_geweke_decomposition(SpectralGranger* sg,
    const double* x, const double* y, int length, int n_freqs, int p);
void sg_print(const SpectralGranger* sg);
#endif
