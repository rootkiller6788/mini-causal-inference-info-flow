#ifndef IFN_MUTUAL_H
#define IFN_MUTUAL_H
#include "ifn_core.h"

/* Mutual information: I(X;Y) = H(X) + H(Y) - H(X,Y)
 * On networks, MI decomposes into direct and indirect flows. */

IFN_MIResult ifn_mutual_information(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, int n_bins);

double ifn_mi_continuous(const double* x, const double* y, int n, int n_bins);
double ifn_mi_knn(const double* x, const double* y, int n, int k);
double ifn_mi_kernel(const double* x, const double* y, int n, double bandwidth);

void ifn_mi_conditioned(const IFN_TimeSeries* multi_ts,
    int n_vars, int idx_x, int idx_y, int n_bins, double* cmi_out);

void ifn_mi_matrix(const IFN_TimeSeries* multi_ts, int n_vars,
    int n_bins, double* mi_matrix);

double ifn_mi_normalized(const IFN_TimeSeries* x, const IFN_TimeSeries* y,
    int n_bins);

double ifn_partial_mutual_information(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, const IFN_TimeSeries* z, int n_bins);

/* Interaction information: I(X;Y;Z) = I(X;Y) - I(X;Y|Z) */
double ifn_interaction_information(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, const IFN_TimeSeries* z, int n_bins);

/* Redundancy: R = I(X;Y) + I(X;Z) - I(X; Y,Z) */
double ifn_redundancy(const IFN_TimeSeries* x, const IFN_TimeSeries* y,
    const IFN_TimeSeries* z, int n_bins);

#endif
