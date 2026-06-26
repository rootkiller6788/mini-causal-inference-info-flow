#ifndef TE_KSG_H
#define TE_KSG_H
#include "te_core.h"

/* ================================================================
 * KSG (Kraskov-Stögbauer-Grassberger, 2004) k-Nearest-Neighbor
 * Continuous estimator for mutual information and transfer entropy.
 * Avoids binning artifacts via Kozachenko-Leonenko entropy estimation.
 * ================================================================ */

/* ---- Core KSG Estimators ---- */
TEResult te_ksg_compute(const TETimeSeries *x, const TETimeSeries *y, int k_embed, int l_embed, int k_nn);
double   te_ksg_estimate(const double *x, const double *y, int n, int k, int l, int knn);
double   te_ksg_entropy_knn(const double *data, int n, int dim, int knn);
double   te_ksg_mutual_info_knn(const double *x, const double *y, int n, int knn);
int      te_ksg_range_search(const double *points, int n, int dim, int idx, int k, int *indices, double *distances);

/* ---- Adaptive KSG Variants ---- */
double te_ksg_adaptive_k(const double *x, const double *y, int n);
double te_ksg_auto_knn(const double *x, const double *y, int n, int k_min, int k_max);
double te_ksg_entropy_loov(const double *data, int n, int dim, int k_min, int k_max);
double te_ksg_mi_manhattan(const double *x, const double *y, int n, int knn);

/* ---- Conditional and Transfer Entropy via CMI ---- */
double te_ksg_conditional_mi(const double *x, const double *y, const double *z, int n, int knn);
double te_ksg_te_via_cmi(const double *x, const double *y, int n, int k, int knn);

#endif
