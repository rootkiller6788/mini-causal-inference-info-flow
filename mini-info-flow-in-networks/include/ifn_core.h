#ifndef IFN_CORE_H
#define IFN_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * ifn_core.h — Information Flow in Networks Core
 *
 * Transfer entropy, Granger causality, mutual information on
 * networks — quantifying directional information flow.
 *
 * References:
 *   Schreiber (2000) "Measuring Information Transfer"
 *   Granger (1969) "Investigating Causal Relations"
 *   Runge et al. (2019) "Causal Network Reconstruction"
 * ============================================================== */

typedef enum { IFN_DISCRETE=0, IFN_CONTINUOUS=1, IFN_MIXED=2 } IFN_DataType;

/* ---- Time series ---- */
typedef struct {
    double* data;
    int length;
    int dim;
    double mean;
    double variance;
    double* histogram;
    int n_bins;
} IFN_TimeSeries;

/* ---- Transfer entropy result ---- */
typedef struct {
    double te_value;
    double te_surrogate_mean;
    double te_surrogate_std;
    double p_value;
    double effective_transfer;
    int source_id;
    int target_id;
    int k_history;
    int l_history;
} IFN_TransferResult;

/* ---- Granger causality result ---- */
typedef struct {
    double F_statistic;
    double p_value;
    double rss_restricted;
    double rss_unrestricted;
    int df1; int df2;
    bool is_causal;
    double* coefficients;
    int lag_order;
} IFN_GrangerResult;

/* ---- Causal graph (DAG) ---- */
typedef struct {
    int n_nodes;
    int* adjacency;      /* n×n matrix */
    int* moral_graph;    /* undirected moral graph */
    int* edges;          /* edge list: edges[2*e] */
    int n_edges;
    double* edge_weights; /* information flow weights */
} IFN_CausalGraph;

/* ---- Mutual information on network ---- */
typedef struct {
    double mi_value;
    double normalized_mi;  /* MI / min(H(X), H(Y)) */
    double* conditional_mi; /* MI conditioned on each other node */
    int n_nodes;
} IFN_MIResult;

/* ---- Network information dynamics ---- */
typedef struct {
    double* node_entropy;
    double* node_inflow;
    double* node_outflow;
    double* node_centrality;
    int n_nodes;
    int n_steps;
    double total_flow;
} IFN_NetworkDynamics;

/* Time series API */
IFN_TimeSeries* ifn_ts_create(int length, int dim);
void ifn_ts_free(IFN_TimeSeries* ts);
void ifn_ts_normalize(IFN_TimeSeries* ts);
double ifn_ts_mean(const IFN_TimeSeries* ts, int dim_idx);
double ifn_ts_variance(const IFN_TimeSeries* ts, int dim_idx);
void ifn_ts_histogram(IFN_TimeSeries* ts, int dim_idx, int n_bins);
double ifn_ts_entropy(IFN_TimeSeries* ts, int dim_idx);
void ifn_ts_delay_embed(const IFN_TimeSeries* ts, int dim_idx, int tau, int m, double** embedded);

/* Probability estimation from time series */
double ifn_estimate_pdf_bins(const double* x, int n, int n_bins, double* bins, double* pdf);
double ifn_estimate_joint_pdf(const double* x, const double* y, int n, int nx, int ny, double* joint);
double ifn_estimate_conditional_pdf(const double* x, const double* y, const double* z, int n, int nx, int ny, int nz, double* cond);

/* Internal MI computation helpers (used across modules) */
double ifn_mutual_information_binned(const double* x, const double* y, int n, int nx, int ny);
double ifn_conditional_entropy_binned(const double* x, const double* y, int n, int nx, int ny);

#endif
