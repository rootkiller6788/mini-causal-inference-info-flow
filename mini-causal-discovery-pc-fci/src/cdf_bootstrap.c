/* cdf_bootstrap.c - Bootstrap-Based Edge Stability Assessment
 *
 * Implements non-parametric bootstrap resampling for assessing
 * the stability of edges discovered by PC/FCI algorithms.
 *
 * Bootstrap procedure (Efron, 1979; Friedman et al., 1999):
 *   1. Draw B bootstrap samples by resampling N rows with replacement
 *   2. Run PC (or FCI) on each sample
 *   3. For each pair (i,j), record the edge type frequency
 *   4. Compute stability score = appearance_count / B
 *
 * Stability selection (Meinshausen & Buhlmann, 2010):
 *   P(edge selected) = sum I(edge in G_b) / B
 *   An edge is stable if its selection probability >= pi_thr.
 */

#include "cdf_bootstrap.h"
#include "cdf_citest.h"
#include "cdf_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* Simple LCG random number generator (portable, deterministic) */
typedef struct {
    unsigned long state;
} cdf_rng;

static cdf_rng* cdf_rng_create(int seed) {
    cdf_rng *r = (cdf_rng*)malloc(sizeof(cdf_rng));
    if (!r) return NULL;
    r->state = (unsigned long)(seed > 0 ? seed : 1);
    return r;
}

/* LCG: a=1664525, c=1013904223, m=2^32 (Numerical Recipes) */
static double cdf_rng_uniform(cdf_rng *r) {
    r->state = r->state * 1664525UL + 1013904223UL;
    return (double)(r->state & 0x7FFFFFFF) / 2147483648.0;
}

static int cdf_rng_int_range(cdf_rng *r, int lo, int hi) {
    return lo + (int)(cdf_rng_uniform(r) * (double)(hi - lo));
}

/* --------------- Bootstrap Resampling --------------- */

CdfDataset* cdf_bootstrap_resample(const CdfDataset *ds)
{
    if (!ds || ds->N < 1) return NULL;

    int N = ds->N, p = ds->p;
    size_t n_elem = (size_t)p * (size_t)N;
    double *new_data = (double*)malloc(n_elem * sizeof(double));
    if (!new_data) return NULL;

    cdf_rng *rng = cdf_rng_create((int)time(NULL) + N + p);
    if (!rng) { free(new_data); return NULL; }

    for (int i = 0; i < p; i++) {
        double *src_col = ds->data + (size_t)i * N;
        double *dst_col = new_data + (size_t)i * N;
        for (int j = 0; j < N; j++) {
            int idx = cdf_rng_int_range(rng, 0, N);
            dst_col[j] = src_col[idx];
        }
    }

    free(rng);
    CdfDataset *res = cdf_dataset_create(new_data, N, p);
    free(new_data);
    return res;
}

/* --------------- Bootstrap Configuration --------------- */

CdfBootstrapConfig cdf_bootstrap_config_default(void)
{
    CdfBootstrapConfig cfg;
    cfg.n_bootstrap = 100;
    cfg.stability_thresh = 0.7;
    cfg.random_seed = 42;
    cfg.use_pc = true;
    cfg.verbose = false;
    return cfg;
}

/* Edge orientation recording helper */
static void record_edge_orientation(CdfEdgeType type, int *u2v_dir,
                                     int *v2u_dir, int *undir,
                                     int *bidir, int *appear)
{
    (*appear)++;
    switch (type) {
    case CDF_EDGE_DIRECTED:   (*u2v_dir)++; break;
    case CDF_EDGE_UNDIRECTED: (*undir)++;   break;
    case CDF_EDGE_BIDIRECTED: (*bidir)++;   break;
    case CDF_EDGE_PARTIAL_I:  (*u2v_dir)++; break;
    case CDF_EDGE_PARTIAL_J:  (*v2u_dir)++; break;
    case CDF_EDGE_NONDIR:     (*undir)++;   break;
    default: break;
    }
}

/* --------------- Main Bootstrap Algorithm --------------- */

CdfBootstrapResult* cdf_bootstrap_edges(const CdfDataset *ds,
                                         const CdfBootstrapConfig *config)
{
    if (!ds || !config || config->n_bootstrap < 1) return NULL;

    CdfBootstrapConfig cfg = *config;
    int B = cfg.n_bootstrap, p = ds->p;
    int max_pairs = p * (p - 1) / 2;

    CdfBootstrapResult *res = (CdfBootstrapResult*)calloc(1, sizeof(CdfBootstrapResult));
    if (!res) return NULL;

    res->edges = (CdfBootstrapEdge*)calloc((size_t)max_pairs, sizeof(CdfBootstrapEdge));
    if (!res->edges) { free(res); return NULL; }

    int idx = 0;
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            res->edges[idx].u = i;
            res->edges[idx].v = j;
            idx++;
        }
    }
    res->n_edges = max_pairs;
    res->n_bootstrap = B;
    res->p = p;

    clock_t t_start = clock();

    for (int b = 0; b < B; b++) {
        CdfDataset *bs_ds = cdf_bootstrap_resample(ds);
        if (!bs_ds) continue;

        CdfGraph *bs_graph = NULL;
        if (cfg.use_pc) {
            CdfPCConfig pc_cfg = cdf_pc_config_default();
            pc_cfg.alpha = CDF_ALPHA_DEFAULT;
            pc_cfg.max_cond_size = -1;
            pc_cfg.verbose = false;
            CdfPCResult *pc_res = cdf_pc_run(bs_ds, &pc_cfg);
            if (pc_res) { bs_graph = pc_res->graph; free(pc_res); }
        } else {
            CdfFCIConfig fci_cfg = cdf_fci_config_default();
            fci_cfg.alpha = CDF_ALPHA_DEFAULT;
            fci_cfg.max_cond_size = -1;
            fci_cfg.verbose = false;
            CdfFCIResult *fci_res = cdf_fci_run(bs_ds, &fci_cfg);
            if (fci_res) { bs_graph = fci_res->graph; free(fci_res); }
        }

        if (bs_graph) {
            idx = 0;
            for (int i = 0; i < p; i++) {
                for (int j = i + 1; j < p; j++) {
                    CdfEdgeType e_ij = bs_graph->edges[i * p + j];
                    if (e_ij != CDF_EDGE_NONE) {
                        record_edge_orientation(e_ij,
                            &res->edges[idx].n_directed_uv,
                            &res->edges[idx].n_directed_vu,
                            &res->edges[idx].n_undirected,
                            &res->edges[idx].n_bidirected,
                            &res->edges[idx].n_appearances);
                    }
                    idx++;
                }
            }
            cdf_graph_free(bs_graph);
        }

        cdf_dataset_free(bs_ds);

        if (cfg.verbose && (b + 1) % 10 == 0)
            printf("  Bootstrap: %d/%d done\n", b + 1, B);
    }

    for (int k = 0; k < max_pairs; k++) {
        res->edges[k].freq_stable = (double)res->edges[k].n_appearances / (double)B;
        res->edges[k].is_stable = (res->edges[k].freq_stable >= cfg.stability_thresh);
    }

    res->elapsed_sec = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    return res;
}

void cdf_bootstrap_result_free(CdfBootstrapResult *res)
{
    if (!res) return;
    free(res->edges);
    free(res);
}

void cdf_bootstrap_print(const CdfBootstrapResult *res)
{
    if (!res) return;
    printf("=== Bootstrap Edge Stability (B=%d) ===\n", res->n_bootstrap);
    printf("  p=%d, edge pairs=%d\n", res->p, res->n_edges);
    int n_stable = 0;
    for (int k = 0; k < res->n_edges; k++)
        if (res->edges[k].is_stable) n_stable++;
    printf("  Stable edges: %d\n", n_stable);

    if (n_stable > 0) {
        printf("  Edge  Freq    Dir(u->v) Dir(v->u) Undir Bidir\n");
        for (int k = 0; k < res->n_edges; k++) {
            if (!res->edges[k].is_stable) continue;
            printf("  %d-%d  %.3f  %-9d %-9d %-5d %d\n",
                   res->edges[k].u, res->edges[k].v,
                   res->edges[k].freq_stable,
                   res->edges[k].n_directed_uv,
                   res->edges[k].n_directed_vu,
                   res->edges[k].n_undirected,
                   res->edges[k].n_bidirected);
        }
    }
    printf("  Time: %.3f sec\n", res->elapsed_sec);
}

/* --------------- FDR Threshold --------------- */

double cdf_bootstrap_fdr_threshold(const CdfBootstrapEdge *edges,
                                    int n_edges, double fdr_target)
{
    if (!edges || n_edges < 1) return 1.0;
    if (fdr_target <= 0.0) return 1.0;
    if (fdr_target >= 1.0) return 0.5;

    /* Stability selection bound:
     *   E[V] <= 1/(2*pi - 1) * q_pi^2 / p
     * Heuristic threshold:
     *   pi_thr = 0.5 + 0.5 * sqrt(1 - alpha * avg_q)
     */
    double total_freq_sum = 0.0;
    for (int k = 0; k < n_edges; k++)
        total_freq_sum += edges[k].freq_stable;

    double avg_q = total_freq_sum / (double)n_edges;
    double factor = 0.5 * (1.0 + sqrt(1.0 - fdr_target * fmin(1.0, avg_q)));
    double pi_thr = fmax(0.5, fmin(factor, 0.95));

    return pi_thr;
}

/* --------------- Stable Graph Extraction --------------- */

CdfGraph* cdf_bootstrap_stable_graph(const CdfBootstrapEdge *edges,
                                      int n_edges, int p, double thresh)
{
    if (!edges || n_edges < 1 || p < 1) return NULL;

    CdfGraph *g = cdf_graph_create(p);
    if (!g) return NULL;

    for (int k = 0; k < n_edges; k++) {
        if (edges[k].freq_stable < thresh) continue;

        int u = edges[k].u, v = edges[k].v;
        int uv = edges[k].n_directed_uv;
        int vu = edges[k].n_directed_vu;
        int un = edges[k].n_undirected;
        int bi = edges[k].n_bidirected;

        if (uv >= vu && uv >= un && uv >= bi && uv > 0)
            cdf_graph_add_edge(g, u, v, CDF_EDGE_DIRECTED);
        else if (vu >= uv && vu >= un && vu >= bi && vu > 0)
            cdf_graph_add_edge(g, v, u, CDF_EDGE_DIRECTED);
        else if (bi >= uv && bi >= vu && bi >= un && bi > 0)
            cdf_graph_add_edge(g, u, v, CDF_EDGE_BIDIRECTED);
        else if (un > 0)
            cdf_graph_add_edge(g, u, v, CDF_EDGE_UNDIRECTED);
    }

    return g;
}
