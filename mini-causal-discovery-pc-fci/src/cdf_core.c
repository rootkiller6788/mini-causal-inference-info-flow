/* cdf_core.c — Core data structures for PC/FCI causal discovery
 *
 * Implements graph representation, dataset management, and basic
 * operations needed by the PC and FCI algorithms.
 */

#include "cdf_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ──────────────── Dataset Lifecycle ──────────────────────────────── */

CdfDataset* cdf_dataset_create(const double *data, int N, int p)
{
    if (!data || N < 2 || p < 1 || p > CDF_MAX_NODES || N > CDF_MAX_SAMPLES)
        return NULL;

    CdfDataset *ds = (CdfDataset*)calloc(1, sizeof(CdfDataset));
    if (!ds) return NULL;

    size_t n_elem = (size_t)p * (size_t)N;
    ds->data = (double*)malloc(n_elem * sizeof(double));
    ds->var_names = (char**)calloc((size_t)p, sizeof(char*));

    if (!ds->data) { free(ds->var_names); free(ds); return NULL; }

    memcpy(ds->data, data, n_elem * sizeof(double));
    ds->N = N;
    ds->p = p;

    /* Default variable names */
    char buf[16];
    for (int i = 0; i < p && ds->var_names; i++) {
        snprintf(buf, sizeof(buf), "X%d", i + 1);
        ds->var_names[i] = strdup(buf);
    }

    return ds;
}

void cdf_dataset_free(CdfDataset *ds)
{
    if (!ds) return;
    if (ds->var_names) {
        for (int i = 0; i < ds->p; i++)
            free(ds->var_names[i]);
        free(ds->var_names);
    }
    free(ds->data);
    free(ds);
}

const char* cdf_dataset_var_name(const CdfDataset *ds, int i)
{
    if (!ds || !ds->var_names || i < 0 || i >= ds->p) return "?";
    return ds->var_names[i];
}

/* ──────────────── Graph Lifecycle ────────────────────────────────── */

CdfGraph* cdf_graph_create(int p)
{
    if (p < 1 || p > CDF_MAX_NODES) return NULL;

    CdfGraph *g = (CdfGraph*)calloc(1, sizeof(CdfGraph));
    if (!g) return NULL;

    int p2 = p * p;
    g->edges        = (CdfEdgeType*)calloc((size_t)p2, sizeof(CdfEdgeType));
    g->sepsets      = (CdfSepSet*)calloc((size_t)p2, sizeof(CdfSepSet));
    g->sepset_count = (int*)malloc((size_t)p2 * sizeof(int));

    if (!g->edges || !g->sepsets || !g->sepset_count) {
        cdf_graph_free(g); return NULL;
    }

    g->p = p;

    /* Initialize sepset_count to -1 (no sepset) */
    for (int i = 0; i < p2; i++)
        g->sepset_count[i] = -1;

    /* Initialize as fully disconnected */
    for (int i = 0; i < p2; i++)
        g->edges[i] = CDF_EDGE_NONE;

    return g;
}

void cdf_graph_free(CdfGraph *g)
{
    if (!g) return;
    /* Free sepset arrays */
    if (g->sepsets) {
        int p2 = g->p * g->p;
        for (int i = 0; i < p2; i++) {
            free(g->sepsets[i].set);
        }
        free(g->sepsets);
    }
    free(g->sepset_count);
    free(g->edges);
    free(g);
}

void cdf_graph_add_edge(CdfGraph *g, int u, int v, CdfEdgeType type)
{
    if (!g || u < 0 || v < 0 || u >= g->p || v >= g->p || u == v) return;

    int p = g->p;
    g->edges[u * p + v] = type;

    /* For undirected/bidirected/nondirected edges, mirror */
    if (type == CDF_EDGE_UNDIRECTED || type == CDF_EDGE_NONDIR ||
        type == CDF_EDGE_BIDIRECTED)
        g->edges[v * p + u] = type;
}

void cdf_graph_remove_edge(CdfGraph *g, int u, int v)
{
    if (!g || u < 0 || v < 0 || u >= g->p || v >= g->p) return;

    int p = g->p;
    g->edges[u * p + v] = CDF_EDGE_NONE;
    g->edges[v * p + u] = CDF_EDGE_NONE;
}

bool cdf_graph_has_edge(const CdfGraph *g, int u, int v)
{
    if (!g || u < 0 || v < 0 || u >= g->p || v >= g->p) return false;
    return g->edges[u * g->p + v] != CDF_EDGE_NONE;
}

CdfEdgeType cdf_graph_edge_type(const CdfGraph *g, int u, int v)
{
    if (!g || u < 0 || v < 0 || u >= g->p || v >= g->p) return CDF_EDGE_NONE;
    return g->edges[u * g->p + v];
}

/* ──────────────── Reachability (Directed Path) ───────────────────── */

bool cdf_graph_reachable(const CdfGraph *g, int u, int v)
{
    if (!g || u < 0 || v < 0 || u >= g->p || v >= g->p) return false;
    if (u == v) return true;

    int p = g->p;
    bool *visited = (bool*)calloc((size_t)p, sizeof(bool));
    int *queue = (int*)malloc((size_t)p * sizeof(int));
    if (!visited || !queue) { free(visited); free(queue); return false; }

    int head = 0, tail = 0;
    queue[tail++] = u;
    visited[u] = true;

    while (head < tail) {
        int cur = queue[head++];
        if (cur == v) { free(visited); free(queue); return true; }

        for (int w = 0; w < p; w++) {
            if (!visited[w] && g->edges[cur * p + w] == CDF_EDGE_DIRECTED) {
                visited[w] = true;
                queue[tail++] = w;
            }
        }
    }

    free(visited); free(queue);
    return false;
}

/* ──────────────── Configuration Defaults ─────────────────────────── */

CdfPCConfig cdf_pc_config_default(void)
{
    CdfPCConfig cfg;
    cfg.alpha = CDF_ALPHA_DEFAULT;
    cfg.max_cond_size = -1;
    cfg.use_fisher_z = true;
    cfg.conservative = false;
    cfg.verbose = false;
    return cfg;
}

CdfFCIConfig cdf_fci_config_default(void)
{
    CdfFCIConfig cfg;
    cfg.alpha = CDF_ALPHA_DEFAULT;
    cfg.max_cond_size = -1;
    cfg.max_path_length = 4;
    cfg.use_fisher_z = true;
    cfg.conservative = false;
    cfg.complete_rule_set = true;
    cfg.verbose = false;
    return cfg;
}

int cdf_version(void)
{
    return 100;  /* version 1.0.0 */
}