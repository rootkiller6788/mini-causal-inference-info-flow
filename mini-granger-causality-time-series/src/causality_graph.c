/* causality_graph.c -- Directed causality graph from pairwise Granger tests.
 *
 * Builds a weighted directed graph where edges represent Granger-causal
 * relationships between multiple time series. This generalizes bivariate
 * Granger causality to the multivariate network setting.
 *
 * Knowledge points:
 * L2: Granger causality graph as causal network representation
 * L5: Kahn topological sort for DAG detection
 * L5: Brandes algorithm for betweenness centrality
 * L6: Multivariate causal discovery via pairwise Granger
 * L8: Hub detection in causal networks
 */

#include "causality_graph.h"
#include "granger_test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

CausalityGraph* graph_build(const TimeSeries* ts, int max_lag, double alpha) {
    if (!ts || ts->dim < 2 || max_lag < 1) return NULL;
    int M = ts->dim, N = ts->length;
    CausalityGraph* g = calloc(1, sizeof(CausalityGraph));
    if (!g) return NULL;
    g->n_vars = M;
    g->n_lags = max_lag;
    g->significance_level = alpha;
    g->adjacency = malloc((size_t)M * sizeof(bool*));
    g->edge_weights = malloc((size_t)M * sizeof(double*));
    g->f_statistics = malloc((size_t)M * sizeof(double*));
    g->p_values = malloc((size_t)M * sizeof(double*));
    for (int i = 0; i < M; i++) {
        g->adjacency[i] = calloc((size_t)M, sizeof(bool));
        g->edge_weights[i] = calloc((size_t)M, sizeof(double));
        g->f_statistics[i] = calloc((size_t)M, sizeof(double));
        g->p_values[i] = calloc((size_t)M, sizeof(double));
    }
    /* Extract each dimension */
    double** series = malloc((size_t)M * sizeof(double*));
    for (int d = 0; d < M; d++) {
        series[d] = malloc((size_t)N * sizeof(double));
        for (int t = 0; t < N; t++)
            series[d][t] = ts->values[t * M + d];
    }
    /* Perform M*(M-1) pairwise Granger tests */
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            if (i == j) continue;
            int best_lag = gt_optimal_lag(series[i], series[j], N, max_lag);
            GrangerTestResult* gt = gt_create();
            gt_test(gt, series[i], series[j], N, best_lag, alpha);
            g->f_statistics[i][j] = gt->f_statistic;
            g->p_values[i][j] = gt->p_value;
            if (gt->is_significant) {
                g->adjacency[i][j] = true;
                g->edge_weights[i][j] = (gt->p_value > 0) ? -log10(gt->p_value) : 16.0;
                g->total_edges++;
            }
            gt_free(gt);
        }
    }
    for (int d = 0; d < M; d++) free(series[d]);
    free(series);
    return g;
}

void graph_free(CausalityGraph* g) {
    if (!g) return;
    for (int i = 0; i < g->n_vars; i++) {
        free(g->adjacency[i]);
        free(g->edge_weights[i]);
        free(g->f_statistics[i]);
        free(g->p_values[i]);
    }
    free(g->adjacency); free(g->edge_weights);
    free(g->f_statistics); free(g->p_values);
    free(g);
}

void graph_print(const CausalityGraph* g) {
    if (!g) return;
    printf("Causality Graph: %d vars, %d edges (alpha=%.3f)\n",
        g->n_vars, g->total_edges, g->significance_level);
    for (int i = 0; i < g->n_vars; i++)
        for (int j = 0; j < g->n_vars; j++)
            if (g->adjacency[i][j])
                printf("  %d -> %d: F=%.2f p=%.4e w=%.2f\n",
                    i, j, g->f_statistics[i][j], g->p_values[i][j], g->edge_weights[i][j]);
}

/* Betweenness centrality via Brandes (2001) algorithm.
 * C_B(v) = sum_{s!=v!=t} sigma_{st}(v) / sigma_{st}
 * where sigma_{st} is number of shortest paths s->t and
 * sigma_{st}(v) counts those passing through v.
 * Complexity: O(V * E) for unweighted graphs. */
GraphMetrics* graph_metrics(const CausalityGraph* g) {
    if (!g) return NULL;
    int n = g->n_vars;
    GraphMetrics* m = calloc(1, sizeof(GraphMetrics));
    m->n_vars = n;
    m->in_degree = calloc((size_t)n, sizeof(int));
    m->out_degree = calloc((size_t)n, sizeof(int));
    m->betweenness = calloc((size_t)n, sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (g->adjacency[i][j]) {
                m->out_degree[i]++;
                m->in_degree[j]++;
            }
    /* Brandes betweenness centrality using BFS from each source */
    for (int s = 0; s < n; s++) {
        int* dist = malloc((size_t)n * sizeof(int));
        double* sigma = calloc((size_t)n, sizeof(double));
        double* delta = calloc((size_t)n, sizeof(double));
        int* stack = malloc((size_t)n * sizeof(int));
        int** predecessors = malloc((size_t)n * sizeof(int*));
        for (int i = 0; i < n; i++) {
            dist[i] = -1;
            predecessors[i] = NULL;
        }
        dist[s] = 0;
        sigma[s] = 1.0;
        int* queue = malloc((size_t)n * sizeof(int));
        int qh = 0, qt = 0, sp = 0;
        queue[qt++] = s;
        while (qh < qt) {
            int v = queue[qh++];
            stack[sp++] = v;
            for (int w = 0; w < n; w++) {
                if (!g->adjacency[v][w]) continue;
                if (dist[w] < 0) {
                    dist[w] = dist[v] + 1;
                    queue[qt++] = w;
                }
                if (dist[w] == dist[v] + 1) {
                    sigma[w] += sigma[v];
                    int* new_pred = realloc(predecessors[w], (size_t)(1 + (predecessors[w] ? 1 : 0)) * sizeof(int));
                    if (predecessors[w]) {
                        new_pred[1] = v;
                        predecessors[w] = new_pred;
                    } else {
                        new_pred[0] = v;
                        predecessors[w] = new_pred;
                    }
                }
            }
        }
        /* Back-propagation */
        while (sp > 0) {
            int w = stack[--sp];
            if (predecessors[w]) {
                double coeff = (1.0 + delta[w]) / sigma[w];
                int* preds = predecessors[w];
                int v = preds[0];
                delta[v] += sigma[v] * coeff;
            }
            if (w != s) m->betweenness[w] += delta[w];
        }
        free(dist); free(sigma); free(delta); free(stack); free(queue);
        for (int i = 0; i < n; i++) free(predecessors[i]);
        free(predecessors);
    }
    /* Normalize */
    double norm = (n > 2) ? 2.0 / ((n - 1.0) * (n - 2.0)) : 1.0;
    for (int i = 0; i < n; i++) m->betweenness[i] *= norm;
    return m;
}

void gm_free(GraphMetrics* m) {
    if (!m) return;
    free(m->in_degree); free(m->out_degree);
    free(m->betweenness); free(m);
}

void gm_print(const GraphMetrics* m) {
    if (!m) return;
    printf("Graph Metrics (%d vars):\n", m->n_vars);
    for (int i = 0; i < m->n_vars; i++)
        printf("  Node %d: in=%d out=%d bc=%.4f\n",
            i, m->in_degree[i], m->out_degree[i], m->betweenness[i]);
}

/* Find causal hubs: nodes with high out-degree relative to network size.
 * A hub is a node that Granger-causes many others. */
int* graph_find_hubs(const CausalityGraph* g, double tau_fraction) {
    if (!g || tau_fraction < 0.0 || tau_fraction > 1.0) return NULL;
    int n = g->n_vars;
    int threshold = (int)(tau_fraction * (n - 1));
    if (threshold < 1) threshold = 1;
    int* hubs = malloc((size_t)(n + 1) * sizeof(int));
    int nh = 0;
    for (int i = 0; i < n; i++) {
        int od = 0;
        for (int j = 0; j < n; j++) if (g->adjacency[i][j]) od++;
        if (od >= threshold) hubs[nh++] = i;
    }
    hubs[nh] = -1;
    return hubs;
}

/* Check if causality graph is a DAG using Kahn's algorithm.
 * Maintains in-degree counter and repeatedly removes nodes with zero in-degree.
 * If all nodes removed, the graph has no cycles. O(V+E). */
bool graph_is_dag(const CausalityGraph* g) {
    if (!g) return true;
    int n = g->n_vars;
    int* indeg = calloc((size_t)n, sizeof(int));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (g->adjacency[i][j]) indeg[j]++;
    int* queue = malloc((size_t)n * sizeof(int));
    int qh = 0, qt = 0, cnt = 0;
    for (int i = 0; i < n; i++)
        if (indeg[i] == 0) queue[qt++] = i;
    while (qh < qt) {
        int v = queue[qh++]; cnt++;
        for (int w = 0; w < n; w++) {
            if (g->adjacency[v][w]) {
                indeg[w]--;
                if (indeg[w] == 0) queue[qt++] = w;
            }
        }
    }
    free(indeg); free(queue);
    return cnt == n;
}