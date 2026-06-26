/* granger_app_neuroscience.c -- L7 Application: Neural Granger causality.
 *
 * Knowledge points:
 * L7: Neuroscience application - fMRI/EEG effective connectivity
 *
 * Granger causality is widely used in neuroscience to infer effective
 * connectivity from fMRI BOLD signals, EEG/MEG recordings, and
 * multi-unit neural spike trains.
 *
 * This implements a simulated 5-region brain network with known
 * ground-truth connectivity:
 *   Visual(0) -> Parietal(1) [dorsal stream]
 *   Visual(0) -> Temporal(2) [ventral stream]
 *   Parietal(1) -> Frontal(3) [attention]
 *   Temporal(2) -> Frontal(3) [recognition]
 *   Frontal(3) -> Motor(4) [action]
 *
 * Reference: Seth, Barrett, Barnett (2015), J. Neurosci. 35(8):3293-3297.
 *            Friston et al. (2013), NeuroImage, 80:144-168.
 */

#include "ts_core.h"
#include "granger_test.h"
#include "causality_graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Simulate 5-region brain network with known directed connectivity.
 * Hemodynamic response is approximated by AR(2) filtering of neural activity.
 * Each region receives input from its upstream regions plus noise. */
static TimeSeries* simulate_brain_network(int length) {
    int n_regions = 5;
    TimeSeries* ts = ts_create(length, n_regions);
    if (!ts) return NULL;
    /* Initial conditions (resting state) */
    for (int r = 0; r < n_regions; r++)
        ts->values[r] = 0.1 * ((double)rand() / RAND_MAX - 0.5);
    /* Connectivity matrix (ground truth):
     *       V  P  T  F  M
     * V [ .9 0  0  0  0 ]  visual: strong AR
     * P [ .3 .8 0  0  0 ]  parietal: V->P
     * T [ .3 0 .85 0  0 ]  temporal: V->T
     * F [ 0 .2 .2 .8 0 ]  frontal: P->F, T->F
     * M [ 0  0  0 .3 .9 ]  motor: F->M
     * coupling[i][j] = influence from j to i
     */
    double coupling[5][5] = {
        {0.90, 0.00, 0.00, 0.00, 0.00},
        {0.30, 0.80, 0.00, 0.00, 0.00},
        {0.30, 0.00, 0.85, 0.00, 0.00},
        {0.00, 0.20, 0.20, 0.80, 0.00},
        {0.00, 0.00, 0.00, 0.30, 0.90}
    };
    double noise_std[5] = {0.5, 0.3, 0.3, 0.25, 0.2};
    for (int t = 1; t < length; t++) {
        for (int r = 0; r < n_regions; r++) {
            double val = 0.0;
            for (int src = 0; src < n_regions; src++)
                val += coupling[r][src] * ts->values[(t - 1) * n_regions + src];
            val += noise_std[r] * ((double)rand() / RAND_MAX - 0.5);
            ts->values[t * n_regions + r] = val;
        }
    }
    return ts;
}

int main(void) {
    printf("=== Granger Causality: Neuroscience Application ===\n\n");
    srand(12345);
    int T = 300;
    const char* region_names[] = {"Visual","Parietal","Temporal","Frontal","Motor"};
    TimeSeries* ts = simulate_brain_network(T);
    printf("Simulated 5-region brain network (%d time points).\n", T);
    printf("Regions: Visual, Parietal, Temporal, Frontal, Motor.\n\n");
    /* Extract all series */
    double** signals = malloc(5 * sizeof(double*));
    for (int r = 0; r < 5; r++)
        signals[r] = ts_extract_dim(ts, r);
    /* Build causality graph */
    printf("--- Causality Graph Reconstruction ---\n");
    CausalityGraph* cg = graph_build(ts, 5, 0.05);
    if (cg) {
        graph_print(cg);
        printf("\nGround truth edges: V->P, V->T, P->F, T->F, F->M\n\n");
        /* Graph metrics */
        GraphMetrics* gm = graph_metrics(cg);
        if (gm) {
            gm_print(gm);
            printf("\nInterpretation:\n");
            printf("  Visual cortex is a causal source (high out-degree).\n");
            printf("  Frontal cortex integrates information (high in-degree).\n");
            printf("  Betweenness centrality identifies relay regions.\n");
            gm_free(gm);
        }
        /* Hub detection */
        int* hubs = graph_find_hubs(cg, 0.3);
        if (hubs) {
            printf("\nCausal hubs (tau=0.3): ");
            for (int i = 0; hubs[i] >= 0; i++)
                printf("%s ", region_names[hubs[i]]);
            printf("\n");
            free(hubs);
        }
        /* DAG check */
        printf("Is DAG? %s\n", graph_is_dag(cg) ? "Yes" : "No");
        graph_free(cg);
    }
    /* Detailed pairwise tests */
    printf("\n--- Detailed Pairwise Tests ---\n");
    GrangerTestResult* gt = gt_create();
    int test_pairs[6][2] = {{0,1},{0,2},{1,3},{2,3},{3,4},{2,1}};
    for (int tp = 0; tp < 6; tp++) {
        int src = test_pairs[tp][0], dst = test_pairs[tp][1];
        gt_test(gt, signals[src], signals[dst], T, 4, 0.05);
        printf("%-10s -> %-10s: F=%.2f p=%.4f %s\n",
            region_names[src], region_names[dst],
            gt->f_statistic, gt->p_value,
            gt->is_significant ? "***" : "");
    }
    gt_free(gt);
    for (int r = 0; r < 5; r++) free(signals[r]);
    free(signals);
    ts_free(ts);
    printf("\n=== Application Complete ===\n");
    return 0;
}