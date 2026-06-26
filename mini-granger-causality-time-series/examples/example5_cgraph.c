/* Example 5: Causality Graph Discovery.
 * Demonstrates building a directed causal graph from multiple
 * time series using pairwise Granger causality tests. */
#include "ts_core.h"
#include "causality_graph.h"
#include <stdio.h>
#include <stdlib.h>
int main(void) {
    printf("=== Causality Graph Discovery ===\n");
    srand(42);
    int T = 200, M = 4;
    /* Simulate 4-variable system with known structure:
     * 0 -> 1, 0 -> 2, 0 -> 3 (a hub causing all others)
     * 1 -> 2 (a chain link)
     */
    TimeSeries* ts = ts_create(T, M);
    ts->values[0]=1.0; ts->values[1]=0.5; ts->values[2]=0.0; ts->values[3]=0.0;
    for (int t = 1; t < T; t++) {
        double v0 = ts->values[(t-1)*M+0], v1 = ts->values[(t-1)*M+1];
        double v2 = ts->values[(t-1)*M+2], v3 = ts->values[(t-1)*M+3];
        ts->values[t*M+0] = 0.8*v0 + 0.1*((double)rand()/RAND_MAX-0.5);
        ts->values[t*M+1] = 0.3*v0 + 0.6*v1 + 0.1*((double)rand()/RAND_MAX-0.5);
        ts->values[t*M+2] = 0.25*v0 + 0.1*v1 + 0.5*v2 + 0.1*((double)rand()/RAND_MAX-0.5);
        ts->values[t*M+3] = 0.2*v0 + 0.7*v3 + 0.1*((double)rand()/RAND_MAX-0.5);
    }
    printf("Ground truth: 0->1, 0->2, 0->3, 1->2\n\n");
    CausalityGraph* cg = graph_build(ts, 4, 0.05);
    if (cg) {
        graph_print(cg);
        GraphMetrics* gm = graph_metrics(cg);
        if (gm) { gm_print(gm); gm_free(gm); }
        int* hubs = graph_find_hubs(cg, 0.4);
        printf("\nCausal hubs: ");
        for (int i = 0; hubs && hubs[i] >= 0; i++) printf("%d ", hubs[i]);
        printf("\nIs DAG: %s\n", graph_is_dag(cg) ? "Yes" : "No");
        if (hubs) free(hubs);
        graph_free(cg);
    }
    ts_free(ts);
    return 0;
}