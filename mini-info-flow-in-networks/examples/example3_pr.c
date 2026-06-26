/* Example: Network Dynamics — PageRank, Centrality, and SIR Model.
 * Demonstrates information flow measures on directed graphs:
 * PageRank for information ranking, SIR epidemic spread,
 * and network entropy measures.
 * Ref: Page et al. (1999); Kermack & McKendrick (1927) */
#include <stdio.h>
#include "ifn_core.h"
#include "ifn_causal.h"
#include "ifn_dynamic.h"
#include "ifn_network.h"
int main(void) {
    printf("=== Network Dynamics Demo ===\n");
    /* Build web-like graph with cycle 0→1→2→0 and dangling node 3 */
    IFN_CausalGraph* g=ifn_causal_create(4);
    ifn_causal_add_edge(g,0,1,1.0); ifn_causal_add_edge(g,1,2,1.0);
    ifn_causal_add_edge(g,2,0,1.0); ifn_causal_add_edge(g,2,3,1.0);

    /* PageRank */
    double pr[4]; ifn_pagerank(g,0.85,50,1e-6,pr);
    printf("PageRank:\n");
    for(int i=0;i<4;i++) printf("  PR[%d]=%.4f\n",i,pr[i]);

    /* Network entropy */
    double H=ifn_network_entropy(g);
    printf("Network entropy H=%.4f\n", H);

    /* Epidemic threshold */
    double tau=ifn_epidemic_threshold(g);
    printf("Epidemic threshold tau_c=%.4f\n", tau);

    /* Clustering coefficient */
    double cc=ifn_clustering_coefficient(g);
    printf("Clustering coefficient CC=%.4f\n", cc);

    /* SIR model simulation */
    int S[10]={4,0}, I[10]={0,0}, R[10]={0,0};
    ifn_sir_model(g,0.3,0.1,10,S,I,R);
    printf("SIR after 10 steps: S=%d I=%d R=%d\n", S[9], I[9], R[9]);

    /* ER random graph generation */
    IFN_CausalGraph* er=ifn_network_random_erdos(8,0.3);
    printf("ER graph: n=%d  edges=%d\n", er->n_nodes, er->n_edges);
    double er_cc=ifn_clustering_coefficient(er);
    printf("ER clustering=%.4f\n", er_cc);
    ifn_causal_free(er);

    ifn_causal_free(g);
    printf("Example 3 PASSED\n");
    return 0;
}
