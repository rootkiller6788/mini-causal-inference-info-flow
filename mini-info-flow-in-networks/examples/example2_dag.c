/* Example: Causal Graphs & d-Separation on DAGs.
 * Demonstrates DAG construction, edge manipulation, d-separation testing,
 * and causal flow computation.
 * Ref: Pearl (2009) Causality, 2nd ed. */
#include <stdio.h>
#include "ifn_core.h"
#include "ifn_causal.h"
int main(void) {
    printf("=== Causal Graph Demo ===\n");
    /* Build chain: 0→1→2, with fork 0→3 */
    IFN_CausalGraph* g=ifn_causal_create(4);
    ifn_causal_add_edge(g,0,1,0.8); ifn_causal_add_edge(g,1,2,0.5);
    ifn_causal_add_edge(g,0,3,0.3);
    printf("DAG=%s  edges=%d\n", ifn_is_dag(g)?"YES":"NO", g->n_edges);
    /* d-separation: 0 and 2 are d-separated by {1} */
    int z[]={1};
    printf("d-sep(0,2 | {1})=%s\n", ifn_d_separated(g,0,2,z,1)?"YES":"NO");
    /* d-separation without conditioning: 0 and 2 connected */
    printf("d-sep(0,2 | {})=%s\n", ifn_d_separated(g,0,2,NULL,0)?"YES":"NO");
    /* Causal flow from 0 to 2 */
    double flow=ifn_causal_flow(g,0,2);
    printf("causal_flow(0->2)=%.4f\n", flow);
    /* Topological sort */
    int order[4];
    ifn_topological_sort(g,order);
    printf("Topological order: ");
    for(int i=0;i<4;i++) printf("%d ",order[i]);
    printf("\n");

    ifn_causal_free(g);
    printf("Example 2 PASSED\n");
    return 0;
}
