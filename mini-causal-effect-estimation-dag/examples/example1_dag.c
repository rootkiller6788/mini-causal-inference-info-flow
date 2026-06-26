#include <stdio.h>
#include <stdlib.h>
#include "dag_representation.h"
#include "causal_identification.h"

int main(void) {
    printf("=== Example 1: DAG Construction and Back-Door Criterion ===\n");
    printf("Causal DAG: X0 -> X1, X0 -> X2, X1 -> X3, X2 -> X3\n");
    printf("Question: is the effect of X1 on X3 identifiable?\n\n");

    /* Build DAG: 4 variables, X0 is common cause of X1 and X2,
       both affect X3. This creates a fork structure:
       X1 <-- X0 --> X2 --> X3, with X1 --> X3. */
    DAG* dag = dag_create(4);
    dag_add_edge(dag, 0, 1);  /* X0 -> X1 */
    dag_add_edge(dag, 0, 2);  /* X0 -> X2 */
    dag_add_edge(dag, 1, 3);  /* X1 -> X3 */
    dag_add_edge(dag, 2, 3);  /* X2 -> X3 */

    /* Verify DAG structure */
    printf("Node count: %d\n", dag->n);
    printf("Edge X0->X1: %s\n", dag_has_edge(dag, 0, 1) ? "yes" : "no");
    printf("Edge X1->X0: %s\n", dag_has_edge(dag, 1, 0) ? "yes" : "no");
    printf("Is acyclic: %s\n\n", dag_is_acyclic(dag) ? "yes" : "no");

    /* Topological sort: verify causal ordering */
    dag_topological_sort(dag);
    printf("Causal order (topological):");
    for (int i = 0; i < dag->n; i++) printf(" X%d", dag->topo_order[i]);
    printf("\n\n");

    /* Find back-door paths from X1 to X3:
       Back-door path: X1 <- X0 -> X2 -> X3 */
    Path* bd_paths[10];
    int n_bd = dag_find_backdoor_paths(dag, 1, 3, bd_paths, 10);
    printf("Back-door paths X1->X3: %d found\n", n_bd);
    for (int p = 0; p < n_bd; p++) {
        printf("  Path %d:", p);
        for (int k = 0; k < bd_paths[p]->length; k++)
            printf(" X%d", bd_paths[p]->nodes[k]);
        printf("\n");
        free(bd_paths[p]->nodes);
        free(bd_paths[p]);
    }
    printf("\n");

    /* Back-door adjustment: X0 blocks the back-door path */
    int adj[10];
    int na = causal_find_backdoor_adjustment(dag, 1, 3, adj, 10);
    printf("Back-door adjustment set: %d variables\n", na);
    printf("Variables to adjust for:");
    for (int i = 0; i < na; i++) printf(" X%d", adj[i]);
    printf("\n");
    printf("Effect of X1 on X3: %s\n\n",
           na > 0 ? "IDENTIFIABLE via back-door" : "NOT identifiable");

    /* Full identification result */
    IdentificationResult* ir = causal_identify(dag, 1, 3);
    printf("Identification result:\n");
    printf("  Identifiable: %s\n", ir->identifiable ? "YES" : "NO");
    printf("  Method: %d (1=Backdoor, 2=Frontdoor, 3=IV, 4=Do-Calc)\n",
           ir->method);
    if (ir->explanation) printf("  %s\n", ir->explanation);
    ident_result_free(ir);

    dag_free(dag);
    printf("\n=== Example 1 Complete ===\n");
    return 0;
}
