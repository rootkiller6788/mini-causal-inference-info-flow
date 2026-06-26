#include "scm_core.h"
#include "scm_graph.h"
#include <stdio.h>
int main(void) {
    printf("========================================\n");
    printf("  Pearl SCM - Demo 1: Basic SCM\n");
    printf("========================================\n\n");
    SCM_Model* m=scm_create();
    int z=scm_add_variable(m,"Z",SCM_EXOGENOUS);
    int x=scm_add_variable(m,"X",SCM_ENDOGENOUS);
    int y=scm_add_variable(m,"Y",SCM_ENDOGENOUS);
    scm_add_edge(m,z,x,0.5); scm_add_edge(m,z,y,0.3); scm_add_edge(m,x,y,0.7);
    m->vars[z].noise=1.0; m->vars[x].noise=0.2; m->vars[y].noise=0.1;
    scm_set_equation(m,x,scm_linear_eq); scm_set_equation(m,y,scm_linear_eq);
    scm_topological_sort(m); scm_sample(m);
    scm_print(m);
    printf("\nGraph: Z->X, Z->Y, X->Y\n");
    printf("DAG: %s\n", scm_is_dag(m)?"YES":"NO");
    SCM_VarSet pa; scm_get_parents(m,y,&pa);
    printf("Parents of Y: %d vars\n", pa.n);
    SCM_VarSet Zset; scm_varset_init(&Zset); scm_varset_add(&Zset,z);
    printf("X d-sep Y | {Z}: %s\n", scm_d_separated(m,x,y,&Zset)?"YES":"NO");
    /* Markov blanket: parents + children + children's other parents */
    SCM_VarSet mb; scm_markov_blanket(m, x, &mb);
    printf("Markov blanket of X: {");
    for (int i = 0; i < mb.n; i++) printf("%s%d", i?",":"", mb.nodes[i]);
    printf("} (%d vars)\n", mb.n);
    /* Find all colliders in the graph */
    SCM_Collider colliders[16]; int nc;
    scm_find_colliders(m, colliders, &nc);
    printf("Colliders in graph: %d\n", nc);
    for (int i = 0; i < nc; i++)
        printf("  %d->%d<-%d\n", colliders[i].parent1, colliders[i].child, colliders[i].parent2);
    /* Ancestors and descendants */
    SCM_VarSet an; scm_get_ancestors(m, y, &an);
    printf("Ancestors of Y: %d vars\n", an.n);
    /* Monte Carlo sampling demo */
    printf("\n--- Monte Carlo (5 samples) ---\n");
    double data[15];
    scm_sample_n(m, data, 5);
    for (int i = 0; i < 5; i++)
        printf("  sample %d: X=%.3f Y=%.3f\n", i+1, data[i*3+1], data[i*3+2]);
    /* Verify it's still a DAG */
    printf("\nDAG integrity check: %s\n", scm_is_dag(m)?"PASS":"FAIL");
    scm_free(m);
    printf("\n========================================\n");
    return 0;
}
/*
 * Pearl SCM Demo 1: Basic Structural Causal Model
 * Demonstrates: SCM creation, graph structure, d-separation,
 *   Markov blanket, colliders, ancestors/descendants, sampling.
 * References:
 *   Pearl (2009) Causality: Models, Reasoning, and Inference, 2nd ed.
 *   Morgan & Winship (2015) Counterfactuals and Causal Inference, 2nd ed.
 *   Pearl, Glymour & Jewell (2016) Causal Inference in Statistics: A Primer.
 */
