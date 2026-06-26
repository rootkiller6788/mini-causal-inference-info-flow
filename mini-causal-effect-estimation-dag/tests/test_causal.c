#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "dag_representation.h"
#include "scm_model.h"
#include "causal_identification.h"
#include "adjustment_sets.h"
#include "effect_estimators.h"

int main(void) {
    printf("=== Causal Effect Estimation Test Suite ===\n");
    /* Test 1: DAG creation */
    printf("  Test 1: DAG creation...\n");
    DAG* dag = dag_create(5);
    assert(dag); assert(dag->n == 5);

    /* Test 2: Add edges */
    printf("  Test 2: Edge addition...\n");
    assert(dag_add_edge(dag, 0, 1) == 0);
    assert(dag_add_edge(dag, 0, 2) == 0);
    assert(dag_add_edge(dag, 1, 3) == 0);
    assert(dag_add_edge(dag, 2, 3) == 0);
    assert(dag_add_edge(dag, 3, 4) == 0);
    assert(dag_has_edge(dag, 0, 1));
    assert(!dag_has_edge(dag, 1, 0));

    /* Test 3: Cycle rejection */
    printf("  Test 3: Cycle detection...\n");
    assert(dag_add_edge(dag, 3, 0) == -1);
    assert(dag_is_acyclic(dag));

    /* Test 4: Topological sort */
    printf("  Test 4: Topological sort...\n");
    assert(dag_topological_sort(dag) == 0);
    assert(dag->topo_order[0] == 0);

    /* Test 5: Path finding */
    printf("  Test 5: Path operations...\n");
    assert(dag_has_path(dag, 0, 4));
    assert(!dag_has_path(dag, 4, 0));

    /* Test 6: Ancestors and descendants */
    printf("  Test 6: Ancestors/descendants...\n");
    int anc[5]; int na = dag_ancestors(dag, 3, anc, 5);
    assert(na >= 2);
    int desc[5]; int nd = dag_descendants(dag, 0, desc, 5);
    assert(nd >= 3);

    /* Test 7: Back-door paths */
    printf("  Test 7: Back-door paths...\n");
    Path* bd[20]; int nb = dag_find_backdoor_paths(dag, 1, 3, bd, 20);
    (void)nb;

    /* Test 8: SCM creation */
    printf("  Test 8: SCM lifecycle...\n");
    SCModel* scm = scm_create(3);
    assert(scm); assert(scm->n_vars == 3);
    dag_add_edge(scm->dag, 0, 1);
    dag_add_edge(scm->dag, 0, 2);
    dag_add_edge(scm->dag, 1, 2);
    scm->noise[0] = 1.0; scm->noise[1] = 0.0; scm->noise[2] = 0.0;
    scm_set_equation(scm, 0, scm_linear_eq, 0);
    scm_set_equation(scm, 1, scm_linear_eq, 1);
    scm_set_equation(scm, 2, scm_linear_eq, 2);
    assert(scm_simulate(scm) == 0);

    /* Test 9: Do-intervention */
    printf("  Test 9: Do-intervention...\n");
    SCModel* scm_do = scm_do_intervention(scm, 1, 10.0);
    assert(scm_do); scm_simulate(scm_do);
    assert(fabs(scm_do->values[1] - 10.0) < 0.01);
    scm_free(scm_do);

    /* Test 10: Average Causal Effect */
    printf("  Test 10: ACE estimation...\n");
    double ace = scm_average_causal_effect(scm, 0, 1.0, 0.0, 2, 100);
    (void)ace;

    /* Test 11: Back-door criterion */
    printf("  Test 11: Back-door criterion...\n");
    int adj2[5]; int na2 = causal_find_backdoor_adjustment(dag, 1, 3, adj2, 5);
    (void)na2;

    /* Test 12: Causal identification */
    printf("  Test 12: Causal identification...\n");
    IdentificationResult* ir = causal_identify(dag, 1, 3);
    assert(ir); ident_result_free(ir);

    /* Test 13: Adjustment sets */
    printf("  Test 13: Adjustment sets...\n");
    AdjustmentSet* as = adjustment_parent_set(dag, 1);
    assert(as); assert(as->n_vars >= 1);
    adjustment_set_free(as);

    /* Test 14: Observational data */
    printf("  Test 14: Observational data...\n");
    ObservationalData* od = obs_data_create(10, 2);
    assert(od);
    for (int i = 0; i < 10; i++) {
        od->X[i*2]=0.5; od->X[i*2+1]=0.3;
        od->T[i]=i%2; od->Y[i]=(double)(i%2)*2.0+1.0;
    }

    /* Test 15: Propensity score */
    printf("  Test 15: Propensity score...\n");
    double* ps = propensity_score_estimate(od);
    assert(ps);

    /* Test 16: IPW */
    printf("  Test 16: IPW estimator...\n");
    CausalEffect* ipw = ipw_estimate(od, ps);
    assert(ipw);
    causal_effect_free(ipw);

    /* Test 17: G-Computation */
    printf("  Test 17: G-Computation...\n");
    CausalEffect* gc = gcomp_estimate(od);
    assert(gc); causal_effect_free(gc);

    /* Test 18: Doubly Robust */
    printf("  Test 18: Doubly Robust...\n");
    CausalEffect* dr = doubly_robust_estimate(od, ps);
    assert(dr); causal_effect_free(dr);

    /* Test 19: Stratification */
    printf("  Test 19: Stratification...\n");
    CausalEffect* st = stratification_estimate(od, ps, 3);
    assert(st); causal_effect_free(st);

    /* Cleanup */
    propensity_free(ps); obs_data_free(od);
    scm_free(scm); dag_free(dag);

    printf("\n=== ALL TESTS PASSED ===\n");
    return 0;
}
