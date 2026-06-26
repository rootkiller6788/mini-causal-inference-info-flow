#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "dci_core.h"
#include "dci_graph.h"
#include "dci_do.h"
#include "dci_backdoor.h"
#include "dci_effect.h"
#include "dci_counterfactual.h"

#define EPS 1e-6

static double simple_eq(const double* p, const double* u) {
    return (p ? p[0] : 0.0) + (u ? *u : 0.0);
}

static void t1(void) { /* SCM create */
    DCI_SCM* s = dci_scm_create("test");
    assert(s != NULL); dci_scm_free(s);
}

static void t2(void) { /* Add variables */
    DCI_SCM* s = dci_scm_create("t");
    int p[] = {};
    int id0 = dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, simple_eq, false);
    int id1 = dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, p, 0, simple_eq, false);
    assert(id0 == 0); assert(id1 == 1); assert(s->n_vars == 2);
    dci_scm_free(s);
}

static void t3(void) { /* DAG check */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, NULL, true);
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, NULL, 0, simple_eq, false);
    dci_scm_add_edge(s, 0, 1);
    assert(dci_scm_is_dag(s));
    dci_scm_free(s);
}

static void t4(void) { /* Graph operations */
    DCI_Graph* g = dci_graph_create(3);
    dci_graph_add_edge(g, 0, 1, 1.0);
    dci_graph_add_edge(g, 1, 2, 1.0);
    dci_graph_add_edge(g, 0, 2, 1.0);
    assert(g->adjacency[0][1] == 1.0);
    assert(dci_graph_is_dag(g));
    dci_graph_free(g);
}

static void t5(void) { /* d-separation — chain X→M→Y */
    DCI_Graph* g = dci_graph_create(3);
    dci_graph_add_edge(g, 0, 1, 1.0);
    dci_graph_add_edge(g, 1, 2, 1.0);
    int Z[] = {1};
    assert(dci_is_d_separated(g, 0, 2, Z, 1)); /* conditioned on M */
    assert(!dci_is_d_separated(g, 0, 2, NULL, 0)); /* unconditioned */
    dci_graph_free(g);
}

static void t6(void) { /* d-separation — fork X←U→Y */
    DCI_Graph* g = dci_graph_create(3);
    dci_graph_add_edge(g, 1, 0, 1.0);
    dci_graph_add_edge(g, 1, 2, 1.0);
    int Z[] = {1};
    assert(dci_is_d_separated(g, 0, 2, Z, 1));
    dci_graph_free(g);
}

static void t7(void) { /* d-separation — collider X→M←Y */
    DCI_Graph* g = dci_graph_create(3);
    dci_graph_add_edge(g, 0, 1, 1.0);
    dci_graph_add_edge(g, 2, 1, 1.0);
    assert(dci_is_d_separated(g, 0, 2, NULL, 0));
    int Z[] = {1};
    assert(!dci_is_d_separated(g, 0, 2, Z, 1)); /* conditioning opens */
    dci_graph_free(g);
}

static void t8(void) { /* Back-door criterion */
    DCI_Graph* g = dci_graph_create(3);
    dci_graph_add_edge(g, 0, 1, 1.0); /* X→Y */
    dci_graph_add_edge(g, 2, 0, 1.0); /* Z→X (confounder) */
    dci_graph_add_edge(g, 2, 1, 1.0); /* Z→Y */
    DCI_BackdoorSet bd = dci_backdoor_find(g, 0, 1);
    assert(bd.is_valid);
    dci_graph_free(g);
}

static void t9(void) { /* Intervention and truncated SCM */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, simple_eq, true);
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, NULL, 0, simple_eq, false);
    int vars[] = {0}; double vals[] = {1.0};
    DCI_Intervention iv = dci_intervention_create(vars, vals, 1);
    assert(iv.n_vars == 1);
    DCI_TruncatedSCM* ts = dci_truncated_scm_create(s, &iv);
    assert(ts != NULL);
    dci_truncated_scm_free(ts);
    dci_scm_free(s);
}

static void t10(void) { /* ATE computation */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, NULL, true);
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, NULL, 0, simple_eq, false);
    DCI_CausalEffect ce = dci_compute_ate(s, 0, 1, 100);
    assert(ce.is_identifiable);
    dci_scm_free(s);
}

static void t11(void) { /* Find directed paths */
    DCI_Graph* g = dci_graph_create(4);
    dci_graph_add_edge(g, 0, 1, 1.0);
    dci_graph_add_edge(g, 1, 2, 1.0);
    dci_graph_add_edge(g, 1, 3, 1.0);
    dci_graph_add_edge(g, 2, 3, 1.0);
    DCI_Path paths[16];
    int n = dci_find_directed_paths(g, 0, 3, paths, 16);
    assert(n > 0);
    dci_graph_free(g);
}

static void t12(void) { /* Back-door adjustment */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "Z", DCI_CONTINUOUS, NULL, 0, NULL, true);
    int px[] = {0};
    dci_scm_add_variable(s, "X", DCI_BINARY, px, 1, simple_eq, false);
    int py[] = {1, 0};
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, py, 2, simple_eq, false);
    DCI_Graph g;
    dci_graph_from_scm(s, &g);
    DCI_BackdoorSet bd = dci_backdoor_find(&g, 1, 2);
    double eff = dci_backdoor_adjust(s, &g, 1, 2, 1.0, 0.0, &bd, 100);
    assert(eff >= 0);
    dci_scm_free(s);
}

static void t13(void) { /* Counterfactual */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, NULL, true);
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, NULL, 0, simple_eq, false);
    DCI_Counterfactual cf = {.antecedent_var=0, .antecedent_value=1.0,
        .consequent_var=1, .consequent_value=0.0};
    DCI_CounterfactualResult r = dci_counterfactual_probability(s, &cf, NULL, NULL, 0, 50);
    assert(r.n_samples == 50);
    dci_scm_free(s);
}

static void t14(void) { /* do-calculus rule application */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, NULL, true);
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, NULL, 0, simple_eq, false);
    DCI_Graph g;
    dci_graph_from_scm(s, &g);
    int x[]={0}; int z[]={1};
    DCI_Derivation der = dci_apply_rule1(s, &g, 1, x, 1, z, 1, NULL, 0);
    assert(der.n_steps >= 0);
    dci_scm_free(s);
}

static void t15(void) { /* SCM factories */
    DCI_SCM* s1 = dci_scm_simpsons_paradox();
    assert(s1 != NULL && s1->n_vars >= 3);
    dci_scm_free(s1);
    DCI_SCM* s2 = dci_scm_iv(100);
    assert(s2 != NULL && s2->n_vars >= 4);
    dci_scm_free(s2);
}

static void t16(void) { /* Graph surgery */
    DCI_Graph* g = dci_graph_create(3);
    dci_graph_add_edge(g, 0, 1, 1.0);
    dci_graph_add_edge(g, 2, 1, 1.0);
    DCI_Graph g_mut;
    dci_graph_mutilate(g, (int[]){0}, 1, &g_mut);
    assert(g_mut.adjacency[0][1] > 0.5); /* outgoing preserved */
    dci_graph_free(g);
}

static void t17(void) { /* Ancestor check */
    DCI_Graph* g = dci_graph_create(3);
    dci_graph_add_edge(g, 0, 1, 1.0);
    dci_graph_add_edge(g, 1, 2, 1.0);
    assert(dci_is_ancestor(g, 0, 2));
    assert(!dci_is_ancestor(g, 2, 0));
    dci_graph_free(g);
}

static void t18(void) { /* Matching estimator */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, NULL, true);
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, NULL, 0, simple_eq, false);
    double m = dci_matching_estimator(s, 0, 1, NULL, 0, 100);
    assert(isfinite(m));
    dci_scm_free(s);
}

static void t19(void) { /* Bootstrap CI */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, NULL, true);
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, NULL, 0, simple_eq, false);
    DCI_ConfidenceInterval ci = dci_ate_bootstrap(s, 0, 1, 100, 20, 0.05);
    assert(ci.estimate > -1e6);
    dci_scm_free(s);
}

static void t20(void) { /* V-structure detection */
    DCI_Graph* g = dci_graph_create(3);
    dci_graph_add_edge(g, 0, 1, 1.0); /* X→M */
    dci_graph_add_edge(g, 2, 1, 1.0); /* Y→M (collider at M) */
    int triples[8][3];
    int n = dci_find_v_structures(g, triples, 8);
    assert(n == 1); /* One v-structure: X→M←Y */
    dci_graph_free(g);
}

static void t21(void) { /* dci_scm_var_count */
    DCI_SCM* s = dci_scm_create("t");
    dci_scm_add_variable(s, "X", DCI_BINARY, NULL, 0, NULL, true);
    dci_scm_add_variable(s, "Y", DCI_CONTINUOUS, NULL, 0, NULL, true);
    assert(dci_scm_var_count(s) == 2);
    dci_scm_free(s);
}

int main(void) {
    t1(); t2(); t3(); t4(); t5(); t6(); t7(); t8(); t9(); t10();
    t11(); t12(); t13(); t14(); t15(); t16(); t17(); t18(); t19(); t20();
    t21();
    printf("All 21 tests passed\n");
    return 0;
}
