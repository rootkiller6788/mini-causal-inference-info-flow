/* ============================================================
 * Example 1: Structural Causal Models & Counterfactuals
 *
 * Demonstrates:
 *   - SCM creation and evaluation
 *   - Hard/soft interventions
 *   - Pearl's 3-step counterfactual procedure
 *   - Counterfactual comparisons
 *   - DAG analysis (parents, children, DAG check)
 *
 * Refs: Pearl (2009) Causality, Ch.7.
 * ============================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cfr_core.h"
#include "cfr_counterfactual.h"

int main(void) {
    printf("=== Example 1: SCM & Counterfactual Reasoning ===\n\n");

    /* 1. Create a simple treatment SCM */
    CFRSCM* scm = cfr_create_simple_treatment_scm();
    printf("1. Created SCM: %s (V=%d, U=%d)\n", scm->name,
           scm->n_endogenous, scm->n_exogenous);

    /* 2. Set exogenous and compute factual */
    double u[] = {0.8, 1.5};
    cfr_scm_set_exogenous(scm, u);
    cfr_scm_compute(scm);
    printf("2. Factual: %s=%.1f, %s=%.2f\n",
           scm->vars[0].name, scm->vars[0].value,
           scm->vars[1].name, scm->vars[1].value);

    /* 3. Hard intervention: do(T=0) */
    CFRIntervention* intv_t0 = cfr_int_create_hard(0, 0.0);
    CFRSCM* scm_t0 = cfr_scm_apply_intervention(scm, intv_t0);
    cfr_scm_compute(scm_t0);
    printf("3. After do(T=0): %s=%.2f\n", scm_t0->vars[1].name,
           scm_t0->vars[1].value);

    /* 4. Counterfactual: what if T=0 when we observed T=1? */
    CFREvidence* ev = cfr_ev_create(2);
    cfr_ev_add_observation(ev, 0, 1.0, 0); /* T=1 observed */
    cfr_ev_add_observation(ev, 1, 5.0, 1); /* Y=5.0 observed */
    CFRCounterfactual* cf = cfr_cf_create(scm);
    cfr_cf_compute_full(cf, ev, intv_t0, 1);
    printf("4. Counterfactual Y_T=0 given observed T=1,Y=5: %.2f\n", cf->mean);

    /* 5. Compare two interventions */
    CFRIntervention* intv_t1 = cfr_int_create_hard(0, 1.0);
    double diff;
    cfr_cf_compare_interventions(scm, 1, intv_t1, intv_t0, ev, &diff);
    printf("5. ATE from CF comparison: %.2f\n", diff);

    /* 6. DAG analysis */
    bool is_dag = cfr_scm_is_dag(scm);
    printf("6. Is DAG: %s\n", is_dag ? "yes" : "no");

    int parents[10];
    int np = cfr_scm_parents(scm, 1, parents, 10);
    printf("   Parents of Y: %d\n", np);

    int children[10];
    int nc = cfr_scm_children(scm, 0, children, 10);
    printf("   Children of T: %d (Y[%d])\n", nc, children[0]);

    /* 7. Collider detection */
    CFRSCM* fdoor = cfr_create_frontdoor_scm();
    int children_x[10];
    cfr_scm_children(fdoor, 0, children_x, 10);
    bool is_collider = cfr_scm_is_collider(fdoor, 1, 0, 3);
    printf("7. M is collider (X->M<-U)? %s\n", is_collider ? "yes" : "no");

    /* 8. Information flow */
    double flow = cfr_scm_information_flow(scm, 0, 1);
    printf("8. Information flow T->Y: %.2f\n", flow);

    /* Cleanup */
    cfr_scm_free(scm); cfr_scm_free(scm_t0); cfr_scm_free(fdoor);
    cfr_int_free(intv_t0); cfr_int_free(intv_t1);
    cfr_ev_free(ev); cfr_cf_free(cf);

    printf("\n=== Example 1 PASSED ===\n");
    return 0;
}
