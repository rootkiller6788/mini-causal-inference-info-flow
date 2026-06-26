#include "scm_core.h"
#include "scm_intervention.h"
#include "scm_identification.h"
#include <stdio.h>
int main(void) {
    printf("========================================\n");
    printf("  Pearl SCM - Demo 2: Back-Door Adj.\n");
    printf("========================================\n\n");
    SCM_Model* m=scm_create();
    int z=scm_add_variable(m,"Z",SCM_EXOGENOUS);
    int x=scm_add_variable(m,"X",SCM_ENDOGENOUS);
    int y=scm_add_variable(m,"Y",SCM_ENDOGENOUS);
    scm_add_edge(m,z,x,0.6); scm_add_edge(m,z,y,0.4); scm_add_edge(m,x,y,0.5);
    scm_set_equation(m,x,scm_linear_eq); scm_set_equation(m,y,scm_linear_eq);
    m->vars[z].noise=0.5; scm_sample(m);
    SCM_CausalEffect* ce=scm_causal_effect_compute(m,x,y);
    scm_causal_effect_print(ce);
    printf("Linear ACE: %.4f\n", scm_ace_linear(m,x,y));
    SCM_BackDoorResult* br=scm_back_door(m,x,y);
    scm_back_door_print(br);
    printf("Identifiable: %s\n", scm_is_identifiable(m,x,y)?"YES":"NO");
    /* Find all valid adjustment sets */
    SCM_VarSet sets[8]; int ns;
    scm_adjustment_sets(m, x, y, sets, &ns, 8);
    printf("Adjustment sets found: %d\n", ns);
    /* Check minimal adjustment */
    if (ns > 0) {
        printf("Is minimal adjustment: %s\n",
            scm_is_minimal_adjustment(m, x, y, &sets[0])?"YES":"NO");
    }
    /* do-Calculus rules verification */
    SCM_VarSet empty; scm_varset_init(&empty);
    bool rule1, rule2, rule3;
    scm_do_calculus_rule1(m, x, y, &empty, &empty, &rule1);
    scm_do_calculus_rule2(m, x, y, &empty, &rule2);
    scm_do_calculus_rule3(m, x, y, &empty, &empty, &rule3);
    printf("do-Calculus Rule1: %s  Rule2: %s  Rule3: %s\n",
        rule1?"YES":"NO", rule2?"YES":"NO", rule3?"YES":"NO");
    /* Find instrumental variable */
    int instr;
    scm_find_instrument(m, x, y, &instr);
    printf("Instrument found: %s (var=%d)\n", instr>=0?"YES":"NO", instr);
    /* Test for unobserved confounder */
    bool has_uc;
    scm_unobserved_confounder_test(m, x, y, &has_uc);
    printf("Unobserved confounder: %s\n", has_uc?"DETECTED":"none");
    /* IPW effect estimation */
    printf("Linear ACE (from path coeffs): %.4f\n", scm_ace_linear(m, x, y));
    scm_causal_effect_free(ce); scm_back_door_free(br); scm_free(m);
    printf("\n========================================\n");
    return 0;
}
/*
 * Pearl SCM Demo 2: Back-Door Adjustment & Identification
 * Demonstrates: causal effect computation, back-door adjustment,
 *   do-calculus rules, instrumental variables, unobserved confounding.
 * References:
 *   Pearl (2009) Causality, Ch. 3 (do-calculus).
 *   Hernan & Robins (2020) Causal Inference: What If.
 *   Angrist & Pischke (2009) Mostly Harmless Econometrics.
 */
