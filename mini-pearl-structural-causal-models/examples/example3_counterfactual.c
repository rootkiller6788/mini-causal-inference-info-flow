#include "scm_core.h"
#include "scm_counterfactual.h"
#include <stdio.h>
int main(void) {
    printf("========================================\n");
    printf("  Pearl SCM - Demo 3: Counterfactuals\n");
    printf("========================================\n\n");
    SCM_Model* m=scm_create();
    int z=scm_add_variable(m,"Z",SCM_EXOGENOUS);
    int x=scm_add_variable(m,"X",SCM_ENDOGENOUS);
    int y=scm_add_variable(m,"Y",SCM_ENDOGENOUS);
    scm_add_edge(m,z,x,0.5); scm_add_edge(m,z,y,0.3); scm_add_edge(m,x,y,0.7);
    scm_set_equation(m,x,scm_linear_eq); scm_set_equation(m,y,scm_linear_eq);
    m->vars[z].noise=0.5; scm_sample(m);
    scm_print(m);
    int ev[]={x}; double evv[]={2.0};
    int iv[]={x}; double ivv[]={1.0};
    SCM_Counterfactual* cf=scm_counterfactual_compute(m,y,ev,evv,1,iv,ivv,1);
    scm_counterfactual_print(cf);
    SCM_ProbNecessitySuff* pns=scm_prob_necessity_sufficiency(m,x,y,1.0,1.0);
    scm_pns_print(pns);
    double nde,nie; scm_mediation_analysis(m,x,z,y,&nde,&nie);
    printf("NDE=%.4f NIE=%.4f\n", nde, nie);
    /* Effect of Treatment on the Treated (ETT) */
    double ett = scm_effect_of_treatment_on_treated(m, x, y);
    printf("ETT: %.4f\n", ett);
    /* 3-step counterfactual: abduction, action, prediction */
    printf("\n--- 3-Step Counterfactual Process ---\n");
    double u_post[SCM_MAX_VARS];
    int ev2[] = {y}; double evv2[] = {1.5};
    scm_abduction_step(m, ev2, evv2, 1, u_post);
    printf("Step 1 (Abduction): U posterior = {");
    for (int i = 0; i < m->n_vars; i++) printf("%s%.3f", i?",":"", u_post[i]);
    printf("}\n");
    SCM_Model mutilated;
    int iv2[] = {x}; double ivv2[] = {0.0};
    scm_action_step(m, iv2, ivv2, 1, &mutilated);
    printf("Step 2 (Action): mutilated graph, do(X=0)\n");
    double pred;
    scm_prediction_step(&mutilated, y, &pred);
    printf("Step 3 (Prediction): Y(U|do(X=0)) = %.4f\n", pred);
    /* Monotonicity check */
    printf("Monotonic X->Y: %s\n", scm_monotonicity_check(m, x, y)?"YES":"NO");
    scm_counterfactual_free(cf); scm_pns_free(pns); scm_free(m);
    printf("\n========================================\n");
    return 0;
}
/*
 * Pearl SCM Demo 3: Counterfactuals & Mediation
 * Demonstrates: counterfactual computation, 3-step process
 *   (abduction/action/prediction), PNS, mediation (NDE/NIE),
 *   effect on the treated (ETT), monotonicity.
 * References:
 *   Pearl (2000) Causality, Ch. 7 (Counterfactuals).
 *   Pearl (2001) Direct and indirect effects, UAI-01.
 *   VanderWeele (2015) Explanation in Causal Inference.
 */
