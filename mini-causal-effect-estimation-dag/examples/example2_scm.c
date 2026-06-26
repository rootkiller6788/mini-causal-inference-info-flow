#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scm_model.h"

int main(void) {
    printf("=== Example 2: Structural Causal Model Simulation ===\n");
    printf("SCM: X0 -> X1 -> X2 (linear with Gaussian noise)\n");
    printf("Question: what is the ACE of X1 on X2?\n\n");

    /* Build SCM with 3 variables:
       X0: exogenous (root)
       X1 = X0 + noise_1 (depends on X0)
       X2 = X1 + noise_2 (depends on X1) */
    SCModel* scm = scm_create(3);
    if (!scm) { printf("Failed to create SCM\n"); return 1; }

    dag_add_edge(scm->dag, 0, 1);  /* X0 -> X1 */
    dag_add_edge(scm->dag, 1, 2);  /* X1 -> X2 */

    scm_set_equation(scm, 0, scm_linear_eq, 0);  /* X0 = noise (root) */
    scm_set_equation(scm, 1, scm_linear_eq, 1);  /* X1 = X0 + noise */
    scm_set_equation(scm, 2, scm_linear_eq, 2);  /* X2 = X1 + noise */

    /* Set noise terms: X0=2.0, X1_noise=0.0, X2_noise=0.0 */
    scm->noise[0] = 2.0;
    scm->noise[1] = 0.0;
    scm->noise[2] = 0.0;

    /* Simulate the SCM */
    if (scm_simulate(scm) == 0) {
        printf("Simulation results:\n");
        printf("  X0 = %.2f (root, exogenous)\n", scm->values[0]);
        printf("  X1 = %.2f (= X0 + noise_1)\n", scm->values[1]);
        printf("  X2 = %.2f (= X1 + noise_2)\n\n", scm->values[2]);
    }

    /* Estimate ACE: E[X2 | do(X1=1)] - E[X2 | do(X1=0)]
       With 500 Monte Carlo samples and Gaussian noise. */
    double ace = scm_average_causal_effect(scm, 1, 1.0, 0.0, 2, 500);
    printf("ACE of X1 on X2: %.3f\n", ace);
    printf("  Interpretation: Setting X1 from 0 to 1 changes X2 by %.3f units\n\n",
           ace);

    /* Verify with do-intervention directly */
    SCModel* scm_do1 = scm_do_intervention(scm, 1, 10.0);
    scm_simulate(scm_do1);
    printf("After do(X1=10.0): X2 = %.2f\n", scm_do1->values[2]);

    SCModel* scm_do0 = scm_do_intervention(scm, 1, 0.0);
    scm_simulate(scm_do0);
    printf("After do(X1=0.0):  X2 = %.2f\n", scm_do0->values[2]);

    /* Demonstrate d-separation in DAG
       In this chain X0->X1->X2:
       - X0 and X2 are d-separated by X1 (conditioning on X1 blocks the path)
       - Without conditioning, X0 and X2 are d-connected */
    int z[1] = { 1 };
    DSepResult* ds1 = dag_d_separation(scm->dag, 0, 2, z, 1);
    printf("\nd-separation test: X0 and X2 given X1...\n");
    printf("  d-separated: %s\n", ds1->is_separated ? "YES" : "NO");
    if (ds1->blocking_path) {
        free(ds1->blocking_path->nodes);
        free(ds1->blocking_path);
    }
    free(ds1);

    DSepResult* ds0 = dag_d_separation(scm->dag, 0, 2, NULL, 0);
    printf("d-separation test: X0 and X2 unconditionally...\n");
    printf("  d-separated: %s\n", ds0->is_separated ? "YES" : "NO");
    if (ds0->blocking_path) {
        free(ds0->blocking_path->nodes);
        free(ds0->blocking_path);
    }
    free(ds0);

    scm_free(scm_do1); scm_free(scm_do0); scm_free(scm);
    printf("\n=== Example 2 Complete ===\n");
    return 0;
}
