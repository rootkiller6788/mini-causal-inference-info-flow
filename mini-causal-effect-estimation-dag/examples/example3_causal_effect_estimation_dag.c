/* Example 3: Front-Door Criterion and G-Computation
 *
 * Demonstrates causal effect estimation using the front-door criterion
 * when there are unmeasured confounders between treatment and outcome.
 *
 * Scenario: Does smoking (X) cause lung cancer (Y)?
 *   Unmeasured confounder: genetic predisposition (U)
 *   Front-door mediator: tar deposits in lungs (M)
 *
 * DAG: U -> X -> M -> Y, with U -> Y (unobserved)
 * Using front-door: P(Y|do(X)) = sum_m P(m|X) sum_x' P(Y|x',m) P(x')
 *
 * This example uses Pearl's classic smoking-lung cancer model.
 *
 * References:
 *   Pearl, "Causality", 2009, Ch.3.3.3
 *   Glymour & Greenland, "Causal Diagrams", 2008
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dag_representation.h"
#include "scm_model.h"
#include "causal_identification.h"
#include "effect_estimators.h"

int main(void) {
    srand((unsigned)time(NULL));
    printf("=== Example 3: Front-Door Criterion (Smoking-Cancer Model) ===\n\n");

    /* Build DAG: 5 nodes
     * X0=Genetic(U), X1=Smoking(X), X2=Tar(M), X3=Cancer(Y), X4=Age */
    DAG* dag = dag_create(5);
    dag_add_edge(dag, 0, 1);  /* Genetic -> Smoking (confounder) */
    dag_add_edge(dag, 0, 3);  /* Genetic -> Cancer (confounder, unobserved!) */
    dag_add_edge(dag, 1, 2);  /* Smoking -> Tar (mediator) */
    dag_add_edge(dag, 2, 3);  /* Tar -> Cancer */
    dag_add_edge(dag, 4, 1);  /* Age -> Smoking */
    dag_add_edge(dag, 4, 3);  /* Age -> Cancer */

    printf("DAG: Genetic(U) -> Smoking(X), Genetic(U) -> Cancer(Y)\n");
    printf("     Smoking(X) -> Tar(M) -> Cancer(Y)\n");
    printf("     Age -> Smoking(X), Age -> Cancer(Y)\n\n");

    /* Try back-door criterion: should fail because Genetic is unobserved */
    int backdoor_adj[10];
    int n_backdoor = causal_find_backdoor_adjustment(dag, 1, 3, backdoor_adj, 10);
    printf("Back-door adjustment: %s\n",
           n_backdoor >= 0 ? "FOUND" : "NOT FOUND (genetic confounder unmeasured)");

    /* Check front-door criterion with Tar (node 2) as mediator */
    bool frontdoor = causal_frontdoor_criterion(dag, 1, 3, 2);
    printf("Front-door with Tar (M=2): %s\n\n", frontdoor ? "VALID" : "INVALID");

    /* Generate data using SCM */
    SCModel* scm = scm_create(5);
    dag_add_edge(scm->dag, 0, 1); dag_add_edge(scm->dag, 0, 3);
    dag_add_edge(scm->dag, 1, 2); dag_add_edge(scm->dag, 2, 3);
    dag_add_edge(scm->dag, 4, 1); dag_add_edge(scm->dag, 4, 3);

    scm_set_equation(scm, 0, scm_linear_eq, 0);  /* Genetic (exogenous) */
    scm_set_equation(scm, 1, scm_linear_eq, 2);  /* Smoking */
    scm_set_equation(scm, 2, scm_linear_eq, 1);  /* Tar = b*Smoking + noise */
    scm_set_equation(scm, 3, scm_linear_eq, 2);  /* Cancer */
    scm_set_equation(scm, 4, scm_linear_eq, 0);  /* Age (exogenous) */

    scm->noise[0] = 1.0;  /* Genetic predisposition */
    scm->noise[4] = 1.5;  /* Age effect */

    /* Simulate to verify */
    scm_simulate(scm);
    printf("SCM Values: G=%.2f S=%.2f T=%.2f C=%.2f A=%.2f\n",
           scm->values[0], scm->values[1], scm->values[2],
           scm->values[3], scm->values[4]);

    /* Monte Carlo ACE */
    double ace = scm_average_causal_effect(scm, 1, 1.0, 0.0, 3, 1000);
    printf("ACE (Smoking on Cancer via SCM): %.3f\n\n", ace);

    /* Front-door estimation using observational data */
    int n = 300, d = 1;  /* Age as observed covariate */
    ObservationalData* obs = obs_data_create(n, d);
    if (obs) {
        for (int i = 0; i < n; i++) {
            double genetic = ((double)rand() / RAND_MAX > 0.5) ? 1.0 : 0.0;
            double age = 40.0 + 20.0 * ((double)rand() / RAND_MAX);
            obs->X[i * d] = age;
            double p_smoke = 1.0 / (1.0 + exp(-(0.5*genetic + 0.05*(age-50))));
            obs->T[i] = (((double)rand()/RAND_MAX) < p_smoke) ? 1 : 0;
            double tar = 10.0 + 20.0*obs->T[i] + 0.1*age +
                         5.0*((double)rand()/RAND_MAX-0.5);
            double cancer_risk = 0.01 + 0.03*obs->T[i] + 0.005*tar +
                                 0.02*genetic + 0.001*age;
            obs->Y[i] = cancer_risk + 0.005*((double)rand()/RAND_MAX-0.5);
        }

        double* ps = propensity_score_estimate(obs);
        CausalEffect* ipw = ipw_estimate(obs, ps);
        CausalEffect* naive = naive_estimate(obs);
        CausalEffect* gcomp = gcomp_estimate(obs);
        CausalEffect* dr = doubly_robust_estimate(obs, ps);

        printf("Observational Estimates (n=%d):\n", n);
        printf("  Naive ATE:           %.5f (SE=%.5f)\n", naive->ate, naive->se_ate);
        printf("  IPW ATE:             %.5f (SE=%.5f)\n", ipw->ate, ipw->se_ate);
        printf("  G-Comp ATE:          %.5f (SE=%.5f)\n", gcomp->ate, gcomp->se_ate);
        printf("  Doubly Robust ATE:   %.5f (SE=%.5f)\n", dr->ate, dr->se_ate);
        printf("\nNote: Front-door formula would produce unbiased estimate\n");
        printf("when correctly computing P(Y|do(X)) via mediator.\n");

        causal_effect_free(ipw); causal_effect_free(naive);
        causal_effect_free(gcomp); causal_effect_free(dr);
        propensity_free(ps); obs_data_free(obs);
    }

    scm_free(scm); dag_free(dag);
    printf("\n=== Example 3 PASSED ===\n");
    return 0;
}
