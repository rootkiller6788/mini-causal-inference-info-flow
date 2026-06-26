/*
 * causal_effect_estimation_dag_app1.c ！ L7 Application: Medical
 *
 * Estimates the causal effect of a drug treatment on recovery time
 * using observational clinical data, adjusting for confounders
 * (age, severity) via back-door criterion and multiple estimators.
 *
 * Scenario: Does Drug A reduce recovery time compared to Drug B?
 *   X = treatment (0=Drug B, 1=Drug A)
 *   Y = recovery time (days)
 *   Confounders: age, disease_severity
 *
 * References:
 *   Rosenbaum & Rubin, Biometrika, 1983
 *   Hernan & Robins, "Causal Inference: What If", 2020
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dag_representation.h"
#include "scm_model.h"
#include "causal_identification.h"
#include "adjustment_sets.h"
#include "effect_estimators.h"

int main(void) {
    printf("=== L7 Application 1: Drug Treatment Effect (Medical) ===\n\n");

    /* 1. Build causal DAG: age->treatment, severity->treatment,
     *    treatment->recovery, age->recovery, severity->recovery */
    DAG* dag = dag_create(4);
    dag_add_edge(dag, 0, 1);  /* age(X0) -> treatment(X1) */
    dag_add_edge(dag, 2, 1);  /* severity(X2) -> treatment(X1) */
    dag_add_edge(dag, 1, 3);  /* treatment(X1) -> recovery(X3) */
    dag_add_edge(dag, 0, 3);  /* age(X0) -> recovery(X3) */
    dag_add_edge(dag, 2, 3);  /* severity(X2) -> recovery(X3) */

    printf("Causal DAG: Age + Severity -> Treatment -> Recovery\n");
    printf("Back-door paths: Treatment <- Severity -> Recovery\n");
    printf("                 Treatment <- Age -> Recovery\n\n");

    /* 2. Identify causal effect */
    IdentificationResult* ir = causal_identify(dag, 1, 3);
    printf("Identification: %s\n", ir->identifiable ? "IDENTIFIABLE" : "NOT identifiable");
    printf("Method: %d  Adjust set size: %d\n", ir->method, ir->n_adjust);
    if (ir->explanation) printf("Reason: %s\n\n", ir->explanation);
    ident_result_free(ir);

    /* 3. Simulate observational data (n=200) */
    int n = 200, d = 2;  /* 2 covariates: age, severity */
    ObservationalData* obs = obs_data_create(n, d);
    for (int i = 0; i < n; i++) {
        double age = 40.0 + 20.0 * ((double)rand() / RAND_MAX);   /* 40-60 */
        double sev = 1.0 + 4.0 * ((double)rand() / RAND_MAX);     /* 1-5 */
        obs->X[i * d] = age;
        obs->X[i * d + 1] = sev;
        /* Treatment assignment depends on age and severity */
        double p_treat = 1.0 / (1.0 + exp(-(0.1*(age-50) + 0.5*(sev-3))));
        obs->T[i] = (((double)rand() / RAND_MAX) < p_treat) ? 1 : 0;
        /* Outcome: recovery = 30 - 5*T + 0.1*age + 2*sev + noise */
        double noise = ((double)rand()/RAND_MAX - 0.5) * 4.0;
        obs->Y[i] = 30.0 - 5.0*obs->T[i] + 0.1*age + 2.0*sev + noise;
    }

    /* 4. Estimate propensity scores */
    double* ps = propensity_score_estimate(obs);

    /* 5. Apply multiple estimators */
    CausalEffect* naive = naive_estimate(obs);
    printf("Naive ATE:           %.3f (SE=%.3f) ！ BIASED (ignores confounding)\n",
           naive->ate, naive->se_ate);

    CausalEffect* ipw = ipw_estimate(obs, ps);
    printf("IPW ATE:             %.3f (SE=%.3f) ！ Consistent if PS correct\n",
           ipw->ate, ipw->se_ate);

    CausalEffect* gcomp = gcomp_estimate(obs);
    printf("G-Computation ATE:   %.3f (SE=%.3f) ！ Consistent if outcome correct\n",
           gcomp->ate, gcomp->se_ate);

    CausalEffect* dr = doubly_robust_estimate(obs, ps);
    printf("Doubly Robust ATE:   %.3f (SE=%.3f) ！ Consistent if either correct\n",
           dr->ate, dr->se_ate);

    CausalEffect* strat = stratification_estimate(obs, ps, 5);
    printf("Stratification ATE:  %.3f (SE=%.3f) ！ Approx 80%% bias removed (K=5)\n",
           strat->ate, strat->se_ate);

    CausalEffect* match = matching_estimate(obs, ps);
    printf("Matching ATE:        %.3f (SE=%.3f) (1:1 NN PS matching)\n",
           match->ate, match->se_ate);

    printf("\nTrue ATE = -5.0 (drug reduces recovery by 5 days)\n");
    printf("Naive bias = %.3f | IPW bias = %.3f | DR bias = %.3f\n",
           naive->ate - (-5.0), ipw->ate - (-5.0), dr->ate - (-5.0));

    /* 6. Covariate balance check */
    double smd[2];
    balance_check_smd(obs, ps, smd, 2);
    printf("\nCovariate Balance (SMD): Age=%.3f Severity=%.3f (target <0.1)\n",
           smd[0], smd[1]);

    /* Cleanup */
    causal_effect_free(naive); causal_effect_free(ipw);
    causal_effect_free(gcomp); causal_effect_free(dr);
    causal_effect_free(strat); causal_effect_free(match);
    propensity_free(ps); obs_data_free(obs); dag_free(dag);

    printf("\n=== Application 1 Complete ===\n");
    return 0;
}
