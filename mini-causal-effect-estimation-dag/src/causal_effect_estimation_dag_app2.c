/*
 * causal_effect_estimation_dag_app2.c ¡ª L7 Application: Economics / Policy
 *
 * Estimates the causal effect of a job training program (e.g., Detroit
 * workforce development 2020) on subsequent earnings, adjusting for
 * confounding by education level, prior work experience, and motivation.
 *
 * Causal DAG:
 *   Education ©¤©¤¡ú Training ©¤©¤¡ú Earnings
 *       ©¦            ¡ü           ©¦
 *   Experience ©¤©¤©¤©¤¡ú©¦©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤¡ú©¦
 *       ©¦            ©¦           ©¦
 *   Motivation ©¤©¤©¤©¤¡ú©¦©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤¡ú©¦
 *
 * With front-door: Training ¡ú Skills ¡ú Earnings
 *
 * References:
 *   LaLonde, "Evaluating the Econometric Evaluations of Training
 *     Programs", AER, 1986
 *   Dehejia & Wahba, "Causal Effects in Nonexperimental Studies",
 *     JASA, 1999
 *   Pearl, "Causality", 2009, Ch.3 (Front-Door Criterion)
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
    printf("=== L7 Application 2: Job Training Effect (Economics) ===\n\n");

    /* 1. DAG: 5 nodes
     *   X0 = Education, X1 = Training, X2 = Earnings,
     *   X3 = Experience, X4 = Motivation */
    DAG* dag = dag_create(5);
    dag_add_edge(dag, 0, 1);  /* Education -> Training */
    dag_add_edge(dag, 3, 1);  /* Experience -> Training */
    dag_add_edge(dag, 4, 1);  /* Motivation -> Training */
    dag_add_edge(dag, 1, 2);  /* Training -> Earnings */
    dag_add_edge(dag, 0, 2);  /* Education -> Earnings */
    dag_add_edge(dag, 3, 2);  /* Experience -> Earnings */
    dag_add_edge(dag, 4, 2);  /* Motivation -> Earnings */

    printf("Causal DAG: Education+Experience+Motivation -> Training -> Earnings\n");
    printf("            Education+Experience+Motivation -> Earnings (direct)\n\n");

    /* 2. Causal identification */
    IdentificationResult* ir = causal_identify(dag, 1, 2);
    printf("Effect of Training on Earnings:\n");
    printf("  Identifiable: %s\n", ir->identifiable ? "YES" : "NO");
    printf("  Method: %d (1=Backdoor, 2=Frontdoor, 3=IV, 4=Do-Calc)\n", ir->method);
    printf("  Adjustment set size: %d\n", ir->n_adjust);
    if (ir->explanation) printf("  %s\n\n", ir->explanation);
    ident_result_free(ir);

    /* 3. Find minimal adjustment sets */
    AdjustmentSet min_sets[ADJ_MAX_SETS];
    int n_min = adjustment_find_minimal_sets(dag, 1, 2, min_sets, ADJ_MAX_SETS);
    printf("Minimal adjustment sets found: %d\n", n_min);
    for (int s = 0; s < n_min; s++) {
        printf("  Set %d: {", s + 1);
        for (int v = 0; v < min_sets[s].n_vars; v++)
            printf(" X%d%s", min_sets[s].variables[v],
                   v < min_sets[s].n_vars - 1 ? "," : " ");
        printf("} (size=%d, efficiency=%.3f)\n",
               min_sets[s].n_vars, min_sets[s].efficiency);
        free(min_sets[s].variables);
    }
    printf("\n");

    /* 4. Generate observational data (Detroit workforce program simulation) */
    int n = 500, d = 3;  /* Education, Experience, Motivation */
    ObservationalData* obs = obs_data_create(n, d);

    for (int i = 0; i < n; i++) {
        double edu   = 10.0 + 8.0 * ((double)rand() / RAND_MAX);  /* 10-18 years */
        double exper =  0.0 + 15.0 * ((double)rand() / RAND_MAX);  /* 0-15 years */
        double motiv = -2.0 + 4.0 * ((double)rand() / RAND_MAX);  /* standardized */

        obs->X[i * d] = edu;
        obs->X[i * d + 1] = exper;
        obs->X[i * d + 2] = motiv;

        /* Training assignment (confounded) */
        double logit_t = -3.0 + 0.15*edu + 0.1*exper + 0.5*motiv;
        double p_train = 1.0 / (1.0 + exp(-logit_t));
        obs->T[i] = (((double)rand() / RAND_MAX) < p_train) ? 1 : 0;

        /* Earnings: $20k-60k base + $8k training effect + covariates + noise */
        double noise = ((double)rand()/RAND_MAX - 0.5) * 10.0;
        obs->Y[i] = 30.0 + 8.0*obs->T[i] + 1.5*edu + 0.8*exper
                    + 2.0*motiv + noise;
    }

    /* 5. Estimation */
    double* ps = propensity_score_estimate(obs);

    CausalEffect* naive = naive_estimate(obs);
    printf("NAIVE ATE:          $%.1fk (SE=$%.2fk)\n",
           naive->ate, naive->se_ate);

    CausalEffect* ipw = ipw_estimate(obs, ps);
    printf("IPW ATE:            $%.1fk (SE=$%.2fk)\n",
           ipw->ate, ipw->se_ate);

    CausalEffect* gcomp = gcomp_estimate(obs);
    printf("G-Comp ATE:         $%.1fk (SE=$%.2fk)\n",
           gcomp->ate, gcomp->se_ate);

    CausalEffect* dr = doubly_robust_estimate(obs, ps);
    printf("Doubly Robust ATE:  $%.1fk (SE=$%.2fk)\n",
           dr->ate, dr->se_ate);

    CausalEffect* strat = stratification_estimate(obs, ps, 10);
    printf("Stratification ATE: $%.1fk (SE=$%.2fk) K=10 (~90%% bias removed)\n",
           strat->ate, strat->se_ate);

    /* Bootstrap SE for DR estimator */
    double dr_point;
    double dr_boot_se = bootstrap_se(obs, ps, doubly_robust_estimate, 200, &dr_point);
    printf("DR Bootstrap SE:    $%.2fk (200 reps)\n", dr_boot_se);

    printf("\nTrue ATE = $8.0k (training increases earnings by $8,000)\n");
    printf("Naive bias = $%.2fk | IPW bias = $%.2fk | DR bias = $%.2fk\n",
           naive->ate - 8.0, ipw->ate - 8.0, dr->ate - 8.0);

    /* Covariate balance after PS adjustment */
    double smd[3];
    balance_check_smd(obs, ps, smd, 3);
    printf("\nBalance check (SMD, target <0.1):\n");
    printf("  Education:  %.4f\n", smd[0]);
    printf("  Experience: %.4f\n", smd[1]);
    printf("  Motivation: %.4f\n", smd[2]);

    /* Cleanup */
    causal_effect_free(naive); causal_effect_free(ipw);
    causal_effect_free(gcomp); causal_effect_free(dr);
    causal_effect_free(strat);
    propensity_free(ps); obs_data_free(obs); dag_free(dag);

    printf("\n=== Application 2 Complete ===\n");
    return 0;
}
