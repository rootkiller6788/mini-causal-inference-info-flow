#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "effect_estimators.h"

int main(void) {
    printf("=== Example 3: Causal Effect Estimation Comparison ===\n\n");
    printf("Simulated RCT with confounding.\n");
    printf("True ATE = 3.0 (treatment effect).\n");
    printf("Confounder: X[0] affects both treatment and outcome.\n\n");

    /* Generate observational data: n=200 units, d=2 covariates
       Treatment T confounded by X[0]: P(T=1) depends on X[0]
       Outcome Y = 2 + 3*T + 0.5*X[0] + noise */
    int n = 200, d = 2;
    ObservationalData* od = obs_data_create(n, d);
    if (!od) { printf("Failed to create data\n"); return 1; }

    for (int i = 0; i < n; i++) {
        /* Covariate 0: continuous confounder U[0,1] */
        od->X[i * d] = (double)rand() / (double)RAND_MAX;
        /* Covariate 1: binary, independent of treatment */
        od->X[i * d + 1] = ((double)rand() / RAND_MAX > 0.5) ? 1.0 : 0.0;

        /* Treatment: confounded by X[0] */
        od->T[i] = (od->X[i * d] > 0.5) ? 1 : 0;

        /* Outcome: true ATE=3, confounded by X[0] */
        double noise = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        od->Y[i] = 2.0 + 3.0 * od->T[i] + 0.5 * od->X[i * d] + noise;
    }

    /* Estimate propensity scores via logistic regression */
    double* ps = propensity_score_estimate(od);
    if (!ps) { printf("PS estimation failed\n"); obs_data_free(od); return 1; }

    /* Compare all estimators */
    CausalEffect* nv  = naive_estimate(od);
    printf("Naive (difference-in-means):\n");
    printf("  ATE = %.4f (SE = %.4f) - BIASED, ignores confounding\n\n",
           nv->ate, nv->se_ate);

    CausalEffect* ipw = ipw_estimate(od, ps);
    printf("IPW (stabilized Horvitz-Thompson):\n");
    printf("  ATE = %.4f (SE = %.4f) - consistent if PS model correct\n\n",
           ipw->ate, ipw->se_ate);

    CausalEffect* gc  = gcomp_estimate(od);
    printf("G-Computation (OLS outcome model):\n");
    printf("  ATE = %.4f (SE = %.4f) - consistent if outcome model correct\n\n",
           gc->ate, gc->se_ate);

    CausalEffect* dr  = doubly_robust_estimate(od, ps);
    printf("Doubly Robust (IPW + outcome model):\n");
    printf("  ATE = %.4f (SE = %.4f) - consistent if EITHER model correct\n\n",
           dr->ate, dr->se_ate);

    CausalEffect* st  = stratification_estimate(od, ps, 5);
    printf("Stratification (K=5 quantiles):\n");
    printf("  ATE = %.4f (SE = %.4f) - ~80%% bias removal\n\n",
           st->ate, st->se_ate);

    CausalEffect* mt  = matching_estimate(od, ps);
    printf("PS Matching (1:1 NN, caliper=0.25):\n");
    printf("  ATE = %.4f (SE = %.4f)\n\n",
           mt->ate, mt->se_ate);

    /* Bootstrap SE for DR estimator */
    double dr_point;
    double dr_boot_se = bootstrap_se(od, ps, doubly_robust_estimate,
                                      100, &dr_point);
    printf("DR Bootstrap SE: %.4f (100 reps)\n\n", dr_boot_se);

    /* Balance check */
    double smd[2];
    balance_check_smd(od, ps, smd, 2);
    printf("Covariate Balance (target |SMD| < 0.1):\n");
    printf("  SMD[X0] = %.4f %s\n", smd[0],
           fabs(smd[0]) < 0.1 ? "(balanced)" : "(imbalanced!)");
    printf("  SMD[X1] = %.4f %s\n\n", smd[1],
           fabs(smd[1]) < 0.1 ? "(balanced)" : "(imbalanced!)");

    /* Summary: compare bias */
    printf("=== Bias Comparison (true ATE = 3.0) ===\n");
    printf("  Naive bias:        %+.4f\n", nv->ate - 3.0);
    printf("  IPW bias:          %+.4f\n", ipw->ate - 3.0);
    printf("  G-Computation bias:%+.4f\n", gc->ate - 3.0);
    printf("  Doubly Robust bias:%+.4f\n", dr->ate - 3.0);
    printf("  Stratification bias:%+.4f\n", st->ate - 3.0);
    printf("  Matching bias:     %+.4f\n", mt->ate - 3.0);

    /* Cleanup */
    causal_effect_free(nv); causal_effect_free(ipw);
    causal_effect_free(gc); causal_effect_free(dr);
    causal_effect_free(st); causal_effect_free(mt);
    propensity_free(ps); obs_data_free(od);

    printf("\n=== Example 3 Complete ===\n");
    return 0;
}
