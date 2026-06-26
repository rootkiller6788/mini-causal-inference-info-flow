/* ============================================================
 * Example 2: Causal Effect Estimation Methods
 *
 * Demonstrates:
 *   - Oracle ATE (known potential outcomes)
 *   - Difference-in-means
 *   - IPW (Inverse Probability Weighting)
 *   - Doubly Robust estimation
 *   - Stratification
 *   - Matching (nearest neighbor)
 *   - Cohen's d and relative risk
 *   - Bonferroni correction
 *
 * Refs: Imbens & Rubin (2015), Hernan & Robins (2020).
 * ============================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfr_core.h"
#include "cfr_effects.h"
#include "cfr_potential.h"
#include "cfr_analysis.h"

int main(void) {
    printf("=== Example 2: Causal Effect Estimation ===\n\n");
    srand(99);
    int n = 20;

    /* 1. Generate synthetic data */
    double Y0[20], Y1[20], T[20], Y[20], covariates[40];
    for (int i = 0; i < n; i++) {
        Y0[i] = 50.0 + 5.0 * ((double)rand()/RAND_MAX);
        Y1[i] = Y0[i] + 8.0 + 2.0 * ((double)rand()/RAND_MAX);
        T[i] = (rand() % 2) ? 1.0 : 0.0;
        Y[i] = (T[i] > 0.5) ? Y1[i] : Y0[i];
        covariates[i*2] = (double)(i % 5);
        covariates[i*2+1] = Y0[i] / 10.0;
    }

    /* 2. Oracle ATE */
    CFRCausalEffects* eff = cfr_eff_create();
    cfr_eff_ate_oracle(eff, Y0, Y1, n);
    printf("Oracle ATE = %.3f [%.3f, %.3f]\n", eff->ate, eff->ci_lower, eff->ci_upper);

    /* 3. Difference in means */
    cfr_eff_ate_diff_means(eff, T, Y, n);
    printf("Diff-in-Means ATE = %.3f (SE=%.3f)\n", eff->ate, eff->se_ate);

    /* 4. IPW */
    double ps[20];
    for (int i=0;i<n;i++) ps[i]=0.5;
    cfr_eff_ate_ipw(eff, T, Y, ps, n);
    printf("IPW ATE = %.3f\n", eff->ate);

    /* 5. Doubly Robust */
    cfr_eff_ate_doubly_robust(eff, T, Y, ps, Y0, Y1, n);
    printf("Doubly Robust ATE = %.3f\n", eff->ate);

    /* 6. Stratification */
    int strata[20];
    for (int i=0;i<n;i++) strata[i] = (i < 10) ? 0 : 1;
    cfr_eff_ate_stratified(eff, T, Y, strata, n, 2);
    printf("Stratified ATE = %.3f\n", eff->ate);

    /* 7. Matching */
    cfr_eff_ate_matching(eff, T, Y, covariates, n, 2);
    printf("Matching ATE = %.3f\n", eff->ate);

    /* 8. Effect size */
    double d = cfr_eff_cohens_d(eff->ate, 5.0);
    printf("Cohen d = %.2f (%.1s effect)\n", d,
           d>0.8?"large":(d>0.5?"medium":"small"));

    /* 9. Relative risk (binary outcome) */
    double Y_bin[20];
    double median = 55.0;
    for (int i=0;i<n;i++) Y_bin[i] = (Y[i] > median) ? 1.0 : 0.0;
    double rr = cfr_eff_relative_risk(T, Y_bin, n);
    printf("Relative Risk = %.3f\n", rr);

    /* 10. Sensitivity analysis */
    double low, up;
    cfr_eff_sensitivity_analysis(eff->ate, 1.0, 1.5, &low, &up);
    printf("Sensitivity (gamma=1.5): [%.2f, %.2f]\n", low, up);

    /* 11. Bonferroni correction */
    double pvals[] = {0.001, 0.01, 0.03, 0.15, 0.50};
    int sig[5] = {0};
    cfr_eff_fwer_bonferroni(pvals, 5, 0.05, sig);
    printf("Bonferroni significant tests: %d/%d\n", sig[0]+sig[1]+sig[2]+sig[3]+sig[4], 5);

    cfr_eff_free(eff);
    printf("\n=== Example 2 PASSED ===\n");
    return 0;
}
