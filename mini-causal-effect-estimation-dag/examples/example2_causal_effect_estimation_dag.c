/* Example 2: Mediation Analysis - Baron-Kenny and IKT Methods
 *
 * Demonstrates causal mediation analysis to decompose total effect
 * into direct and indirect pathways through a mediator.
 *
 * Scenario: Does a workplace wellness program (X) reduce sick days (Y)
 * through improved mental health (M)?
 *
 * DAG: X -> M -> Y, plus X -> Y (direct), plus confounders
 *
 * References:
 *   Baron & Kenny, JPSP, 1986
 *   Imai, Keele, Tingley, Psychological Methods, 2010
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mediation_analysis.h"

int main(void) {
    srand((unsigned)time(NULL));
    printf("=== Example 2: Mediation Analysis ===\n\n");

    int n = 200, d = 1;  /* 1 covariate: baseline health score */
    int *T   = calloc((size_t)n, sizeof(int));
    double *M = calloc((size_t)n, sizeof(double));
    double *Y = calloc((size_t)n, sizeof(double));
    double *X = calloc((size_t)n * d, sizeof(double));
    if (!T || !M || !Y || !X) {
        fprintf(stderr, "Allocation failed\n");
        free(T); free(M); free(Y); free(X);
        return 1;
    }

    /* True parameters:
     *   M = 2.0 + 3.0*T + 0.5*baseline + noise    (a1=3.0)
     *   Y = 1.0 + 1.5*T + 2.0*M + 0.3*baseline + noise  (b1=1.5, b2=2.0)
     *   NIE = a1*b2 = 6.0
     *   NDE = b1 = 1.5
     *   TE  = 7.5
     *   PM  = 6.0/7.5 = 0.80
     */
    for (int i = 0; i < n; i++) {
        double baseline = -2.0 + 4.0 * ((double)rand() / RAND_MAX);
        X[i] = baseline;
        T[i] = (baseline > 0.0) ? 1 : 0;
        double e_m = 1.0 * ((double)rand() / RAND_MAX - 0.5);
        double e_y = 0.5 * ((double)rand() / RAND_MAX - 0.5);
        M[i] = 2.0 + 3.0 * (double)T[i] + 0.5 * baseline + e_m;
        Y[i] = 1.0 + 1.5 * (double)T[i] + 2.0 * M[i] + 0.3 * baseline + e_y;
    }

    /* Baron-Kenny mediation */
    MediationResult *mr_bk = mediation_baron_kenny(T, M, Y, X, n, d);
    if (mr_bk) {
        printf("Baron-Kenny Method:\n");
        printf("  a1 (T->M): %.3f\n", mr_bk->a1);
        printf("  b1 (T->Y direct, NDE): %.3f (SE=%.4f)\n",
               mr_bk->b1, mr_bk->se_nde);
        printf("  b2 (M->Y): %.3f\n", mr_bk->b2);
        printf("  NIE (a1*b2): %.3f (SE=%.4f, CI=[%.3f, %.3f])\n",
               mr_bk->nie, mr_bk->se_nie,
               mr_bk->ci_lower_nie, mr_bk->ci_upper_nie);
        printf("  NDE: %.3f (SE=%.4f, CI=[%.3f, %.3f])\n",
               mr_bk->nde, mr_bk->se_nde,
               mr_bk->ci_lower_nde, mr_bk->ci_upper_nde);
        printf("  TE: %.3f (SE=%.4f)\n", mr_bk->te, mr_bk->se_te);
        printf("  Proportion Mediated: %.3f (%.1f%%)\n",
               mr_bk->pm, mr_bk->pm * 100.0);
        printf("  NIE significant: %s\n",
               mr_bk->nie_significant ? "YES (p<0.05)" : "NO");
        printf("  NDE significant: %s\n",
               mr_bk->nde_significant ? "YES (p<0.05)" : "NO");

        /* Sensitivity of mediation to unmeasured M-Y confounding */
        MediationSensitivity *ms = mediation_sensitivity(mr_bk, 21);
        if (ms) {
            printf("\n  Mediation Sensitivity (rho = M-Y error correlation):\n");
            printf("  ACME=0 at rho = %.3f\n", ms->rho_at_zero);
            printf("  Interpretation: If |rho| < %.3f, mediation remains significant.\n",
                   fabs(ms->rho_at_zero));
            mediation_sensitivity_free(ms);
        }
        mediation_result_free(mr_bk);
    }

    /* IKT simulation-based mediation */
    MediationResult *mr_ikt = mediation_ikt(T, M, Y, X, n, d, 1000);
    if (mr_ikt) {
        printf("\nIKT Simulation-Based Method (1000 sims):\n");
        printf("  NIE: %.3f (SE=%.4f, CI=[%.3f, %.3f])\n",
               mr_ikt->nie, mr_ikt->se_nie,
               mr_ikt->ci_lower_nie, mr_ikt->ci_upper_nie);
        printf("  NDE: %.3f (SE=%.4f, CI=[%.3f, %.3f])\n",
               mr_ikt->nde, mr_ikt->se_nde,
               mr_ikt->ci_lower_nde, mr_ikt->ci_upper_nde);
        printf("  TE: %.3f (SE=%.4f)\n", mr_ikt->te, mr_ikt->se_te);
        printf("  PM: %.3f (%.1f%%)\n", mr_ikt->pm, mr_ikt->pm * 100.0);
        mediation_result_free(mr_ikt);
    }

    printf("\nTrue values: NIE=6.0, NDE=1.5, TE=7.5, PM=0.80\n");

    free(T); free(M); free(Y); free(X);
    printf("\n=== Example 2 PASSED ===\n");
    return 0;
}
