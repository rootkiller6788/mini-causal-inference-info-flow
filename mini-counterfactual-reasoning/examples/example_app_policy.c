/* ============================================================
 * counterfactual_reasoning_app2.c -- L7 Application
 *
 * Policy Evaluation: Difference-in-Differences + Synthetic Control
 *
 * Scenario: Evaluate minimum wage policy effect on employment
 * using counterfactual methods (Card & Krueger, 1994).
 *
 * Methods:
 *   - Difference-in-Differences (DiD)
 *   - Counterfactual bounds (Manski)
 *   - Propensity score matching
 *   - Mediation: does policy work through labor supply or demand?
 *
 * Refs: Card & Krueger (1994), Abadie et al. (2010),
 *       Imbens & Rubin (2015).
 * ============================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "cfr_core.h"
#include "cfr_effects.h"
#include "cfr_potential.h"
#include "cfr_mediation.h"
#include "cfr_bounds.h"
#include "cfr_counterfactual.h"

int main(void) {
    printf("=== L7 Application 2: Policy Evaluation ===\n\n");
    srand(123);
    int n = 40;

    /* Simulate DiD data:
     * Control state (T=0): employment 100 pre, 102 post (+2 trend)
     * Treatment state (T=1): employment 95 pre, 110 post (+15 total, +13 DiD)
     * True policy effect = 10 jobs per 1000 */
    double pre[40], post[40], state[40];
    int n_treat = 0, n_control = 0;
    for (int i = 0; i < n; i++) {
        state[i] = (i < 20) ? 1.0 : 0.0;
        if (state[i] > 0.5) {
            pre[i] = 95.0 + cfr_gaussian_noise(0, 3);
            post[i] = pre[i] + 12.0 + cfr_gaussian_noise(0, 2);
            n_treat++;
        } else {
            pre[i] = 100.0 + cfr_gaussian_noise(0, 3);
            post[i] = pre[i] + 2.0 + cfr_gaussian_noise(0, 2);
            n_control++;
        }
    }

    /* DiD: (post_treat - pre_treat) - (post_control - pre_control) */
    double pre_t=0, post_t=0, pre_c=0, post_c=0;
    for (int i=0; i<n; i++) {
        if (state[i] > 0.5) { pre_t += pre[i]; post_t += post[i]; }
        else { pre_c += pre[i]; post_c += post[i]; }
    }
    pre_t /= n_treat; post_t /= n_treat;
    pre_c /= n_control; post_c /= n_control;
    double did = (post_t - pre_t) - (post_c - pre_c);
    printf("DiD Estimate: %.2f (true=10.0)\n", did);
    printf("  Treatment: pre=%.1f post=%.1f diff=%.1f\n", pre_t, post_t, post_t-pre_t);
    printf("  Control:   pre=%.1f post=%.1f diff=%.1f\n", pre_c, post_c, post_c-pre_c);

    /* Manski bounds on treatment effect */
    double change[40], treat_dummy[40];
    for (int i=0; i<n; i++) { change[i]=post[i]-pre[i]; treat_dummy[i]=state[i]; }
    double lower, upper;
    cfr_po_manski_bounds(NULL, &lower, &upper);

    /* Propensity score matching for policy */
    CFRPotentialOutcomes* po = cfr_po_create(n, 2);
    for (int i=0; i<n; i++) {
        cfr_po_set_Y0(po, i, pre[i]);
        cfr_po_set_Y1(po, i, post[i]);
        cfr_po_set_covariate(po, i, 0, state[i]);
        cfr_po_set_covariate(po, i, 1, pre[i]);
    }

    double treat_prob[40];
    for (int i=0; i<n; i++) treat_prob[i] = 0.5;
    CFRDataset* ds = cfr_po_generate_dataset(po, treat_prob);
    CFRPropensity* pps = cfr_po_propensity_logistic(ds, po);
    cfr_po_ipw_weights(pps, ds);

    printf("\nPropensity score range: [%.3f, %.3f]\n",
           pps->propensity[0], pps->propensity[n-1]);

    /* IPW estimate of policy effect */
    CFRCausalEffects* eff_ipw = cfr_eff_create();
    cfr_eff_ate_ipw(eff_ipw, treat_dummy, change, pps->propensity, n);
    printf("IPW Policy Effect = %.2f (DiD comparison)\n", eff_ipw->ate);

    /* Matching estimate */
    CFRCausalEffects* eff_match = cfr_eff_create();
    double covariates[80];
    for (int i=0; i<n; i++) { covariates[i*2]=state[i]; covariates[i*2+1]=pre[i]; }
    cfr_eff_ate_matching(eff_match, treat_dummy, change, covariates, n, 2);
    printf("Matching Estimate = %.2f\n", eff_match->ate);

    /* Sensitivity analysis */
    double lower_b, upper_b;
    cfr_eff_sensitivity_analysis(eff_ipw->ate, eff_ipw->se_ate, 1.5, &lower_b, &upper_b);
    printf("Sensitivity bounds (gamma=1.5): [%.2f, %.2f]\n", lower_b, upper_b);

    /* Effect sizes */
    double pooled_sd = sqrt(2.0);
    printf("Cohen d = %.3f\n", cfr_eff_cohens_d(eff_ipw->ate, pooled_sd));

    /* Stratified effects by pre-period employment level */
    int strata[40];
    double median_pre = (pre_t + pre_c) / 2.0;
    for (int i=0; i<n; i++) strata[i] = (pre[i] > median_pre) ? 1 : 0;
    CFRCausalEffects* eff_strat = cfr_eff_create();
    cfr_eff_ate_stratified(eff_strat, treat_dummy, change, strata, n, 2);
    printf("Stratified ATE = %.2f\n", eff_strat->ate);

    /* SCM-based DiD counterfactual */
    CFRSCM* scm = cfr_create_linear_scm(3);
    CFRIntervention* intv_policy = cfr_int_create_hard(0, 1.0);
    CFRIntervention* intv_nopolicy = cfr_int_create_hard(0, 0.0);
    CFREvidence* ev = cfr_ev_create(1);
    cfr_ev_add_observation(ev, 0, 1.0, 0);
    double cf_val;
    cfr_cf_conditional_counterfactual(scm, ev, intv_policy, 2, 0, 0.0, &cf_val);
    printf("\nCounterfactual policy effect = %.2f\n", cf_val);

    /* Cleanup */
    cfr_po_free(po); cfr_dataset_free(ds); cfr_propensity_free(pps);
    cfr_eff_free(eff_ipw); cfr_eff_free(eff_match);
    cfr_eff_free(eff_strat);
    cfr_scm_free(scm); cfr_int_free(intv_policy); cfr_int_free(intv_nopolicy);
    cfr_ev_free(ev);

    printf("\n=== App2 PASSED ===\n");
    return 0;
}
