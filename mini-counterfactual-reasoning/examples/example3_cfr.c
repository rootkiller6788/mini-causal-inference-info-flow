/* ============================================================
 * Example 3: Mediation Analysis & Bounds
 *
 * Demonstrates:
 *   - Barron-Kenny mediation
 *   - Pearl mediation formula
 *   - Controlled direct effects
 *   - Mediation sensitivity
 *   - Probability bounds (PN, PS, PNS)
 *   - Experimental/non-experimental bounds
 *   - Balke-Pearl IV bounds
 *
 * Refs: Pearl (2001), VanderWeele (2015).
 * ============================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cfr_core.h"
#include "cfr_mediation.h"
#include "cfr_bounds.h"
#include "cfr_counterfactual.h"

int main(void) {
    printf("=== Example 3: Mediation & Bounds ===\n\n");

    /* --- Part A: Mediation Analysis --- */

    /* 1. Create mediation SCM and estimate effects */
    CFRSCM* scm = cfr_create_mediation_scm();
    CFRMediationModel* mod = cfr_med_model_create(scm, 0, 1, 2);
    CFRMediationEffects* meff = cfr_med_effects_create();

    /* Barron-Kenny */
    cfr_med_barron_kenny(mod, meff);
    printf("--- Mediation (Barron-Kenny) ---\n");
    printf("  TE  = %.3f\n", meff->te);
    printf("  NDE = %.3f\n", meff->nde0);
    printf("  NIE = %.3f\n", meff->nie0);
    printf("  Proportion mediated = %.1f%%\n", meff->proportion_mediated * 100);

    /* Pearl formula */
    cfr_med_pearl_formula(mod, meff);
    printf("--- Mediation (Pearl Formula) ---\n");
    printf("  NDE = %.3f, NIE = %.3f, TE = %.3f\n", meff->nde0, meff->nie0, meff->te);

    /* Controlled Direct Effect */
    double m_levels[] = {0.0, 0.5, 1.0};
    cfr_med_controlled_de(mod, m_levels, 3, meff);
    printf("  CDE at M=%.1f: %.3f\n", m_levels[1], meff->cde[1]);

    /* Sensitivity */
    cfr_med_sensitivity_rho(mod, 0.2, meff);
    printf("  NIE under rho=0.2: %.3f\n", meff->nie0);

    double required_rho;
    cfr_med_required_rho(mod, 0.5, &required_rho);
    printf("  Required rho to reduce NIE to 0.5: %.3f\n", required_rho);

    /* Cross-world mediation via counterfactuals */
    double nde, nie, te;
    cfr_cf_mediation_effects(scm, 0, 1, 2, &nde, &nie, &te);
    printf("  CF Mediation: NDE=%.3f NIE=%.3f TE=%.3f\n", nde, nie, te);

    /* --- Part B: Probability Bounds --- */

    printf("\n--- Counterfactual Bounds ---\n");
    CFRBoundsData* bdata = cfr_bounds_data_create(6);
    /* Simulated data: 3 treated (all success), 3 control (1 success) */
    bdata->treatments[0]=1; bdata->treatments[1]=1; bdata->treatments[2]=1;
    bdata->treatments[3]=0; bdata->treatments[4]=0; bdata->treatments[5]=0;
    bdata->outcomes[0]=1; bdata->outcomes[1]=1; bdata->outcomes[2]=1;
    bdata->outcomes[3]=0; bdata->outcomes[4]=0; bdata->outcomes[5]=1;
    bdata->randomized = true;

    CFRBounds* bnd = cfr_bounds_create();

    cfr_bounds_trivial(bdata, bnd);
    printf("  Trivial PN bounds:   [%.2f, %.2f]\n", bnd->pn_lower, bnd->pn_upper);

    cfr_bounds_experimental(bdata, bnd);
    printf("  Experimental PN bounds: [%.3f, %.3f]\n", bnd->pn_lower, bnd->pn_upper);
    printf("  Experimental PS bounds: [%.3f, %.3f]\n", bnd->ps_lower, bnd->ps_upper);

    /* Monotonicity bounds */
    CFRBoundsData* bdata_m = cfr_bounds_data_create(6);
    bdata_m->treatments[0]=1; bdata_m->treatments[1]=1; bdata_m->treatments[2]=1;
    bdata_m->treatments[3]=0; bdata_m->treatments[4]=0; bdata_m->treatments[5]=0;
    bdata_m->outcomes[0]=1; bdata_m->outcomes[1]=1; bdata_m->outcomes[2]=1;
    bdata_m->outcomes[3]=0; bdata_m->outcomes[4]=0; bdata_m->outcomes[5]=0;
    bdata_m->monotonicity = true;

    cfr_bounds_with_monotonicity(bdata_m, bnd);
    printf("  Monotonicity PN: %.3f (exact)\n", bnd->pn_lower);
    printf("  Identifiable: %s\n", bnd->identifiable ? "yes" : "no");

    /* Excess fraction */
    double af = cfr_bounds_excess_fraction(bdata);
    printf("  Excess fraction: %.3f\n", af);

    /* Individual bounds */
    double pn_b[2], ps_b[2];
    cfr_bounds_individual(1.0, 0.0, pn_b, ps_b);
    printf("  Individual PN: %.2f, PS: %.2f\n", pn_b[0], ps_b[0]);

    /* Balke-Pearl IV bounds */
    double iv_low, iv_up;
    cfr_bounds_balke_pearl_iv(0.8, 0.4, 0.2, 0.6, &iv_low, &iv_up);
    printf("  Balke-Pearl IV bounds: [%.3f, %.3f]\n", iv_low, iv_up);

    /* Cross-world counterfactual bounds */
    double pn_l, pn_u;
    cfr_bounds_cross_world_counterfactual(0.3, 0.2, &pn_l, &pn_u);
    printf("  Cross-world PN: [%.3f, %.3f]\n", pn_l, pn_u);

    /* Stochastic dominance */
    double y1d[]={2,3,4,5,6}, y0d[]={1,2,3,4,5};
    int dominates;
    cfr_bounds_stochastic_dominance(y1d, 5, y0d, 5, &dominates);
    printf("  Y(1) dominates Y(0): %s\n", dominates==1?"yes":(dominates==-1?"Y(0) dominates":"no"));

    /* Cleanup */
    cfr_scm_free(scm); cfr_med_model_free(mod); cfr_med_effects_free(meff);
    cfr_bounds_free(bnd); cfr_bounds_data_free(bdata); cfr_bounds_data_free(bdata_m);

    printf("\n=== Example 3 PASSED ===\n");
    return 0;
}
