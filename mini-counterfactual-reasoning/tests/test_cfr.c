#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cfr_core.h"
#include "cfr_counterfactual.h"
#include "cfr_potential.h"
#include "cfr_effects.h"
#include "cfr_mediation.h"
#include "cfr_bounds.h"

#define EPS 1e-9
#define TOL 1e-4

int main(void) {
    printf("=== Counterfactual Reasoning Test Suite ===\n");

    printf("  [scm]  SCM create..."); fflush(stdout);
    CFRSCM* scm = cfr_create_simple_treatment_scm();
    assert(scm); assert(scm->n_vars == 4);
    assert(cfr_find_var_index(scm, "T") >= 0);
    cfr_scm_free(scm);
    printf(" OK\n"); fflush(stdout);

    printf("  [intv] Intervention..."); fflush(stdout);
    scm = cfr_create_simple_treatment_scm();
    CFRIntervention* intv = cfr_int_create_hard(0, 1.0);
    assert(intv); assert(intv->var_id == 0);
    CFRSCM* mod_scm = cfr_scm_apply_intervention(scm, intv);
    assert(mod_scm);
    cfr_scm_free(scm); cfr_scm_free(mod_scm); cfr_int_free(intv);
    printf(" OK\n"); fflush(stdout);

    printf("  [po]   Potential outcomes..."); fflush(stdout);
    CFRPotentialOutcomes* po = cfr_po_create(10, 2);
    assert(po); assert(po->n_units == 10);
    cfr_po_set_Y0(po,5,2.5); cfr_po_set_Y1(po,5,7.5);
    assert(fabs(cfr_po_ite(po, 5) - 5.0) < TOL);
    cfr_po_free(po);
    printf(" OK\n"); fflush(stdout);

    printf("  [eff]  Causal effects..."); fflush(stdout);
    double Y0[]={1,2,3,4,5}, Y1[]={2,3,4,5,6};
    CFRCausalEffects* eff = cfr_eff_create(); assert(eff);
    cfr_eff_ate_oracle(eff, Y0, Y1, 5);
    assert(fabs(eff->ate - 1.0) < TOL);
    cfr_eff_free(eff);
    printf(" OK\n"); fflush(stdout);

    printf("  [med]  Mediation..."); fflush(stdout);
    scm = cfr_create_mediation_scm();
    assert(scm);
    CFRMediationModel* mod = cfr_med_model_create(scm, 0, 1, 2);
    assert(mod_scm);
    CFRMediationEffects* meff = cfr_med_effects_create(); assert(meff);
    assert(cfr_med_barron_kenny(mod, meff) == 0);
    assert(fabs(meff->te - 3.5) < TOL);
    cfr_med_effects_free(meff); cfr_med_model_free(mod); cfr_scm_free(scm);
    printf(" OK\n"); fflush(stdout);

    printf("  [bnd]  Bounds..."); fflush(stdout);
    CFRBoundsData* data = cfr_bounds_data_create(4);
    data->treatments[0]=1; data->treatments[1]=1;
    data->treatments[2]=0; data->treatments[3]=0;
    data->outcomes[0]=1; data->outcomes[1]=1;
    data->outcomes[2]=0; data->outcomes[3]=1;
    CFRBounds* b = cfr_bounds_create(); assert(b);
    cfr_bounds_experimental(data, b);
    assert(b->pn_lower >= 0 && b->pn_upper <= 1.0);
    cfr_bounds_free(b); cfr_bounds_data_free(data);
    printf(" OK\n"); fflush(stdout);

    printf("\n=== All tests passed! ===\n");
    return 0;
}
