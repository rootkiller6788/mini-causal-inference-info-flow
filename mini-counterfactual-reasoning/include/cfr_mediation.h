#ifndef CFR_MEDIATION_H
#define CFR_MEDIATION_H

#include "cfr_core.h"

/* ============================================================
 * cfr_mediation.h -- Causal Mediation Analysis
 *
 * Decomposes total effect into direct and indirect effects
 * through a mediator variable M.
 *
 *   TE = NDE + NIE (under certain assumptions)
 *
 * Key estimands:
 *   NDE = E[Y(1, M(0)) - Y(0, M(0))]  (Natural Direct Effect)
 *   NIE = E[Y(1, M(1)) - Y(1, M(0))]  (Natural Indirect Effect)
 *   TE  = E[Y(1) - Y(0)] = NDE(0) + NIE(1) = NDE(1) + NIE(0)
 *
 * Controlled Direct Effect:
 *   CDE(m) = E[Y(1,m) - Y(0,m)]
 *
 * Assumptions:
 *   1. No unmeasured T-Y confounding
 *   2. No unmeasured T-M confounding
 *   3. No unmeasured M-Y confounding
 *   4. No T-induced M-Y confounding
 *
 * Refs: Pearl (2001), VanderWeele (2015),
 *       Imai, Keele & Tingley (2010)
 * ============================================================ */

typedef struct {
    double nde0;     /* NDE at M(0) */
    double nde1;     /* NDE at M(1) */
    double nie0;     /* NIE at Y(0) */
    double nie1;     /* NIE at Y(1) */
    double te;       /* Total Effect */
    double cde[10];  /* Controlled DE at M=m */
    int n_cde;       /* number of CDE levels */
    double proportion_mediated; /* NIE/TE */
} CFRMediationEffects;

typedef struct {
    CFRSCM* scm;
    int treatment_id;
    int mediator_id;
    int outcome_id;
    double* covariates;
    int n_covariates;
} CFRMediationModel;

/* --- Lifecycle --- */
CFRMediationModel* cfr_med_model_create(CFRSCM* scm,
    int treatment_id, int mediator_id, int outcome_id);
void cfr_med_model_free(CFRMediationModel* mod);

/* --- Effect Computation --- */
CFRMediationEffects* cfr_med_effects_create(void);
void cfr_med_effects_free(CFRMediationEffects* eff);

/* Barron-Kenny approach (linear SCM) */
int cfr_med_barron_kenny(CFRMediationModel* mod, CFRMediationEffects* eff);

/* Pearl's mediation formula (non-parametric) */
int cfr_med_pearl_formula(CFRMediationModel* mod, CFRMediationEffects* eff);

/* Linear SEM method */
int cfr_med_linear_sem(CFRMediationModel* mod, CFRMediationEffects* eff);

/* Controlled Direct Effect for specified mediator levels */
int cfr_med_controlled_de(CFRMediationModel* mod,
    double* mediator_levels, int n_levels, CFRMediationEffects* eff);

/* --- Sensitivity Analysis --- */
void cfr_med_sensitivity_rho(CFRMediationModel* mod,
    double rho, CFRMediationEffects* eff);
double cfr_med_required_rho(CFRMediationModel* mod,
    double target_nie, double* achieved_rho);

/* --- Output --- */
void cfr_med_print(CFRMediationEffects* eff);
void cfr_med_print_model(CFRMediationModel* mod);

#endif /* CFR_MEDIATION_H */
