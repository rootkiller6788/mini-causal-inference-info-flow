#ifndef CFR_POTENTIAL_H
#define CFR_POTENTIAL_H

#include "cfr_core.h"

/* ============================================================
 * cfr_potential.h -- Potential Outcomes Framework
 *
 * In the Rubin Causal Model (RCM), each unit i has potential
 * outcomes Y_i(0) and Y_i(1) under control and treatment.
 *
 * The observed outcome is: Y_i = T_i * Y_i(1) + (1-T_i) * Y_i(0)
 *
 * Key assumption: SUTVA (Stable Unit Treatment Value Assumption)
 * - No interference between units
 * - Single version of treatment
 *
 * Refs: Rubin (1974), Imbens & Rubin (2015)
 * ============================================================ */

typedef struct {
    double* Y0;          /* potential outcome under control */
    double* Y1;          /* potential outcome under treatment */
    double* covariates;  /* unit covariates (n_units x n_covars) */
    int n_units;
    int n_covariates;
    bool randomized;     /* treatment was randomized */
} CFRPotentialOutcomes;

typedef struct {
    double* treatment;   /* treatment assignment T_i */
    double* outcome;     /* observed outcome Y_i */
    int n_units;
} CFRDataset;

typedef struct {
    double* propensity;   /* propensity scores e(x) = P(T=1|X=x) */
    double* weights_ipw;  /* inverse probability weights */
    double* weights_att;  /* ATT weights */
    int n_units;
} CFRPropensity;

/* --- Lifecycle --- */
CFRPotentialOutcomes* cfr_po_create(int n_units, int n_covariates);
void                  cfr_po_free(CFRPotentialOutcomes* po);
void                  cfr_po_set_Y0(CFRPotentialOutcomes* po, int i, double v);
void                  cfr_po_set_Y1(CFRPotentialOutcomes* po, int i, double v);
void                  cfr_po_set_covariate(CFRPotentialOutcomes* po,
                                            int unit, int covar, double v);

/* --- Generate Data --- */
CFRDataset* cfr_po_generate_dataset(CFRPotentialOutcomes* po,
                                     double* treatment_prob);
void cfr_dataset_free(CFRDataset* ds);

/* --- Propensity Score --- */
CFRPropensity* cfr_po_propensity_logistic(CFRDataset* ds,
                                           CFRPotentialOutcomes* po);
void cfr_propensity_free(CFRPropensity* ps);
void cfr_po_ipw_weights(CFRPropensity* ps, CFRDataset* ds);

/* --- Individual Treatment Effect --- */
double cfr_po_ite(CFRPotentialOutcomes* po, int unit);
int    cfr_po_ite_all(CFRPotentialOutcomes* po, double* ite_array);
double cfr_po_ite_from_data(CFRPropensity* ps, CFRDataset* ds, int unit);

/* --- Bounds --- */
void cfr_po_manski_bounds(CFRDataset* ds, double* lower, double* upper);
void cfr_po_instrumental_variable_bounds(CFRDataset* ds,
                                          double* instrument,
                                          double* lower, double* upper);

/* --- Output --- */
void cfr_po_print(CFRPotentialOutcomes* po);
void cfr_po_print_ite(CFRPotentialOutcomes* po);

#endif /* CFR_POTENTIAL_H */

