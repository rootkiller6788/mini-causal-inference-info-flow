#ifndef SCM_ADVANCED_H
#define SCM_ADVANCED_H
#include "scm_core.h"

/* ============================================================================
 * Advanced SCM Topics (L8)
 * References:
 *   Bareinboim & Pearl (2016) "Causal inference and the data-fusion problem"
 *   Hernan & Robins (2020) "Causal Inference: What If" (Part III)
 *   VanderWeele & Ding (2017) "Sensitivity analysis in observational research"
 *   Kusner et al. (2017) "Counterfactual Fairness" (NIPS)
 *   Robins (1986) "A new approach to causal inference in mortality studies"
 * ============================================================================ */

/* Transportability result: can we transport effect from source to target population? */
typedef struct {
    bool transportable;           /* Is the effect transportable? */
    double source_effect;         /* Causal effect in source population */
    double target_effect;         /* Estimated effect in target population */
    double selection_diagram_score; /* How many differences between populations */
    char   required_adjustment[256];/* Variables needed for transport formula */
} SCM_TransportResult;

/* Selection bias correction result */
typedef struct {
    double naive_effect;          /* Effect ignoring selection */
    double corrected_effect;      /* Effect after correction */
    double selection_ratio;       /* Ratio: corrected / naive */
    bool   selection_present;     /* Is selection bias detected? */
    double heckman_lambda;        /* Inverse Mills ratio */
} SCM_SelectionCorrection;

/* Sensitivity analysis result for unmeasured confounding */
typedef struct {
    double observed_effect;       /* Effect with measured covariates */
    double e_value;               /* Minimum strength of unmeasured confounder
                                     to explain away the effect */
    double lower_bound;           /* Lower bound of true effect */
    double upper_bound;           /* Upper bound of true effect */
    bool   robust;                /* Whether result is robust to confounding */
} SCM_SensitivityResult;

/* Counterfactual fairness result */
typedef struct {
    double total_effect;          /* Total effect of protected attribute */
    double direct_effect;         /* Direct (discriminatory) effect */
    double indirect_effect;       /* Indirect (explainable) effect */
    double fairness_violation;    /* Measure of unfairness */
    bool   is_fair;               /* Whether decision satisfies counterfactual fairness */
} SCM_FairnessResult;

/* Time-varying treatment / g-formula result */
typedef struct {
    double g_formula_estimate;    /* g-computation estimate */
    double ipw_estimate;          /* IPW estimate */
    double doubly_robust_estimate;/* DR estimate */
    int    n_time_points;         /* Number of time points */
} SCM_LongitudinalResult;

/* === Transportability (L8: Bareinboim & Pearl 2016) === */

/* Test if a causal effect can be transported from source to target population.
 * Uses selection diagrams: S-nodes point to variables whose distributions
 * differ between populations. */
SCM_TransportResult* scm_transportability_test(const SCM_Model* source_model,
                                                const SCM_Model* target_model,
                                                int treatment, int outcome);

/* Compute the transport formula: direct transport without adjustment.
 * P*(y|do(x)) = sum_z P(y|do(x), z) P*(z) */
double scm_direct_transport(const SCM_Model* source, const SCM_Model* target,
                             int treatment, int outcome);

/* Compute transport formula with S-admissibility.
 * When a set Z renders the effect transportable: P*(y|do(x)) = sum_z P(y|do(x),z) P*(z) */
double scm_s_admissible_transport(const SCM_Model* source, const SCM_Model* target,
                                   int treatment, int outcome, const SCM_VarSet* Z);

/* === Selection Bias Correction (L8: Heckman-type models) === */

/* Detect and correct for selection bias using Heckman two-step method.
 * selection_var: whether observation is observed (1) or missing (0).
 * Returns corrected treatment effect. */
SCM_SelectionCorrection* scm_selection_bias_correction(const double* data, int n_rows,
                                                         int n_cols, int treatment,
                                                         int outcome, int selection_var,
                                                         const int* covariates, int n_covariates);

/* Compute Inverse Mills Ratio: lambda = phi(z)/Phi(z) for selected,
 * lambda = -phi(z)/(1-Phi(z)) for non-selected.
 * where z is the linear predictor from the selection equation. */
double scm_inverse_mills_ratio(double z_score, bool selected);

/* === Counterfactual Fairness (L8: Kusner et al. 2017) === */

/* Assess counterfactual fairness of a decision.
 * protected_attr: sensitive attribute (e.g., race, gender).
 * decision: the decision variable.
 * If changing protected_attr (counterfactually) changes the decision,
 * the algorithm is not counterfactually fair. */
SCM_FairnessResult* scm_counterfactual_fairness(const SCM_Model* m, int protected_attr,
                                                  int decision, int outcome);

/* === Sensitivity Analysis (L8: VanderWeele & Ding 2017) === */

/* Compute E-value for unmeasured confounding.
 * E-value: minimum strength of association (risk ratio scale) that an
 * unmeasured confounder would need to have with both treatment and outcome
 * to fully explain away the observed effect. */
SCM_SensitivityResult* scm_sensitivity_analysis_evalue(const SCM_Model* m,
                                                         int treatment, int outcome,
                                                         double observed_rr);

/* Partial identification: compute Manski bounds on treatment effect.
 * Bounds: [E[Y(1)]_L - E[Y(0)]_U, E[Y(1)]_U - E[Y(0)]_L].
 * Uses no assumptions beyond bounded outcomes. */
void scm_manski_bounds(const double* data, int n_rows, int n_cols,
                        int treatment, int outcome,
                        double outcome_min, double outcome_max,
                        double* lower_bound, double* upper_bound);

/* === Time-Varying Treatments (L8: g-formula / MSM) === */

/* g-formula for time-varying treatments.
 * Simulates counterfactual outcomes under a specified treatment regime.
 * treatment_regime: array of length n_time_points with desired treatment values. */
SCM_LongitudinalResult* scm_g_formula(const double* data, int n_rows, int n_cols,
                                        const int* treatment_vars,
                                        const int* outcome_vars,
                                        const int* covariates,
                                        int n_time_points);

/* Marginal Structural Model with IPW for time-varying confounding.
 * Estimates stabilized IPW weights for each time point. */
SCM_LongitudinalResult* scm_marginal_structural_model(const double* data, int n_rows,
                                                        int n_cols,
                                                        const int* treatment_vars,
                                                        const int* outcome_vars,
                                                        const int* confounders,
                                                        int n_time_points);

/* Doubly robust estimator combining g-formula and IPW. */
SCM_LongitudinalResult* scm_doubly_robust_estimator(const double* data, int n_rows,
                                                      int n_cols,
                                                      const int* treatment_vars,
                                                      const int* outcome_vars,
                                                      const int* confounders,
                                                      int n_time_points);

/* === Utility === */

void scm_transport_result_free(SCM_TransportResult* tr);
void scm_transport_result_print(const SCM_TransportResult* tr);
void scm_selection_correction_free(SCM_SelectionCorrection* sc);
void scm_selection_correction_print(const SCM_SelectionCorrection* sc);
void scm_sensitivity_result_free(SCM_SensitivityResult* sr);
void scm_sensitivity_result_print(const SCM_SensitivityResult* sr);
void scm_fairness_result_free(SCM_FairnessResult* fr);
void scm_fairness_result_print(const SCM_FairnessResult* fr);
void scm_longitudinal_result_free(SCM_LongitudinalResult* lr);
void scm_longitudinal_result_print(const SCM_LongitudinalResult* lr);

#endif /* SCM_ADVANCED_H */
