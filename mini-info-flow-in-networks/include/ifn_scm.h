#ifndef IFN_SCM_H
#define IFN_SCM_H
#include "ifn_core.h"

/* Structural Causal Models (SCM): Pearl's framework for causal reasoning.
 *
 * An SCM is a 4-tuple (U, V, F, P(U)):
 *   U = exogenous variables (noise)
 *   V = endogenous variables
 *   F = structural equations: v_i = f_i(PA_i, u_i)
 *   P(U) = distribution over exogenous variables
 *
 * Key operations:
 *   - Intervention do(X=x): replace f_x with constant x
 *   - Counterfactual: what would Y be had X been x, given evidence e?
 *   - Identifiability: can causal effect be estimated from observational data?
 *
 * Ref: Pearl (2009) Causality, 2nd ed.
 *      Peters, Janzing & Schölkopf (2017) Elements of Causal Inference
 */

/* Additive Noise Model (ANM): Y = f(X) + N, N ⟂ X.
 * Causal direction identifiable when noise is independent of cause.
 * Returns confidence score for direction X->Y vs Y->X. */
double ifn_scm_additive_noise_test(const double* x, const double* y,
    int n, int n_bins);

/* Nonlinear ANM fitting with kernel ridge regression.
 * Fits f in Y = f(X) + N, tests independence of residuals from X. */
double ifn_scm_fit_anm_krr(const double* x, const double* y, int n,
    double lambda, double sigma, double* fitted);

/* Back-door adjustment: P(Y|do(X=x)) = sum_z P(Y|X=x,Z=z)P(Z=z)
 * when Z satisfies back-door criterion relative to (X,Y).
 * g: causal DAG, cause/sink: node indices, covariates: back-door set. */
double ifn_scm_backdoor_adjustment(const IFN_CausalGraph* g,
    int cause, int effect, const int* z, int n_z,
    const double* data, int n_samples);

/* Front-door adjustment using intermediate variable M:
 * P(Y|do(X=x)) = sum_m P(m|do(x))*sum_x' P(y|do(x'),m)P(x'). */
double ifn_scm_frontdoor_adjustment(const IFN_CausalGraph* g,
    int cause, int effect, int mediator,
    const double* data, int n_samples);

/* Counterfactual inference via three-step process (Pearl 2009):
 * 1. Abduction: update P(U) given evidence
 * 2. Action: apply do(X=x)
 * 3. Prediction: compute resulting distribution of Y */
void ifn_scm_counterfactual(const IFN_CausalGraph* g,
    int cause, int effect, double cause_value,
    const double* evidence, int n_evidence,
    double* counterfactual_mean, double* counterfactual_var);

/* LiNGAM (Linear Non-Gaussian Acyclic Model):
 * Discovers causal order using ICA on linear structural equations.
 * Returns causal ordering of variables. */
void ifn_lingam_discover(const double* data, int n_samples, int n_vars,
    int* causal_order);

#endif
