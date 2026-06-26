#ifndef SCM_MODEL_H
#define SCM_MODEL_H
#include "dag_representation.h"

/*
 * scm_model.h — Structural Causal Models (SCM)
 *
 * An SCM is a triplet (U, V, F) where:
 *   U = exogenous (unobserved) variables
 *   V = endogenous (observed) variables
 *   F = set of structural equations: v_i = f_i(pa_i, u_i)
 *
 * The causal DAG is defined by: edge v_j -> v_i iff v_j in pa_i.
 */

typedef struct { int n; double *U; } ExogenousVars;
typedef double (*StructuralEq)(const double *parents, int n_parents, double u);

typedef struct {
    int             n_vars;       /* number of endogenous variables   */
    double         *values;       /* current values of V              */
    StructuralEq   *equations;    /* structural equation for each var */
    int            *eq_n_parents; /* number of parents for each eq    */
    double         *noise;        /* exogenous noise U                */
    DAG            *dag;          /* causal DAG                       */
} SCModel;

SCModel* scm_create(int n_vars);
void scm_free(SCModel* scm);
int scm_set_equation(SCModel* scm, int i, StructuralEq eq, int n_parents);
double scm_evaluate(const SCModel* scm, int i);

/* Simulate the SCM by evaluating all equations in topological order. */
int scm_simulate(SCModel* scm);

/* Intervention: set X = x and remove all incoming edges to X.
 * Returns a new SCM representing the post-intervention distribution. */
SCModel* scm_do_intervention(const SCModel* scm, int x, double value);

/* Compute the Average Causal Effect (ACE):
 *   ACE = E[Y | do(X=x1)] - E[Y | do(X=x0)]
 * via Monte Carlo simulation of the intervened SCM. */
double scm_average_causal_effect(SCModel* scm, int x, double x1, double x0,
                                  int y, int n_samples);

/* Linear SCM: Y = beta*X + gamma*Z + noise */
double scm_linear_eq(const double *parents, int np, double u);

/* Noisy-OR: P(Y=1) = 1 - prod_i (1 - p_i)^X_i */
double scm_noisy_or(const double *parents, int np, double u);

#define SCM_MAX_VARS 100
#define SCM_DEFAULT_SAMPLES 5000

#endif
