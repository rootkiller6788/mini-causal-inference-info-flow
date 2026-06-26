#ifndef ADJUSTMENT_SETS_H
#define ADJUSTMENT_SETS_H
#include "causal_identification.h"

/*
 * adjustment_sets.h — Finding Valid Covariate Adjustment Sets
 *
 * For estimating P(Y | do(X)), we need to find a set of covariates Z
 * such that conditioning on Z identifies the causal effect:
 *   P(Y | do(X)) = sum_z P(Y | X, Z=z) * P(Z=z)
 *
 * Algorithms: Pearl's back-door, Generalized Adjustment Criterion,
 * IDA (Intervention calculus when DAG is Absent),
 * Covariate selection for efficiency.
 */

typedef struct {
    int    *variables;   /* indices of adjustment variables      */
    int     n_vars;      /* number of adjustment variables       */
    double  efficiency;  /* asymptotic variance measure (lower better) */
    bool    is_valid;    /* whether this set satisfies criterion */
} AdjustmentSet;

/* Find all minimal adjustment sets for (X, Y).
 * Minimal = removing any variable breaks the criterion. */
int adjustment_find_minimal_sets(const DAG* dag, int x, int y,
                                  AdjustmentSet *sets, int max_sets);

/* Find the optimal (minimum asymptotic variance) adjustment set.
 * Among valid sets, choose the one with smallest variance. */
AdjustmentSet* adjustment_optimal_set(const DAG* dag, int x, int y,
                                       const double *data, int n, int d);

/* Generalized Adjustment Criterion (Perkovic et al., 2015):
 * Z is valid if: (1) Z blocks all "proper causal paths" and
 * (2) Z contains no "forbidden nodes" (descendants of certain paths). */
bool adjustment_generalized_criterion(const DAG* dag, int x, int y,
                                       const int *z, int n_z);

/* Find the parent adjustment set: Z = parents(X).
 * Always valid under the causal Markov condition. */
AdjustmentSet* adjustment_parent_set(const DAG* dag, int x);

/* Covariate selection via Lasso/regression for efficiency.
 * Given candidate set, select variables that reduce outcome variance. */
int adjustment_select_efficient(const DAG* dag, int x, int y,
                                 const double *data, int n, int d,
                                 int *selected, int max_sel);

/* Check if Z is a valid adjustment set for (X,Y) using the
 * generalized adjustment criterion. */
bool adjustment_is_valid(const DAG* dag, int x, int y,
                          const int *z, int n_z);

void adjustment_set_free(AdjustmentSet* as);

#define ADJ_MAX_SETS    20
#define ADJ_MAX_VARS    50

#endif
