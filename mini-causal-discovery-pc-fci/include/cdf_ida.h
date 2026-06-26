#ifndef CDF_IDA_H
#define CDF_IDA_H

/*
 * cdf_ida.h — IDA: Intervention-calculus when DAG is Absent
 *
 * IDA (Maathuis, Kalisch & Buhlmann, 2009; Maathuis, Colombo,
 * Kalisch & Buhlmann, 2010) estimates the causal effect of X on Y
 * from observational data when the true DAG is unknown.
 *
 * The PC algorithm outputs a CPDAG (Markov equivalence class).
 * IDA enumerates all valid DAGs within this equivalence class,
 * estimates the causal effect for each, and returns the multiset
 * of possible effects as bounds on the true effect.
 *
 * Causal effect definition (Pearl's do-calculus):
 *   θ_{X→Y} = E[Y | do(X = x+1)] - E[Y | do(X = x)]
 *
 * Under linear Gaussian models:
 *   Y = β_{YX}·X + Σ_{j∈pa(Y)\{X}} β_{Yj}·pa_j(Y) + ε_Y
 *   θ_{X→Y} = β_{YX|pa(Y)} (regression coefficient adjusted for parents)
 *
 * The set of possible causal effects is:
 *   Θ = { θ_{X→Y}(G) : G ∈ Markov equivalence class of CPDAG }
 *
 * Reference: Maathuis et al. (2009), "Estimating high-dimensional
 *   intervention effects from observational data". Annals of Stats.
 * Reference: Maathuis et al. (2010), "Predicting causal effects in
 *   large-scale systems". JRSS-B.
 */

#include <stddef.h>
#include <stdbool.h>
#include "cdf_core.h"
#include "cdf_pc.h"

/* ── IDA Result ────────────────────────────────────────────────────── */

/** A single causal effect estimate from one valid DAG. */
typedef struct {
    int     x;              /* cause variable */
    int     y;              /* effect variable */
    double  effect_coef;    /* regression coefficient beta_{YX|pa(Y)} */
    int     *parents_used;  /* parent set used for adjustment */
    int     n_parents;      /* size of parent set */
} CdfIDAEffect;

/** The multiset of possible causal effects across the equivalence class. */
typedef struct {
    CdfIDAEffect *effects;  /* [n_effects] array of effect estimates */
    int           n_effects;/* number of valid DAGs enumerated */
    double        min_effect; /* minimum estimated causal effect */
    double        max_effect; /* maximum estimated causal effect */
    double        median_effect; /* median across equivalence class */
    bool          all_same_sign; /* true if all effects have same sign */
} CdfIDAResult;

/* ── IDA API ───────────────────────────────────────────────────────── */

/**
 * Estimate the causal effect of X on Y using IDA.
 *
 * Given a learned CPDAG and the dataset, IDA:
 *   1. Enumerates all valid parent sets of Y consistent with the CPDAG
 *   2. For each parent set, computes the regression coefficient
 *      of X in {X} ∪ pa(Y) predicting Y
 *   3. Returns the multiset of effects
 *
 * Theorem (Maathuis et al., 2009):
 *   Under linear Gaussian SEM with equal error variances,
 *   the set Θ contains the true causal effect θ* with probability → 1
 *   as N → ∞, assuming PC consistently recovers the CPDAG.
 *
 * Complexity: O(|parents_Y| * p * N) per DAG; total depends on
 * number of valid parent configurations.
 *
 * @param cp  learned CPDAG (from PC algorithm)
 * @param ds  observational dataset
 * @param x   cause variable index
 * @param y   effect variable index
 * @return    IDA result (caller must cdf_ida_result_free)
 */
CdfIDAResult* cdf_ida_estimate(const CdfGraph *cp,
                                const CdfDataset *ds,
                                int x, int y);

/** Free IDA result. */
void cdf_ida_result_free(CdfIDAResult *res);

/** Print IDA result summary. */
void cdf_ida_print(const CdfIDAResult *res);

/* ── IDA Utilities ─────────────────────────────────────────────────── */

/**
 * Enumerate all valid parent sets of a node consistent with a CPDAG.
 * A parent set is valid if it does not violate the CPDAG edge constraints.
 *
 * @param cp       learned CPDAG
 * @param y        target node
 * @param parent_sets  [max_sets][p] output: enumerated parent sets
 * @param n_sets   output: number of valid parent sets found
 * @param max_sets maximum number of parent sets to enumerate
 */
void cdf_ida_enumerate_parent_sets(const CdfGraph *cp, int y,
                                    int **parent_sets, int *n_sets,
                                    int max_sets);

/**
 * Estimate the causal effect using the back-door criterion.
 * For the effect of X on Y, a set Z satisfies the back-door criterion
 * if Z blocks all back-door paths from X to Y and no node in Z is a
 * descendant of X.
 *
 * @param ds  dataset
 * @param cp  CPDAG
 * @param x   cause
 * @param y   effect
 * @param Z   adjustment set
 * @param nZ  size of Z
 * @return    causal effect estimate
 */
double cdf_ida_effect_estimate(const CdfDataset *ds, int x, int y,
                                const int *Z, int nZ);

/**
 * Back-door admissible set check: does Z satisfy the back-door
 * criterion for (X→Y) in the given graph?
 *
 * @param g   causal graph (DAG or CPDAG)
 * @param x   cause
 * @param y   effect
 * @param Z   candidate adjustment set
 * @param nZ  size of Z
 * @return    true if Z is back-door admissible
 */
bool cdf_ida_is_backdoor_admissible(const CdfGraph *g, int x, int y,
                                     const int *Z, int nZ);

/**
 * Front-door criterion check: does Z satisfy the front-door
 * criterion for the causal effect of X on Y?
 *
 * Z satisfies the front-door criterion for (X→Y) iff:
 *   1. Z intercepts ALL directed paths from X to Y
 *   2. There is no unblocked back-door path from X to Z
 *   3. All back-door paths from Z to Y are blocked by {X}
 *
 * This permits estimation when back-door is not available
 * (e.g., unobserved confounders on X-Y edge).
 *
 * @param g   DAG
 * @param x   cause
 * @param y   effect
 * @param Z   candidate set (typically a mediator)
 * @param nZ  size of Z
 * @return    true if Z is front-door admissible
 */
bool cdf_ida_is_frontdoor_admissible(const CdfGraph *g, int x, int y,
                                      const int *Z, int nZ);

#endif /* CDF_IDA_H */
