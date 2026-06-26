#ifndef CAUSAL_IDENTIFICATION_H
#define CAUSAL_IDENTIFICATION_H
#include "scm_model.h"

/*
 * Causal Identification: determining whether a causal effect can be
 * estimated from observational data given a DAG.
 *
 * Key criteria:
 *   - Back-door criterion: adjust for Z that blocks all back-door paths
 *   - Front-door criterion: use mediator M to identify effect
 *   - do-calculus rules (Pearl, 1995)
 */

typedef struct {
    bool identifiable;    int method;
    int *adjustment_set;  int n_adjust;
    char *explanation;
} IdentificationResult;

#define ID_METHOD_NONE       0
#define ID_METHOD_BACKDOOR   1
#define ID_METHOD_FRONTDOOR  2
#define ID_METHOD_IV         3
#define ID_METHOD_DO_CALC    4

/* Check if Z satisfies the back-door criterion for (X, Y).
 * Z blocks all back-door paths from X to Y and contains no descendants of X. */
bool causal_backdoor_criterion(const DAG* dag, int x, int y,
                                const int *z, int n_z);

/* Check if M satisfies the front-door criterion for (X, Y).
 * 1. All directed paths from X to Y go through M.
 * 2. There is no unblocked back-door path from X to M.
 * 3. All back-door paths from M to Y are blocked by X. */
bool causal_frontdoor_criterion(const DAG* dag, int x, int y, int m);

/* Identify causal effect of X on Y using the DAG.
 * Returns identification result with recommended adjustment set. */
IdentificationResult* causal_identify(const DAG* dag, int x, int y);
void ident_result_free(IdentificationResult* r);

/* do-calculus Rule 2 (action/observation exchange):
 * P(y | do(x), z, w) = P(y | x, z, w) if Y is d-separated from X by Z, W
 * in the graph with incoming edges to X removed. */
bool do_calculus_rule2(const DAG* dag, int x, int y,
                        const int *z, int n_z, const int *w, int n_w);

/* Find all variables that satisfy the back-door criterion for (X,Y). */
int causal_find_backdoor_adjustment(const DAG* dag, int x, int y,
                                     int *adjustment_set, int max_n);

#endif
