/*
 * scm_model.c �� Structural Causal Models (SCM)
 *
 * Implements Pearl's SCM framework: a triplet (U, V, F) where U are
 * exogenous variables, V are endogenous, and F are structural equations.
 * Supports linear equations, noisy-OR gates, do-intervention, and
 * Monte Carlo Average Causal Effect (ACE) estimation.
 *
 * References:
 *   Pearl, "Causality", 2009, Ch.7 (Counterfactuals and SCM)
 *   Pearl, Glymour, Jewell, "Causal Inference in Statistics", 2016, Ch.3
 *   Hernan & Robins, "Causal Inference: What If", 2020, Ch.2
 */

#include "scm_model.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * SCM Lifecycle
 * ================================================================ */

SCModel* scm_create(int n_vars) {
    if (n_vars <= 0 || n_vars > SCM_MAX_VARS) return NULL;
    SCModel* scm = calloc(1, sizeof(SCModel));
    if (!scm) return NULL;
    scm->n_vars = n_vars;
    scm->values = calloc((size_t)n_vars, sizeof(double));
    scm->equations = calloc((size_t)n_vars, sizeof(StructuralEq));
    scm->eq_n_parents = calloc((size_t)n_vars, sizeof(int));
    scm->noise = calloc((size_t)n_vars, sizeof(double));
    scm->dag = dag_create(n_vars);
    if (!scm->values || !scm->equations || !scm->eq_n_parents ||
        !scm->noise || !scm->dag) {
        scm_free(scm);
        return NULL;
    }
    return scm;
}

void scm_free(SCModel* scm) {
    if (!scm) return;
    dag_free(scm->dag);
    free(scm->values);
    free(scm->equations);
    free(scm->eq_n_parents);
    free(scm->noise);
    free(scm);
}

/* ================================================================
 * Structural Equation Management
 * ================================================================ */

int scm_set_equation(SCModel* scm, int i, StructuralEq eq, int n_parents) {
    if (!scm || i < 0 || i >= scm->n_vars || !eq || n_parents < 0)
        return -1;
    scm->equations[i] = eq;
    scm->eq_n_parents[i] = n_parents;
    return 0;
}

/*
 * Evaluate the structural equation for variable i.
 * The equation receives: parents[], n_parents, noise[i].
 */
double scm_evaluate(const SCModel* scm, int i) {
    if (!scm || i < 0 || i >= scm->n_vars || !scm->equations[i])
        return 0.0;

    int np = scm->eq_n_parents[i];
    double *parents_buf = calloc((size_t)(np > 0 ? np : 1), sizeof(double));
    if (!parents_buf) return 0.0;

    /* Collect parent values by scanning the DAG */
    int pcount = 0;
    for (int j = 0; j < scm->n_vars && pcount < np; j++) {
        if (dag_has_edge(scm->dag, j, i)) {
            parents_buf[pcount++] = scm->values[j];
        }
    }

    double result = scm->equations[i](parents_buf, np, scm->noise[i]);
    free(parents_buf);
    return result;
}

/* ================================================================
 * SCM Simulation �� evaluate all equations in topological order
 * ================================================================ */

int scm_simulate(SCModel* scm) {
    if (!scm) return -1;

    /* Ensure topological order is up-to-date */
    if (dag_topological_sort(scm->dag) != 0) return -1;

    for (int idx = 0; idx < scm->n_vars; idx++) {
        int i = scm->dag->topo_order[idx];
        if (i < 0) return -1;
        /* Skip intervened variables (equations set to NULL by do-intervention).
         * Their values are fixed externally and must not be overwritten. */
        if (scm->equations[i] == NULL) continue;
        scm->values[i] = scm_evaluate(scm, i);
    }
    return 0;
}

/* ================================================================
 * do-Intervention: set X = x, remove all incoming edges to X
 *
 * Creates a new SCM for the post-intervention world where the
 * structural equation of X is replaced by the constant x.
 * ================================================================ */

SCModel* scm_do_intervention(const SCModel* scm, int x, double value) {
    if (!scm || x < 0 || x >= scm->n_vars) return NULL;

    SCModel* scm_do = scm_create(scm->n_vars);
    if (!scm_do) return NULL;

    /* Copy DAG but remove edges into x */
    for (int i = 0; i < scm->n_vars; i++) {
        for (int j = 0; j < scm->n_vars; j++) {
            if (j == x) continue;  /* skip edges into x */
            if (dag_has_edge(scm->dag, i, j))
                dag_add_edge(scm_do->dag, i, j);
        }
    }

    /* Copy equations except for x */
    for (int i = 0; i < scm->n_vars; i++) {
        if (i == x) {
            /* X is fixed �� use a constant equation (defined below) */
            scm_set_equation(scm_do, i, NULL, 0);
            scm_do->values[i] = value;
        } else {
            scm_set_equation(scm_do, i, scm->equations[i],
                             scm->eq_n_parents[i]);
        }
    }

    /* Copy noise */
    memcpy(scm_do->noise, scm->noise, (size_t)scm->n_vars * sizeof(double));

    return scm_do;
}

/* ================================================================
 * Monte Carlo Average Causal Effect (ACE)
 *
 * ACE = E[Y | do(X=x1)] - E[Y | do(X=x0)]
 *
 * For linear SCMs with Gaussian noise, we sample noise repeatedly,
 * simulate the intervened model, and average the outcome.
 * ================================================================ */

/* Simple Box-Muller transform for Gaussian samples */
static double rand_gaussian(void) {
    double u1 = (double)rand() / (double)RAND_MAX;
    double u2 = (double)rand() / (double)RAND_MAX;
    /* Avoid log(0) */
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

double scm_average_causal_effect(SCModel* scm, int x, double x1, double x0,
                                  int y, int n_samples) {
    if (!scm || x < 0 || x >= scm->n_vars || y < 0 || y >= scm->n_vars ||
        n_samples <= 0) return 0.0;

    /* Create intervened SCMs */
    SCModel* scm_x1 = scm_do_intervention(scm, x, x1);
    SCModel* scm_x0 = scm_do_intervention(scm, x, x0);
    if (!scm_x1 || !scm_x0) {
        if (scm_x1) scm_free(scm_x1);
        if (scm_x0) scm_free(scm_x0);
        return 0.0;
    }

    /* For constant-intervention SCMs (no randomness), directly simulate */
    /* Check if the SCM has any non-constant equations for non-x variables */
    bool is_deterministic = true;
    for (int i = 0; i < scm->n_vars; i++) {
        if (i != x && scm->equations[i] && scm->eq_n_parents[i] > 0) {
            /* has parents �� might be deterministic or stochastic */
            is_deterministic = false;
            break;
        }
    }

    double sum_y1 = 0.0, sum_y0 = 0.0;

    if (is_deterministic) {
        /* Single evaluation */
        /* Set x1 */
        scm_x1->values[x] = x1;
        for (int idx = 1; idx < scm_x1->n_vars; idx++) {
            int i = scm_x1->dag->topo_order[idx];
            if (i >= 0 && i != x) {
                /* Simple deterministic: value = sum of parents */
                double s = 0.0;
                for (int p = 0; p < scm_x1->n_vars; p++)
                    if (dag_has_edge(scm_x1->dag, p, i))
                        s += scm_x1->values[p];
                scm_x1->values[i] = s + scm_x1->noise[i];
            }
        }
        sum_y1 = scm_x1->values[y];

        /* Set x0 */
        scm_x0->values[x] = x0;
        for (int idx = 1; idx < scm_x0->n_vars; idx++) {
            int i = scm_x0->dag->topo_order[idx];
            if (i >= 0 && i != x) {
                double s = 0.0;
                for (int p = 0; p < scm_x0->n_vars; p++)
                    if (dag_has_edge(scm_x0->dag, p, i))
                        s += scm_x0->values[p];
                scm_x0->values[i] = s + scm_x0->noise[i];
            }
        }
        sum_y0 = scm_x0->values[y];
    } else {
        /* Monte Carlo simulation */
        for (int s = 0; s < n_samples; s++) {
            /* Randomize noise */
            for (int i = 0; i < scm_x1->n_vars; i++) {
                if (i != x) {
                    scm_x1->noise[i] = rand_gaussian() * 0.5;
                    scm_x0->noise[i] = rand_gaussian() * 0.5;
                }
            }
            scm_x1->values[x] = x1;
            scm_x0->values[x] = x0;

            if (scm_simulate(scm_x1) == 0) sum_y1 += scm_x1->values[y];
            if (scm_simulate(scm_x0) == 0) sum_y0 += scm_x0->values[y];
        }
        sum_y1 /= (double)n_samples;
        sum_y0 /= (double)n_samples;
    }

    scm_free(scm_x1);
    scm_free(scm_x0);
    return sum_y1 - sum_y0;
}

/* ================================================================
 * Pre-built Structural Equations
 * ================================================================ */

/*
 * Linear equation: Y = sum_i beta_i * parent_i + noise
 * The coefficients are derived from the noise[] array for simplicity:
 * noise[0] encodes the intercept, or this is a simple linear combiner.
 *
 * For a general linear SCM: parents = [beta_0, beta_1, ..., beta_{k-1}]
 * Y = beta_0*X_0 + beta_1*X_1 + ... + beta_{k-1}*X_{k-1} + u
 *
 * Here we use a simple unit-coefficient version: Y = sum(parents) + u
 * The full coefficient version is provided via scm_linear_eq_coeff below.
 */
double scm_linear_eq(const double *parents, int np, double u) {
    double sum = 0.0;
    for (int i = 0; i < np; i++) sum += parents[i];
    return sum + u;
}

/*
 * Linear equation with explicit coefficients.
 * coeff array = [beta_0, ..., beta_{np-1}] stored in parents as:
 *   first np values = parent values
 * This version interprets each parent value as beta_i * X_i already.
 * For a simple sum: see scm_linear_eq above.
 */
double scm_linear_coeff_eq(const double *parents, int np, double u,
                            const double *coeff) {
    double sum = 0.0;
    for (int i = 0; i < np; i++) sum += coeff[i] * parents[i];
    return sum + u;
}

/*
 * Noisy-OR gate (common in epidemiology / Bayesian networks):
 *   P(Y=1 | X_1,...,X_k) = 1 - prod_i (1 - p_i)^{X_i}
 *
 * Here we implement the continuous relaxation:
 *   Y = 1 - prod_i (1 - parents[i]) + noise
 * where 0 <= parent[i] <= 1 represents probability of activation.
 */
double scm_noisy_or(const double *parents, int np, double u) {
    double prod = 1.0;
    bool any_active = false;
    for (int i = 0; i < np; i++) {
        double p = parents[i];
        if (p < 0.0) p = 0.0;
        if (p > 1.0) p = 1.0;
        if (p > 1e-8) {
            prod *= (1.0 - p);
            any_active = true;
        }
    }
    double base = any_active ? (1.0 - prod) : 0.0;
    /* Bound to [0, 1] */
    double result = base + u;
    if (result < 0.0) result = 0.0;
    if (result > 1.0) result = 1.0;
    return result;
}

/*
 * Threshold (sigmoid-based) structural equation:
 *   Y = sigma(beta_0 + sum_i beta_i * X_i + u) where sigma(z) = 1/(1+e^{-z})
 * Useful for binary outcomes.
 */
double scm_sigmoid_eq(const double *parents, int np, double u) {
    double z = u;  /* u acts as intercept */
    for (int i = 0; i < np; i++) z += parents[i];
    /* Clamp to avoid overflow */
    if (z > 50.0) return 1.0;
    if (z < -50.0) return 0.0;
    return 1.0 / (1.0 + exp(-z));
}

/*
 * Interaction term equation:
 *   Y = beta*X1*X2 + sum_i gamma_i * X_i + u
 * Captures effect modification / interaction.
 */
double scm_interaction_eq(const double *parents, int np, double u) {
    double sum = u;
    for (int i = 0; i < np; i++) sum += parents[i];
    /* Add pairwise interaction of first two parents if present */
    if (np >= 2) sum += 0.5 * parents[0] * parents[1];
    return sum;
}

/*
 * Polynomial equation: Y = a + b*X + c*X^2 + u
 * Captures non-linear dose-response relationships.
 */
double scm_polynomial_eq(const double *parents, int np, double u) {
    double sum = u;
    for (int i = 0; i < np; i++) {
        sum += parents[i] + 0.1 * parents[i] * parents[i];
    }
    return sum;
}
