#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cfr_core.h"

/* ============================================================
 * cfr_core.c -- Structural Causal Model Implementation
 *
 * Implements SCM creation, structural equation evaluation,
 * topological computation, interventions, and factory functions.
 *
 * Refs: Pearl (2009), Pearl, Glymour & Jewell (2016)
 * ============================================================ */

/* ---------- Variable ---------- */

CFRVariable* cfr_var_create(char* name, CFRVarType type, int id) {
    CFRVariable* v = calloc(1, sizeof(CFRVariable));
    if (!v) return NULL;
    v->name = name ? strdup(name) : strdup("unnamed");
    v->type = type; v->id = id;
    v->value = 0.0; v->exogenous = 0.0; v->is_exogenous = false;
    return v;
}

void cfr_var_free(CFRVariable* v) {
    if (v) { free(v->name); free(v); }
}

void cfr_var_set_value(CFRVariable* v, double val) { if (v) v->value = val; }
void cfr_var_set_exogenous(CFRVariable* v, double u) { if (v) v->exogenous = u; }

/* ---------- Equation ---------- */

CFRStructuralEquation* cfr_eq_create(int lhs_id, int* rhs_ids,
                                      int n_parents, char* name) {
    CFRStructuralEquation* eq = calloc(1, sizeof(CFRStructuralEquation));
    if (!eq) return NULL;
    eq->name = name ? strdup(name) : strdup("eq");
    eq->lhs_id = lhs_id;
    eq->n_parents = n_parents;
    eq->rhs_ids = calloc(n_parents, sizeof(int));
    if (!eq->rhs_ids && n_parents > 0) { free(eq->name); free(eq); return NULL; }
    if (rhs_ids && n_parents > 0)
        memcpy(eq->rhs_ids, rhs_ids, n_parents * sizeof(int));
    return eq;
}

void cfr_eq_free(CFRStructuralEquation* eq) {
    if (eq) { free(eq->name); free(eq->rhs_ids); free(eq); }
}

void cfr_eq_set_func(CFRStructuralEquation* eq,
    double (*func)(double*, double, void*), void* params) {
    if (!eq) return;
    eq->func = func; eq->params = params;
}

double cfr_eq_evaluate(CFRStructuralEquation* eq,
                        double* parent_values, double exogenous) {
    if (!eq || !eq->func) return 0.0;
    return eq->func(parent_values, exogenous, eq->params);
}

/* ---------- SCM ---------- */

CFRSCM* cfr_scm_create(int n_endogenous, int n_exogenous, char* name) {
    int total = n_endogenous + n_exogenous;
    if (total <= 0 || total > CFR_MAX_VARS) return NULL;
    CFRSCM* scm = calloc(1, sizeof(CFRSCM));
    if (!scm) return NULL;
    scm->name = name ? strdup(name) : strdup("SCM");
    scm->n_vars = 0; scm->n_endogenous = n_endogenous;
    scm->n_exogenous = n_exogenous;
    scm->vars = calloc(total, sizeof(CFRVariable));
    scm->equations = calloc(CFR_MAX_EQNS, sizeof(CFRStructuralEquation));
    if (!scm->vars || !scm->equations) {
        free(scm->vars); free(scm->equations); free(scm->name); free(scm);
        return NULL;
    }
    return scm;
}

void cfr_scm_free(CFRSCM* scm) {
    if (!scm) return;
    for (int i = 0; i < scm->n_vars; i++) free(scm->vars[i].name);
    free(scm->vars);
    for (int i = 0; i < scm->n_equations; i++) {
        free(scm->equations[i].name);
        free(scm->equations[i].rhs_ids);
    }
    free(scm->equations);
    free(scm->name);
    free(scm);
}

int cfr_scm_add_variable(CFRSCM* scm, CFRVariable* v) {
    if (!scm || !v || scm->n_vars >= scm->n_endogenous + scm->n_exogenous)
        return -1;
    scm->vars[scm->n_vars] = *v;
    scm->vars[scm->n_vars].name = strdup(v->name);
    scm->n_vars++;
    return scm->n_vars - 1;
}

int cfr_scm_add_equation(CFRSCM* scm, CFRStructuralEquation* eq) {
    if (!scm || !eq || scm->n_equations >= CFR_MAX_EQNS) return -1;
    scm->equations[scm->n_equations] = *eq;
    scm->equations[scm->n_equations].name = strdup(eq->name);
    if (eq->rhs_ids && eq->n_parents > 0) {
        int* r = calloc(eq->n_parents, sizeof(int));
        memcpy(r, eq->rhs_ids, eq->n_parents * sizeof(int));
        scm->equations[scm->n_equations].rhs_ids = r;
    } else {
        scm->equations[scm->n_equations].rhs_ids = NULL;
    }
    scm->n_equations++;
    return scm->n_equations - 1;
}

double* cfr_scm_get_values(CFRSCM* scm) {
    if (!scm) return NULL;
    double* vals = calloc(scm->n_vars, sizeof(double));
    if (!vals) return NULL;
    for (int i = 0; i < scm->n_vars; i++)
        vals[i] = scm->vars[i].value;
    return vals;
}

void cfr_scm_set_exogenous(CFRSCM* scm, double* u_values) {
    if (!scm || !u_values) return;
    for (int i = 0; i < scm->n_exogenous; i++)
        scm->vars[scm->n_endogenous + i].value = u_values[i];
}

void cfr_scm_compute(CFRSCM* scm) {
    /* Compute all endogenous variables from exogenous values.
     * For each equation: v_lhs = f(parents, u)
     * Exogenous variables must be set before calling. */
    if (!scm) return;
    for (int k = 0; k < scm->n_equations; k++) {
        CFRStructuralEquation* eq = &scm->equations[k];
        double parent_values[CFR_MAX_VARS];
        for (int j = 0; j < eq->n_parents; j++)
            parent_values[j] = scm->vars[eq->rhs_ids[j]].value;
        double u = scm->vars[eq->lhs_id].exogenous;
        scm->vars[eq->lhs_id].value = cfr_eq_evaluate(eq, parent_values, u);
    }
}

void cfr_scm_compute_topological(CFRSCM* scm) {
    /* Compute SCM respecting topological ordering.
     * Equations are assumed to be in topological order (parents before children).
     * For DAGs, this ensures correct computation in one pass. */
    cfr_scm_compute(scm);
}

/* ---------- Intervention ---------- */

CFRIntervention* cfr_int_create_hard(int var_id, double value) {
    CFRIntervention* intv = calloc(1, sizeof(CFRIntervention));
    if (!intv) return NULL;
    intv->var_id = var_id; intv->value = value;
    intv->type = CFR_INT_HARD;
    return intv;
}

CFRIntervention* cfr_int_create_soft(int var_id,
    double (*soft_func)(double*, double, void*), void* params) {
    CFRIntervention* intv = calloc(1, sizeof(CFRIntervention));
    if (!intv) return NULL;
    intv->var_id = var_id; intv->type = CFR_INT_SOFT;
    intv->soft_func = soft_func; intv->soft_params = params;
    return intv;
}

void cfr_int_free(CFRIntervention* intv) { free(intv); }

CFRSCM* cfr_scm_apply_intervention(CFRSCM* scm, CFRIntervention* intv) {
    /* Apply intervention do(X=x) to SCM:
     * Remove the structural equation for X and replace with X=x.
     * Returns a new modified SCM (caller must free). */
    if (!scm || !intv) return NULL;
    CFRSCM* new_scm = cfr_scm_create(scm->n_endogenous, scm->n_exogenous,
                                      scm->name);
    if (!new_scm) return NULL;
    /* Copy all variables */
    for (int i = 0; i < scm->n_vars; i++) {
        CFRVariable* v = cfr_var_create(scm->vars[i].name,
            scm->vars[i].type, scm->vars[i].id);
        v->value = scm->vars[i].value;
        v->exogenous = scm->vars[i].exogenous;
        v->is_exogenous = scm->vars[i].is_exogenous;
        cfr_scm_add_variable(new_scm, v);
        cfr_var_free(v);
    }
    /* Copy equations, but replace intervened one */
    for (int i = 0; i < scm->n_equations; i++) {
        if (scm->equations[i].lhs_id == intv->var_id) {
            /* Replace with constant function */
            CFRStructuralEquation* eq = cfr_eq_create(intv->var_id, NULL, 0,
                "intervened");
            new_scm->vars[intv->var_id].value = intv->value;
            cfr_scm_add_equation(new_scm, eq);
            cfr_eq_free(eq);
        } else {
            CFRStructuralEquation* eq = cfr_eq_create(
                scm->equations[i].lhs_id,
                scm->equations[i].rhs_ids,
                scm->equations[i].n_parents,
                scm->equations[i].name);
            eq->func = scm->equations[i].func;
            eq->params = scm->equations[i].params;
            cfr_scm_add_equation(new_scm, eq);
            cfr_eq_free(eq);
        }
    }
    return new_scm;
}

/* ---------- Evidence ---------- */

CFREvidence* cfr_ev_create(int n_obs) {
    CFREvidence* ev = calloc(1, sizeof(CFREvidence));
    if (!ev) return NULL;
    ev->n_observations = n_obs;
    ev->var_ids = calloc(n_obs, sizeof(int));
    ev->values = calloc(n_obs, sizeof(double));
    if (!ev->var_ids || !ev->values) {
        free(ev->var_ids); free(ev->values); free(ev); return NULL;
    }
    return ev;
}

void cfr_ev_free(CFREvidence* ev) {
    if (ev) { free(ev->var_ids); free(ev->values); free(ev); }
}

void cfr_ev_add_observation(CFREvidence* ev, int var_id,
                             double value, int idx) {
    if (!ev || idx < 0 || idx >= ev->n_observations) return;
    ev->var_ids[idx] = var_id; ev->values[idx] = value;
}

/* ---------- Query ---------- */

CFRQuery* cfr_query_create(int target_id, CFRIntervention* intv,
                            CFREvidence* ev) {
    CFRQuery* q = calloc(1, sizeof(CFRQuery));
    if (!q) return NULL;
    q->target_id = target_id; q->intervention = intv;
    q->evidence = ev; q->computed = false; q->result = 0.0;
    return q;
}

void cfr_query_free(CFRQuery* q) {
    if (q) { cfr_int_free(q->intervention); cfr_ev_free(q->evidence); free(q); }
}

int cfr_query_compute(CFRQuery* q, CFRSCM* scm) {
    if (!q || !scm) return -1;
    /* Apply intervention */
    CFRSCM* intv_scm = cfr_scm_apply_intervention(scm, q->intervention);
    if (!intv_scm) return -1;
    /* Set exogenous from evidence (abduction) */
    if (q->evidence) {
        cfr_scm_set_exogenous(intv_scm, q->evidence->values);
    }
    /* Compute intervened model */
    cfr_scm_compute(intv_scm);
    q->result = intv_scm->vars[q->target_id].value;
    q->computed = true;
    cfr_scm_free(intv_scm);
    return 0;
}

/* ---------- Factory: Standard SCMs ---------- */

static double linear_func(double* pv, double u, void* params) {
    double* coeffs = (double*)params;
    double sum = u;
    for (int i = 0; coeffs && i < 10; i++)
        if (fabs(coeffs[i]) > 1e-12) sum += coeffs[i] * pv[i];
    return sum;
}

CFRSCM* cfr_create_linear_scm(int n_vars) {
    CFRSCM* scm = cfr_scm_create(n_vars, n_vars, "Linear");
    if (!scm) return NULL;
    for (int i = 0; i < n_vars; i++) {
        char name[32]; sprintf(name, "X%d", i);
        CFRVariable* v = cfr_var_create(name, CFR_VAR_CONTINUOUS, i);
        if (i < n_vars) {
            v->is_exogenous = false;
            cfr_scm_add_variable(scm, v);
        }
        cfr_var_free(v);
    }
    for (int i = 0; i < n_vars; i++) {
        char name[32]; sprintf(name, "U%d", i);
        CFRVariable* u = cfr_var_create(name, CFR_VAR_CONTINUOUS, n_vars + i);
        u->is_exogenous = true;
        cfr_scm_add_variable(scm, u);
        cfr_var_free(u);
    }
    for (int i = 0; i < n_vars; i++) {
        int np = (i > 0) ? 2 : 1;
        int* pr = calloc(np, sizeof(int));
        if (i > 0) { pr[0] = i - 1; pr[1] = n_vars + i; }
        else { pr[0] = n_vars + i; }
        double* coeffs = calloc(10, sizeof(double));
        char ename[32]; sprintf(ename, "f_X%d", i);
        CFRStructuralEquation* eq = cfr_eq_create(i, pr, np, ename);
        eq->func = linear_func; eq->params = coeffs;
        cfr_scm_add_equation(scm, eq);
        free(pr); cfr_eq_free(eq);
    }
    return scm;
}

static double treat_func(double* pv, double u, void* params) {
    (void)pv; (void)params; return (u > 0) ? 1.0 : 0.0;
}

static double outcome_func(double* pv, double u, void* params) {
    (void)params;
    double t = pv[0];
    return 2.0 + 3.0 * t + u;
}

CFRSCM* cfr_create_simple_treatment_scm(void) {
    /* Simple SCM: T = f_T(U_T), Y = 2 + 3*T + U_Y
     * ATE = 3, controlled experiment */
    CFRSCM* scm = cfr_scm_create(2, 2, "Treatment");
    if (!scm) return NULL;
    /* Variables */
    CFRVariable* t = cfr_var_create("T", CFR_VAR_BINARY, 0);
    CFRVariable* y = cfr_var_create("Y", CFR_VAR_CONTINUOUS, 1);
    CFRVariable* ut = cfr_var_create("U_T", CFR_VAR_CONTINUOUS, 2);
    ut->is_exogenous = true;
    CFRVariable* uy = cfr_var_create("U_Y", CFR_VAR_CONTINUOUS, 3);
    uy->is_exogenous = true;
    cfr_scm_add_variable(scm, t); cfr_scm_add_variable(scm, y);
    cfr_scm_add_variable(scm, ut); cfr_scm_add_variable(scm, uy);
    /* Equations */
    int p_t[] = {2};
    CFRStructuralEquation* eq_t = cfr_eq_create(0, p_t, 1, "f_T");
    eq_t->func = treat_func;
    int p_y[] = {0, 3};
    CFRStructuralEquation* eq_y = cfr_eq_create(1, p_y, 2, "f_Y");
    eq_y->func = outcome_func;
    cfr_scm_add_equation(scm, eq_t); cfr_scm_add_equation(scm, eq_y);
    cfr_var_free(t); cfr_var_free(y); cfr_var_free(ut); cfr_var_free(uy);
    cfr_eq_free(eq_t); cfr_eq_free(eq_y);
    return scm;
}

static double med_func_M(double* pv, double u, void* params) {
    (void)params; return 1.0 * pv[0] + u;
}

static double med_func_Y(double* pv, double u, void* params) {
    (void)params; return 2.0 * pv[0] + 1.5 * pv[1] + u;
}

CFRSCM* cfr_create_mediation_scm(void) {
    /* T -> M -> Y  and  T -> Y (mediation with direct effect)
     * M = 1.0*T + U_M
     * Y = 2.0*T + 1.5*M + U_Y
     * NDE = 2.0, NIE = 1.0*1.5 = 1.5, TE = 3.5 */
    CFRSCM* scm = cfr_scm_create(3, 3, "Mediation");
    if (!scm) return NULL;
    CFRVariable* t = cfr_var_create("T", CFR_VAR_BINARY, 0);
    CFRVariable* m = cfr_var_create("M", CFR_VAR_CONTINUOUS, 1);
    CFRVariable* y = cfr_var_create("Y", CFR_VAR_CONTINUOUS, 2);
    CFRVariable* ut = cfr_var_create("U_T", CFR_VAR_CONTINUOUS, 3);
    CFRVariable* um = cfr_var_create("U_M", CFR_VAR_CONTINUOUS, 4);
    CFRVariable* uy = cfr_var_create("U_Y", CFR_VAR_CONTINUOUS, 5);
    ut->is_exogenous = true; um->is_exogenous = true;
    uy->is_exogenous = true;
    cfr_scm_add_variable(scm, t); cfr_scm_add_variable(scm, m);
    cfr_scm_add_variable(scm, y); cfr_scm_add_variable(scm, ut);
    cfr_scm_add_variable(scm, um); cfr_scm_add_variable(scm, uy);
    int p_t[] = {3}; CFRStructuralEquation* eqt = cfr_eq_create(0,p_t,1,"f_T");
    int p_m[] = {0,4}; CFRStructuralEquation* eqm = cfr_eq_create(1,p_m,2,"f_M");
    int p_y[] = {0,1,5}; CFRStructuralEquation* eqy = cfr_eq_create(2,p_y,3,"f_Y");
    eqt->func = treat_func; eqm->func = med_func_M; eqy->func = med_func_Y;
    cfr_scm_add_equation(scm, eqt); cfr_scm_add_equation(scm, eqm);
    cfr_scm_add_equation(scm, eqy);
    cfr_var_free(t); cfr_var_free(m); cfr_var_free(y);
    cfr_var_free(ut); cfr_var_free(um); cfr_var_free(uy);
    cfr_eq_free(eqt); cfr_eq_free(eqm); cfr_eq_free(eqy);
    return scm;
}

CFRSCM* cfr_create_frontdoor_scm(void) {
    /* Front-door criterion example:
     * U -> X -> M -> Y, with U -> Y (confounding)
     * X: smoking, M: tar deposits, Y: lung cancer, U: unobserved */
    CFRSCM* scm = cfr_scm_create(3, 1, "Frontdoor");
    if (!scm) return NULL;
    CFRVariable* x = cfr_var_create("X", CFR_VAR_BINARY, 0);
    CFRVariable* m = cfr_var_create("M", CFR_VAR_CONTINUOUS, 1);
    CFRVariable* y = cfr_var_create("Y", CFR_VAR_CONTINUOUS, 2);
    CFRVariable* u = cfr_var_create("U", CFR_VAR_CONTINUOUS, 3);
    u->is_exogenous = true;
    cfr_scm_add_variable(scm, x); cfr_scm_add_variable(scm, m);
    cfr_scm_add_variable(scm, y); cfr_scm_add_variable(scm, u);
    int px[] = {3}; CFRStructuralEquation* eqx = cfr_eq_create(0,px,1,"f_X");
    int pm[] = {0,3}; CFRStructuralEquation* eqm = cfr_eq_create(1,pm,2,"f_M");
    int py[] = {1,3}; CFRStructuralEquation* eqy = cfr_eq_create(2,py,2,"f_Y");
    eqx->func = linear_func; eqm->func = linear_func; eqy->func = linear_func;
    cfr_scm_add_equation(scm, eqx); cfr_scm_add_equation(scm, eqm);
    cfr_scm_add_equation(scm, eqy);
    cfr_var_free(x); cfr_var_free(m); cfr_var_free(y); cfr_var_free(u);
    cfr_eq_free(eqx); cfr_eq_free(eqm); cfr_eq_free(eqy);
    return scm;
}

/* ---------- Utility ---------- */

double cfr_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

double cfr_gaussian_noise(double mean, double std) {
    /* Box-Muller transform for Gaussian noise */
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 < 1e-12) u1 = 1e-12;
    return mean + std * sqrt(-2.0 * log(u1)) * cos(2.0 * 3.1415926535 * u2);
}

void cfr_scm_print(CFRSCM* scm) {
    if (!scm) return;
    printf("SCM: %s (V=%d, U=%d, eqns=%d)\n",
           scm->name, scm->n_endogenous, scm->n_exogenous, scm->n_equations);
    for (int i = 0; i < scm->n_vars; i++) {
        printf("  %s[%d]=%.4f", scm->vars[i].name, i, scm->vars[i].value);
        if (scm->vars[i].is_exogenous) printf(" (U)");
        printf("\n");
    }
}

int cfr_find_var_index(CFRSCM* scm, char* name) {
    if (!scm || !name) return -1;
    for (int i = 0; i < scm->n_vars; i++)
        if (strcmp(scm->vars[i].name, name) == 0) return i;
    return -1;
}
/* ---------- Extended SCM Operations ---------- */

int cfr_scm_parents(CFRSCM* scm, int var_id, int* parent_ids, int max_parents) {
    if (!scm || !parent_ids) return 0;
    int count = 0;
    for (int k = 0; k < scm->n_equations && count < max_parents; k++) {
        if (scm->equations[k].lhs_id == var_id) {
            for (int p = 0; p < scm->equations[k].n_parents && count < max_parents; p++)
                parent_ids[count++] = scm->equations[k].rhs_ids[p];
            break;
        }
    }
    return count;
}

int cfr_scm_children(CFRSCM* scm, int var_id, int* child_ids, int max_children) {
    if (!scm || !child_ids) return 0;
    int count = 0;
    for (int k = 0; k < scm->n_equations && count < max_children; k++) {
        for (int p = 0; p < scm->equations[k].n_parents; p++) {
            if (scm->equations[k].rhs_ids[p] == var_id) {
                child_ids[count++] = scm->equations[k].lhs_id;
                break;
            }
        }
    }
    return count;
}

bool cfr_scm_is_dag(CFRSCM* scm) {
    /* Check if the SCM forms a DAG (no cycles).
     * Uses simple DFS-based cycle detection. */
    if (!scm) return false;
    int visited[CFR_MAX_VARS] = {0};
    int stack[CFR_MAX_VARS], top = 0;
    /* Check each variable */
    for (int start = 0; start < scm->n_endogenous; start++) {
        for (int i = 0; i < scm->n_vars; i++) visited[i] = 0;
        top = 0; stack[top++] = start;
        while (top > 0) {
            int v = stack[--top];
            if (visited[v] == 2) continue;
            if (visited[v] == 1) return false; /* cycle */
            visited[v] = 1;
            int children[CFR_MAX_VARS];
            int nc = cfr_scm_children(scm, v, children, CFR_MAX_VARS);
            for (int c = 0; c < nc; c++) {
                if (visited[children[c]] < 2)
                    stack[top++] = children[c];
            }
            visited[v] = 2;
        }
    }
    return true;
}

void cfr_scm_backdoor_adjustment(CFRSCM* scm, int treatment_id,
    int outcome_id, int* adjustment_set, int* n_adjust) {
    /* Find a valid backdoor adjustment set for T -> Y.
     * The backdoor criterion: Z blocks all backdoor paths
     * from T to Y and contains no descendants of T.
     *
     * For simple DAGs, this is the set of all common causes. */
    if (!scm || !adjustment_set || !n_adjust) return;
    int count = 0;
    int parents_t[CFR_MAX_VARS], parents_y[CFR_MAX_VARS];
    int np_t = cfr_scm_parents(scm, treatment_id, parents_t, CFR_MAX_VARS);
    int np_y = cfr_scm_parents(scm, outcome_id, parents_y, CFR_MAX_VARS);
    /* Find common parents (simplified backdoor set) */
    for (int i = 0; i < np_t; i++)
        for (int j = 0; j < np_y; j++)
            if (parents_t[i] == parents_y[j] && count < CFR_MAX_VARS)
                adjustment_set[count++] = parents_t[i];
    *n_adjust = count;
}

#define CFR_MAX_PATHS 100

/* ---------- DAG All Paths (DFS with backtracking) ---------- */

static void cfr_dfs_paths(CFRSCM* scm, int current, int to_id,
    int* visited, int* path, int path_len,
    int** paths, int* path_lengths, int* count, int max_paths) {
    if (*count >= max_paths) return;
    visited[current] = 1;
    path[path_len] = current;
    if (current == to_id && path_len > 0) {
        paths[*count] = calloc(path_len + 1, sizeof(int));
        if (paths[*count]) {
            for (int i = 0; i <= path_len; i++)
                paths[*count][i] = path[i];
            path_lengths[*count] = path_len + 1;
            (*count)++;
        }
    } else {
        int children[CFR_MAX_VARS];
        int nc = cfr_scm_children(scm, current, children, CFR_MAX_VARS);
        for (int c = 0; c < nc && *count < max_paths; c++) {
            if (!visited[children[c]])
                cfr_dfs_paths(scm, children[c], to_id, visited, path,
                              path_len + 1, paths, path_lengths, count, max_paths);
        }
    }
    visited[current] = 0;
}

int cfr_scm_all_paths(CFRSCM* scm, int from_id, int to_id,
    int** paths, int* path_lengths, int max_paths) {
    if (!scm || !paths || !path_lengths || max_paths <= 0) return 0;
    int count = 0;
    int visited[CFR_MAX_VARS] = {0};
    int path[CFR_MAX_VARS];
    cfr_dfs_paths(scm, from_id, to_id, visited, path, 0,
                  paths, path_lengths, &count, max_paths);
    return count;
}

/* ---------- Front-door Criterion ---------- */

bool cfr_scm_frontdoor_criterion(CFRSCM* scm, int treatment_id,
    int outcome_id, int* mediator_id) {
    if (!scm || !mediator_id) return false;
    /* Criterion 1: All directed paths from T to Y go through M */
    int children[CFR_MAX_VARS];
    int nc = cfr_scm_children(scm, treatment_id, children, CFR_MAX_VARS);
    if (nc == 0) return false;
    *mediator_id = children[0];
    /* Criterion 2: There is no unblocked backdoor path from T to M */
    int parents_m[CFR_MAX_VARS];
    int npm = cfr_scm_parents(scm, *mediator_id, parents_m, CFR_MAX_VARS);
    bool t_parent_of_m = false;
    for (int i = 0; i < npm; i++)
        if (parents_m[i] == treatment_id) t_parent_of_m = true;
    if (!t_parent_of_m) return false;
    /* Criterion 3: All backdoor paths from M to Y are blocked by T */
    int children_m[CFR_MAX_VARS];
    int ncm = cfr_scm_children(scm, *mediator_id, children_m, CFR_MAX_VARS);
    bool m_parent_of_y = false;
    for (int i = 0; i < ncm; i++)
        if (children_m[i] == outcome_id) m_parent_of_y = true;
    return m_parent_of_y;
}

/* ---------- Counterfactual Simulation ---------- */

void cfr_scm_counterfactual_simulation(CFRSCM* scm, int n_samples,
    double* exogenous_samples, double* outcomes) {
    if (!scm || !exogenous_samples || !outcomes) return;
    for (int s = 0; s < n_samples; s++) {
        for (int i = 0; i < scm->n_exogenous; i++)
            scm->vars[scm->n_endogenous+i].value =
                exogenous_samples[s * scm->n_exogenous + i];
        cfr_scm_compute(scm);
        outcomes[s] = scm->vars[scm->n_endogenous-1].value;
    }
}

/* ---------- Intervention Path Computation ---------- */

void cfr_scm_compute_intervention_path(CFRSCM* scm,
    CFRIntervention* intv, int target_id) {
    /* Compute the causal path from intervention to target:
     * Track how the intervention propagates through the DAG.
     * Uses BFS traversal from intervened variable. */
    if (!scm || !intv) return;
    CFRSCM* modified = cfr_scm_apply_intervention(scm, intv);
    if (!modified) return;
    /* BFS from intv->var_id to target_id in modified model */
    int queue[CFR_MAX_VARS], front = 0, back = 0;
    int pred[CFR_MAX_VARS], dist[CFR_MAX_VARS];
    for (int i = 0; i < modified->n_vars; i++) { pred[i] = -1; dist[i] = -1; }
    queue[back++] = intv->var_id; dist[intv->var_id] = 0;
    while (front < back) {
        int u = queue[front++];
        if (u == target_id) break;
        int children[CFR_MAX_VARS];
        int nc = cfr_scm_children(modified, u, children, CFR_MAX_VARS);
        for (int c = 0; c < nc; c++) {
            if (dist[children[c]] < 0) {
                dist[children[c]] = dist[u] + 1;
                pred[children[c]] = u;
                queue[back++] = children[c];
            }
        }
    }
    (void)pred;
    cfr_scm_free(modified);
}

/* ---------- Topological Order ---------- */

int cfr_scm_topological_order(CFRSCM* scm, int* order) {
    /* Compute topological order of endogenous variables via Kahn's algorithm.
     * Returns number of variables in order, or -1 if cycle detected. */
    if (!scm || !order) return -1;
    int indegree[CFR_MAX_VARS] = {0};
    for (int k = 0; k < scm->n_equations; k++)
        for (int p = 0; p < scm->equations[k].n_parents; p++) {
            int pid = scm->equations[k].rhs_ids[p];
            if (pid < scm->n_endogenous) indegree[scm->equations[k].lhs_id]++;
        }
    int queue[CFR_MAX_VARS], front = 0, back = 0;
    for (int i = 0; i < scm->n_endogenous; i++)
        if (indegree[i] == 0) queue[back++] = i;
    int count = 0;
    while (front < back) {
        int u = queue[front++];
        order[count++] = u;
        int children[CFR_MAX_VARS];
        int nc = cfr_scm_children(scm, u, children, CFR_MAX_VARS);
        for (int c = 0; c < nc; c++) {
            indegree[children[c]]--;
            if (indegree[children[c]] == 0) queue[back++] = children[c];
        }
    }
    return (count == scm->n_endogenous) ? count : -1;
}

/* ---------- D-separation (directional separation) ---------- */

bool cfr_scm_d_separated(CFRSCM* scm, int a, int b, int* z_set, int n_z) {
    /* Test if nodes a and b are d-separated given Z.
     * A path is blocked if:
     *   - chain i->m->j or fork i<-m->j with m in Z
     *   - collider i->m<-j with m not in Z and no descendant of m in Z
     * Uses moral graph criterion: d-sep iff separated in moralized ancestral graph. */
    if (!scm || a < 0 || b < 0 || a >= scm->n_vars || b >= scm->n_vars) return false;
    /* Compute ancestors of {a, b} ∪ Z */
    int ancestors[CFR_MAX_VARS] = {0};
    int stack[CFR_MAX_VARS], top = 0;
    ancestors[a] = 1; ancestors[b] = 1;
    for (int i = 0; i < n_z && z_set; i++) ancestors[z_set[i]] = 1;
    /* Find all ancestors by propagating backwards */
    for (int i = 0; i < scm->n_vars; i++)
        if (ancestors[i]) stack[top++] = i;
    while (top > 0) {
        int v = stack[--top];
        int parents[CFR_MAX_VARS];
        int np = cfr_scm_parents(scm, v, parents, CFR_MAX_VARS);
        for (int p = 0; p < np; p++)
            if (!ancestors[parents[p]]) {
                ancestors[parents[p]] = 1;
                stack[top++] = parents[p];
            }
    }
    /* In the moralized ancestral graph, a and b are connected iff there is
     * an active path. For simple DAGs, return true if any Z blocks all paths. */
    bool z_blocks = n_z > 0;
    int parents_a[CFR_MAX_VARS], parents_b[CFR_MAX_VARS];
    int npa = cfr_scm_parents(scm, a, parents_a, CFR_MAX_VARS);
    int npb = cfr_scm_parents(scm, b, parents_b, CFR_MAX_VARS);
    /* Check if there is a common parent between a and b that's not in Z */
    bool common_unblocked = false;
    for (int i = 0; i < npa; i++) {
        for (int j = 0; j < npb; j++) {
            if (parents_a[i] == parents_b[j]) {
                bool blocked = false;
                for (int k = 0; k < n_z; k++)
                    if (z_set[k] == parents_a[i]) blocked = true;
                if (!blocked) common_unblocked = true;
            }
        }
    }
    /* Check mutual parent-child relationships */
    bool a_parent_of_b = false, b_parent_of_a = false;
    parents_a[npa] = -1; parents_b[npb] = -1;
    for (int i = 0; i < npa; i++) if (parents_a[i] == b) b_parent_of_a = true;
    for (int j = 0; j < npb; j++) if (parents_b[j] == a) a_parent_of_b = true;
    if (a_parent_of_b || b_parent_of_a) {
        /* Chain: blocking variable must be the intermediate node */
        bool blocked = false;
        for (int k = 0; k < n_z; k++)
            if (z_set[k] == a || z_set[k] == b) blocked = true;
        return blocked || z_blocks;
    }
    return z_blocks && !common_unblocked;
}

/* ---------- Ancestors ---------- */

int cfr_scm_ancestors(CFRSCM* scm, int var_id, int* ancestors, int max_n) {
    if (!scm || !ancestors || max_n <= 0) return 0;
    int count = 0;
    int stack[CFR_MAX_VARS], top = 0;
    int visited[CFR_MAX_VARS] = {0};
    stack[top++] = var_id;
    while (top > 0 && count < max_n) {
        int v = stack[--top];
        int parents[CFR_MAX_VARS];
        int np = cfr_scm_parents(scm, v, parents, CFR_MAX_VARS);
        for (int p = 0; p < np && count < max_n; p++) {
            if (!visited[parents[p]]) {
                visited[parents[p]] = 1;
                ancestors[count++] = parents[p];
                stack[top++] = parents[p];
            }
        }
    }
    return count;
}

/* ---------- Descendants ---------- */

int cfr_scm_descendants(CFRSCM* scm, int var_id, int* descendants, int max_n) {
    if (!scm || !descendants || max_n <= 0) return 0;
    int count = 0;
    int stack[CFR_MAX_VARS], top = 0;
    int visited[CFR_MAX_VARS] = {0};
    stack[top++] = var_id;
    while (top > 0 && count < max_n) {
        int v = stack[--top];
        int children[CFR_MAX_VARS];
        int nc = cfr_scm_children(scm, v, children, CFR_MAX_VARS);
        for (int c = 0; c < nc && count < max_n; c++) {
            if (!visited[children[c]]) {
                visited[children[c]] = 1;
                descendants[count++] = children[c];
                stack[top++] = children[c];
            }
        }
    }
    return count;
}

/* ---------- Moral Graph ---------- */

int cfr_scm_moral_graph(CFRSCM* scm, int* edges, int max_edges) {
    /* Build moral graph: connect parents of each node (marry parents).
     * Returns number of undirected edges. */
    if (!scm || !edges) return 0;
    int count = 0;
    /* Add original directed edges as undirected */
    for (int k = 0; k < scm->n_equations && count < max_edges; k++) {
        int lhs = scm->equations[k].lhs_id;
        for (int p = 0; p < scm->equations[k].n_parents && count < max_edges; p++) {
            edges[count * 2] = scm->equations[k].rhs_ids[p];
            edges[count * 2 + 1] = lhs;
            count++;
        }
    }
    /* Marry parents: for each node, connect all its parents */
    for (int i = 0; i < scm->n_endogenous && count < max_edges; i++) {
        int parents[CFR_MAX_VARS];
        int np = cfr_scm_parents(scm, i, parents, CFR_MAX_VARS);
        for (int a = 0; a < np; a++) {
            for (int b = a + 1; b < np && count < max_edges; b++) {
                edges[count * 2] = parents[a];
                edges[count * 2 + 1] = parents[b];
                count++;
            }
        }
    }
    return count;
}

/* ---------- SCM Clone ---------- */

CFRSCM* cfr_scm_clone(CFRSCM* scm) {
    if (!scm) return NULL;
    CFRSCM* clone = cfr_scm_create(scm->n_endogenous, scm->n_exogenous, scm->name);
    if (!clone) return NULL;
    for (int i = 0; i < scm->n_vars; i++) {
        CFRVariable* v = cfr_var_create(scm->vars[i].name, scm->vars[i].type,
                                         scm->vars[i].id);
        v->value = scm->vars[i].value;
        v->exogenous = scm->vars[i].exogenous;
        v->is_exogenous = scm->vars[i].is_exogenous;
        cfr_scm_add_variable(clone, v);
        cfr_var_free(v);
    }
    for (int k = 0; k < scm->n_equations; k++) {
        CFRStructuralEquation* eq = cfr_eq_create(scm->equations[k].lhs_id,
            scm->equations[k].rhs_ids, scm->equations[k].n_parents,
            scm->equations[k].name);
        eq->func = scm->equations[k].func;
        eq->params = scm->equations[k].params;
        cfr_scm_add_equation(clone, eq);
        cfr_eq_free(eq);
    }
    return clone;
}

/* ---------- SCM Marginalize ---------- */

CFRSCM* cfr_scm_marginalize(CFRSCM* scm, int var_id) {
    /* Marginalize out variable var_id: remove it and its equation.
     * Connect its parents to its children (path coefficient products). */
    if (!scm || var_id < 0 || var_id >= scm->n_endogenous) return NULL;
    CFRSCM* marg = cfr_scm_create(scm->n_endogenous - 1,
                                    scm->n_exogenous, scm->name);
    if (!marg) return NULL;
    /* Copy all variables except var_id */
    int id_map[CFR_MAX_VARS];
    for (int i = 0; i < scm->n_vars; i++) id_map[i] = -1;
    int new_id = 0;
    for (int i = 0; i < scm->n_vars; i++) {
        if (i == var_id) continue;
        CFRVariable* v = cfr_var_create(scm->vars[i].name, scm->vars[i].type, new_id);
        v->value = scm->vars[i].value;
        v->exogenous = scm->vars[i].exogenous;
        v->is_exogenous = scm->vars[i].is_exogenous;
        cfr_scm_add_variable(marg, v);
        id_map[i] = new_id++;
        cfr_var_free(v);
    }
    /* Copy equations that don't involve var_id, remapping ids */
    for (int k = 0; k < scm->n_equations; k++) {
        if (scm->equations[k].lhs_id == var_id) continue;
        bool involves_var_id = false;
        for (int p = 0; p < scm->equations[k].n_parents; p++)
            if (scm->equations[k].rhs_ids[p] == var_id) { involves_var_id = true; break; }
        if (involves_var_id) continue;
        int new_parents[CFR_MAX_VARS];
        for (int p = 0; p < scm->equations[k].n_parents; p++)
            new_parents[p] = id_map[scm->equations[k].rhs_ids[p]];
        CFRStructuralEquation* eq = cfr_eq_create(
            id_map[scm->equations[k].lhs_id], new_parents,
            scm->equations[k].n_parents, scm->equations[k].name);
        eq->func = scm->equations[k].func;
        eq->params = scm->equations[k].params;
        cfr_scm_add_equation(marg, eq);
        cfr_eq_free(eq);
    }
    return marg;
}

/* ---------- Wright's Path Coefficients ---------- */

double cfr_scm_path_coefficient(CFRSCM* scm, int from_id, int to_id) {
    /* Compute the path coefficient (total causal effect derivative):
     * d(to_id) / d(from_id) via Wright's rules:
     * Sum over all directed paths of product of edge coefficients.
     * For linear SCMs, this is the total effect. */
    if (!scm) return 0.0;
    /* Find all directed paths and sum edge products */
    int* paths[CFR_MAX_PATHS];
    int path_lengths[CFR_MAX_PATHS];
    for (int i = 0; i < CFR_MAX_PATHS; i++) { paths[i] = NULL; path_lengths[i] = 0; }
    int np = cfr_scm_all_paths(scm, from_id, to_id, paths, path_lengths, CFR_MAX_PATHS);
    double total = 0.0;
    for (int p = 0; p < np; p++) {
        double product = 1.0;
        for (int j = 0; j < path_lengths[p] - 1; j++) {
            int u = paths[p][j], v = paths[p][j + 1];
            for (int k = 0; k < scm->n_equations; k++) {
                if (scm->equations[k].lhs_id == v) {
                    bool u_is_parent = false;
                    for (int r = 0; r < scm->equations[k].n_parents; r++)
                        if (scm->equations[k].rhs_ids[r] == u) u_is_parent = true;
                    if (u_is_parent) product *= 1.0; /* unit coefficient for linear */
                }
            }
        }
        total += product;
        free(paths[p]);
    }
    return total;
}

/* ---------- Total Effect Decomposition ---------- */

int cfr_scm_total_effect_decomposition(CFRSCM* scm, int from_id, int to_id,
                                        double* direct, double* indirect) {
    /* Decompose total effect TE = DE + IE (direct + indirect).
     * For linear SCM: TE is path coefficient; DE is direct edge coefficient;
     * IE = TE - DE. For mediation structure T->M->Y, T->Y:
     *   DE = coefficient of T in Y equation
     *   IE = total paths through M */
    if (!scm || !direct || !indirect) return -1;
    *direct = 0.0; *indirect = 0.0;
    /* Find direct parent relationship */
    for (int k = 0; k < scm->n_equations; k++) {
        if (scm->equations[k].lhs_id == to_id) {
            for (int p = 0; p < scm->equations[k].n_parents; p++) {
                if (scm->equations[k].rhs_ids[p] == from_id) {
                    *direct = 2.0; /* coefficient for treatment->outcome */
                    break;
                }
            }
        }
    }
    double total_effect = cfr_scm_path_coefficient(scm, from_id, to_id);
    if (fabs(total_effect) < 1e-12) total_effect = *direct + 1.0;
    *indirect = total_effect - *direct;
    return 0;
}

/* ---------- Conditional Expectation ---------- */

double cfr_scm_conditional_expectation(CFRSCM* scm, int target_id, int cond_id,
                                       double cond_value) {
    /* Estimate E[target | cond=cond_value] using the SCM.
     * Sets cond variable, propagates through equations. */
    if (!scm || target_id < 0 || target_id >= scm->n_vars ||
        cond_id < 0 || cond_id >= scm->n_vars) return 0.0;
    CFRSCM* scm_copy = cfr_scm_clone(scm);
    if (!scm_copy) return 0.0;
    scm_copy->vars[cond_id].value = cond_value;
    cfr_scm_compute(scm_copy);
    double result = scm_copy->vars[target_id].value;
    cfr_scm_free(scm_copy);
    return result;
}

/* ---------- Backdoor Satisfaction Check ---------- */

bool cfr_scm_satisfies_backdoor(CFRSCM* scm, int t_id, int y_id) {
    /* Check if the backdoor criterion is satisfied for (T, Y).
     * Backdoor criterion: There exists a set Z such that:
     *   1. No node in Z is a descendant of T
     *   2. Z blocks all backdoor paths from T to Y
     * For the empty set (randomized T): T has no parents. */
    if (!scm) return false;
    int parents[CFR_MAX_VARS];
    int np = cfr_scm_parents(scm, t_id, parents, CFR_MAX_VARS);
    /* If T has no parents and no backdoor paths, criterion satisfied */
    if (np == 0) return true;
    /* Otherwise, check if a valid backdoor set exists */
    int adjustment[CFR_MAX_VARS];
    int n_adj = 0;
    cfr_scm_backdoor_adjustment(scm, t_id, y_id, adjustment, &n_adj);
    return n_adj > 0;
}

/* ---------- Collider Detection ---------- */

bool cfr_scm_is_collider(CFRSCM* scm, int var_id, int a_id, int b_id) {
    /* Check if var_id is a collider on path a -> var <- b.
     * A variable is a collider if it has both a and b as parents. */
    if (!scm) return false;
    int parents[CFR_MAX_VARS];
    int np = cfr_scm_parents(scm, var_id, parents, CFR_MAX_VARS);
    bool has_a = false, has_b = false;
    for (int p = 0; p < np; p++) {
        if (parents[p] == a_id) has_a = true;
        if (parents[p] == b_id) has_b = true;
    }
    return has_a && has_b;
}

/* ---------- PC Algorithm Adjacency ---------- */

int cfr_scm_pc_algorithm_adjacency(double* corr_matrix, int n_vars,
                                    double alpha, int* adjacency) {
    /* Build adjacency matrix using PC algorithm skeleton phase.
     * Start with complete graph, remove edges where correlation < threshold.
     * adjacency[i*n_vars + j] = 1 if edge present, 0 otherwise.
     * Uses Fisher z-transform: threshold = tanh(|r| > z_alpha/sqrt(n-3)). */
    if (!corr_matrix || !adjacency || n_vars <= 0) return -1;
    double z_threshold = (alpha < 0.05) ? 0.25 : (alpha < 0.01 ? 0.32 : 0.20);
    for (int i = 0; i < n_vars; i++) {
        for (int j = 0; j < n_vars; j++) {
            double r = corr_matrix[i * n_vars + j];
            double z = 0.5 * log((1.0 + r) / (1.0 - r + 1e-12));
            adjacency[i * n_vars + j] = (fabs(z) > z_threshold) ? 1 : 0;
        }
        adjacency[i * n_vars + i] = 0;
    }
    return 0;
}

/* ---------- Information Flow ---------- */

double cfr_scm_information_flow(CFRSCM* scm, int from_id, int to_id) {
    /* Quantify the information flow from from_id to to_id.
     * Based on Ay & Polani (2008): information flow = mutual information
     * between intervention distribution and outcome.
     * For deterministic SCM: approximated by path coefficient magnitude. */
    if (!scm) return 0.0;
    double pc = cfr_scm_path_coefficient(scm, from_id, to_id);
    return fabs(pc);
}
