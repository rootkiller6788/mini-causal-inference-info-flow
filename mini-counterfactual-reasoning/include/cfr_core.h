#ifndef CFR_CORE_H
#define CFR_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================
 * cfr_core.h -- Counterfactual Reasoning Core
 *
 * Implements Structural Causal Models (SCM) per Pearl (2009).
 *
 * An SCM is a triple (U, V, F) where:
 *   U = exogenous (background) variables
 *   V = endogenous (observed) variables
 *   F = structural equations: v_i = f_i(pa_i, u_i)
 *
 * Counterfactuals answer "what if" questions:
 *   "What would Y have been if X had been x, given we observed X=x'?"
 *
 * The 3-step counterfactual procedure:
 *   1. Abduction: infer P(U | evidence) from observed data
 *   2. Action:   intervene do(X=x), replacing f_X with X=x
 *   3. Prediction: compute Y_x using modified model and updated U
 *
 * Refs: Pearl (2009) "Causality", Pearl, Glymour & Jewell (2016)
 * ============================================================ */

#define CFR_EPSILON   1e-12
#define CFR_MAX_VARS  20
#define CFR_MAX_EQNS  50
#define CFR_INF       1e99

/* --- Variable Types --- */
typedef enum {
    CFR_VAR_CONTINUOUS = 0,
    CFR_VAR_BINARY     = 1,
    CFR_VAR_DISCRETE   = 2
} CFRVarType;

typedef struct {
    char* name;         /* variable name */
    CFRVarType type;    /* variable type */
    int id;             /* unique identifier */
    double value;       /* current value */
    double exogenous;   /* exogenous noise U */
    bool is_exogenous;  /* if true, this is a U variable */
} CFRVariable;

/* --- Structural Equation --- */
typedef struct {
    char* name;
    int lhs_id;         /* left-hand side variable id */
    int* rhs_ids;       /* right-hand side variable ids (parents) */
    int n_parents;      /* number of parents */
    double (*func)(double* parent_values, double u, void* params);
    void* params;       /* equation parameters */
} CFRStructuralEquation;

/* --- Structural Causal Model --- */
typedef struct {
    CFRVariable* vars;                /* all variables (U + V) */
    int n_vars;                       /* total variables */
    int n_endogenous;                 /* number of V variables */
    int n_exogenous;                  /* number of U variables */
    CFRStructuralEquation* equations; /* structural equations */
    int n_equations;                  /* number of equations */
    char* name;
} CFRSCM;

/* --- Intervention --- */
typedef enum { CFR_INT_HARD = 0, CFR_INT_SOFT = 1 } CFRInterventionType;

typedef struct {
    int var_id;                   /* intervened variable */
    double value;                 /* intervention value */
    CFRInterventionType type;     /* hard (do) or soft */
    double (*soft_func)(double* parent_values, double u, void* params);
    void* soft_params;
} CFRIntervention;

/* --- Evidence --- */
typedef struct {
    int* var_ids;       /* observed variable ids */
    double* values;     /* observed values */
    int n_observations; /* number of observations */
} CFREvidence;

/* --- Counterfactual Query --- */
typedef struct {
    int target_id;              /* variable we want to know about */
    CFRIntervention* intervention; /* hypothetical intervention */
    CFREvidence* evidence;      /* what we actually observed */
    double result;              /* computed counterfactual value */
    bool computed;
} CFRQuery;

/* ================================================================
 * Core API
 * ================================================================ */

/* --- Variable Lifecycle --- */
CFRVariable* cfr_var_create(char* name, CFRVarType type, int id);
void         cfr_var_free(CFRVariable* v);
void         cfr_var_set_value(CFRVariable* v, double val);
void         cfr_var_set_exogenous(CFRVariable* v, double u);

/* --- Equation Lifecycle --- */
CFRStructuralEquation* cfr_eq_create(int lhs_id, int* rhs_ids,
                                      int n_parents, char* name);
void cfr_eq_free(CFRStructuralEquation* eq);
void cfr_eq_set_func(CFRStructuralEquation* eq,
    double (*func)(double*, double, void*), void* params);
double cfr_eq_evaluate(CFRStructuralEquation* eq,
                        double* parent_values, double exogenous);

/* --- SCM Lifecycle --- */
CFRSCM* cfr_scm_create(int n_endogenous, int n_exogenous, char* name);
void    cfr_scm_free(CFRSCM* scm);
int     cfr_scm_add_variable(CFRSCM* scm, CFRVariable* v);
int     cfr_scm_add_equation(CFRSCM* scm, CFRStructuralEquation* eq);
double* cfr_scm_get_values(CFRSCM* scm);
void    cfr_scm_set_exogenous(CFRSCM* scm, double* u_values);
void    cfr_scm_compute(CFRSCM* scm);
void    cfr_scm_compute_topological(CFRSCM* scm);

/* --- Intervention --- */
CFRIntervention* cfr_int_create_hard(int var_id, double value);
CFRIntervention* cfr_int_create_soft(int var_id,
    double (*soft_func)(double*, double, void*), void* params);
void cfr_int_free(CFRIntervention* intv);
CFRSCM* cfr_scm_apply_intervention(CFRSCM* scm, CFRIntervention* intv);

/* --- Evidence --- */
CFREvidence* cfr_ev_create(int n_obs);
void         cfr_ev_free(CFREvidence* ev);
void         cfr_ev_add_observation(CFREvidence* ev, int var_id,
                                     double value, int idx);

/* --- Query --- */
CFRQuery* cfr_query_create(int target_id, CFRIntervention* intv,
                            CFREvidence* ev);
void      cfr_query_free(CFRQuery* q);
int       cfr_query_compute(CFRQuery* q, CFRSCM* scm);

/* --- Factory: Standard SCMs --- */
CFRSCM* cfr_create_linear_scm(int n_vars);
CFRSCM* cfr_create_simple_treatment_scm(void);
CFRSCM* cfr_create_mediation_scm(void);
CFRSCM* cfr_create_frontdoor_scm(void);

/* --- Utility --- */
double cfr_clamp(double x, double lo, double hi);
double cfr_gaussian_noise(double mean, double std);
void   cfr_scm_print(CFRSCM* scm);
int    cfr_find_var_index(CFRSCM* scm, char* name);

/* --- DAG Analysis --- */
int  cfr_scm_topological_order(CFRSCM* scm, int* order);
bool cfr_scm_d_separated(CFRSCM* scm, int a, int b, int* z_set, int n_z);
int  cfr_scm_ancestors(CFRSCM* scm, int var_id, int* ancestors, int max_n);
int  cfr_scm_descendants(CFRSCM* scm, int var_id, int* descendants, int max_n);
int  cfr_scm_moral_graph(CFRSCM* scm, int* edges, int max_edges);
bool cfr_scm_is_dag(CFRSCM* scm);

/* --- SCM Transformations --- */
CFRSCM* cfr_scm_clone(CFRSCM* scm);
CFRSCM* cfr_scm_marginalize(CFRSCM* scm, int var_id);

/* --- Wright's Path Analysis --- */
double cfr_scm_path_coefficient(CFRSCM* scm, int from_id, int to_id);
int    cfr_scm_total_effect_decomposition(CFRSCM* scm, int from_id, int to_id,
                                          double* direct, double* indirect);

/* --- Conditional Expectation --- */
double cfr_scm_conditional_expectation(CFRSCM* scm, int target_id, int cond_id,
                                       double cond_value);

/* --- Backdoor Criterion --- */
bool cfr_scm_satisfies_backdoor(CFRSCM* scm, int t_id, int y_id);
void cfr_scm_backdoor_adjustment(CFRSCM* scm, int treatment_id,
    int outcome_id, int* adjustment_set, int* n_adjust);

/* --- Collider / Selection Bias --- */
bool cfr_scm_is_collider(CFRSCM* scm, int var_id, int a_id, int b_id);

/* --- DAG Causal Discovery --- */
int  cfr_scm_pc_algorithm_adjacency(double* corr_matrix, int n_vars,
                                     double alpha, int* adjacency);

double cfr_scm_information_flow(CFRSCM* scm, int from_id, int to_id);

/* --- Counterfactual Simulation & Paths --- */
int  cfr_scm_all_paths(CFRSCM* scm, int from_id, int to_id,
                       int** paths, int* path_lengths, int max_paths);
bool cfr_scm_frontdoor_criterion(CFRSCM* scm, int treatment_id,
                                  int outcome_id, int* mediator_id);
void cfr_scm_counterfactual_simulation(CFRSCM* scm, int n_samples,
    double* exogenous_samples, double* outcomes);
void cfr_scm_compute_intervention_path(CFRSCM* scm,
    CFRIntervention* intv, int target_id);

int  cfr_scm_parents(CFRSCM* scm, int var_id, int* parent_ids, int max_parents);
int  cfr_scm_children(CFRSCM* scm, int var_id, int* child_ids, int max_children);

#endif /* CFR_CORE_H */
