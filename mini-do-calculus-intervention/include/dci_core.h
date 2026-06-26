#ifndef DCI_CORE_H
#define DCI_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * dci_core.h — Core Types for Do-Calculus & Causal Inference
 *
 * Structural Causal Model (SCM): M = <U, V, F, P(u)>
 *   U: exogenous (unobserved) variables
 *   V: endogenous (observed) variables
 *   F: structural equations v_i = f_i(pa_i, u_i)
 *   P(u): distribution over exogenous variables
 *
 * Intervention do(X=x): Replace f_X with X = x (graph surgery)
 *
 * References:
 *   Pearl (2009) Causality: Models, Reasoning, and Inference
 *   Pearl (1995) Causal Diagrams for Empirical Research
 *   Pearl & Mackenzie (2018) The Book of Why
 *   Spirtes, Glymour, Scheines (2000) Causation, Prediction, and Search
 * ============================================================== */

#define DCI_MAX_VARS      32
#define DCI_MAX_PARENTS   16
#define DCI_MAX_PATHS     128
#define DCI_MAX_LEN       1024
#define DCI_MAX_PTS       5000

/* ---- Variable Types ---- */
typedef enum { DCI_CONTINUOUS = 0, DCI_BINARY = 1, DCI_DISCRETE = 2 } DCI_VarType;

/* ---- Structural Equation ---- */
typedef double (*DCI_StructuralEq)(const double* parents, const double* u);
typedef void (*DCI_StructuralEqVec)(const double* parents, const double* u, double* out);

/* ---- Variable in SCM ---- */
typedef struct {
    char name[64];
    int id;
    DCI_VarType type;
    double value;
    int parent_ids[DCI_MAX_PARENTS];
    int n_parents;
    int child_ids[DCI_MAX_VARS];
    int n_children;
    DCI_StructuralEq equation;
    bool is_exogenous;
} DCI_Variable;

/* ---- Structural Causal Model ---- */
typedef struct {
    DCI_Variable vars[DCI_MAX_VARS];
    int n_vars;
    char name[128];
} DCI_SCM;

/* ---- Graph Adjacency ---- */
typedef struct {
    double adjacency[DCI_MAX_VARS][DCI_MAX_VARS];
    int n_nodes;
    char names[DCI_MAX_VARS][64];
} DCI_Graph;

/* ---- SCM Operations ---- */
DCI_SCM* dci_scm_create(const char* name);
void dci_scm_free(DCI_SCM* scm);
int dci_scm_add_variable(DCI_SCM* scm, const char* name, DCI_VarType type,
    int* parent_ids, int n_parents, DCI_StructuralEq eq, bool exogenous);
void dci_scm_add_edge(DCI_SCM* scm, int from_id, int to_id);
bool dci_scm_is_dag(DCI_SCM* scm);
double dci_scm_evaluate(DCI_SCM* scm, const double* exog_vals);
void dci_scm_print(DCI_SCM* scm);

/* ---- Variable Operations ---- */
DCI_Variable* dci_var_by_name(DCI_SCM* scm, const char* name);
DCI_Variable* dci_var_by_id(DCI_SCM* scm, int id);
int dci_var_count_parents(DCI_SCM* scm, int var_id);
int dci_var_count_children(DCI_SCM* scm, int var_id);

/* ---- Graph Operations ---- */
DCI_Graph* dci_graph_create(int n_nodes);
void dci_graph_free(DCI_Graph* g);
void dci_graph_add_edge(DCI_Graph* g, int from, int to, double weight);
void dci_graph_from_scm(DCI_SCM* scm, DCI_Graph* g);
void dci_graph_print(DCI_Graph* g);

double dci_scm_query(DCI_SCM* scm, int query_var, const int* interventions, const double* intervention_values, int n_interventions);
DCI_SCM* dci_scm_copy(DCI_SCM* scm);
DCI_SCM* dci_scm_simpsons_paradox(void);
DCI_SCM* dci_scm_iv(int n_samples);
DCI_SCM* dci_scm_mediation(void);
DCI_SCM* dci_scm_m_bias(void);
DCI_SCM* dci_scm_frontdoor(void);
int dci_scm_var_count(const DCI_SCM* scm);
double dci_scm_get_value(const DCI_SCM* scm, int var_id);

#endif
