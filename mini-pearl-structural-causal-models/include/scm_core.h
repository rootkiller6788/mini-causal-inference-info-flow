#ifndef SCM_CORE_H
#define SCM_CORE_H
#include <stdbool.h>
#include <stddef.h>
#define SCM_MAX_VARS 16
#define SCM_MAX_EDGES 64
#define SCM_MAX_PATHS 128
#define SCM_NAME_LEN 32
#define SCM_EPSILON 1e-12
typedef enum { SCM_ENDOGENOUS=0, SCM_EXOGENOUS=1 } SCM_VarType;
typedef struct { int id; char name[SCM_NAME_LEN]; SCM_VarType type; double value; double noise; } SCM_Variable;
typedef struct { int from; int to; double coef; } SCM_Edge;
typedef struct SCM_Model SCM_Model;
typedef double (*SCM_Equation)(SCM_Model* m, int var_id);
struct SCM_Model { int n_vars, n_edges; SCM_Variable vars[SCM_MAX_VARS]; SCM_Edge edges[SCM_MAX_EDGES]; SCM_Equation equations[SCM_MAX_VARS]; bool adjacency[SCM_MAX_VARS][SCM_MAX_VARS]; int topo_order[SCM_MAX_VARS]; bool intervened[SCM_MAX_VARS]; double intervention_value[SCM_MAX_VARS]; };
typedef struct { int n; int nodes[SCM_MAX_VARS]; } SCM_VarSet;
typedef struct { int nodes[SCM_MAX_VARS]; int len; } SCM_Path;
SCM_Model* scm_create(void);
void scm_free(SCM_Model* m);
int  scm_add_variable(SCM_Model* m, const char* name, SCM_VarType type);
void scm_set_equation(SCM_Model* m, int var_id, SCM_Equation eq);
bool scm_is_exogenous(const SCM_Model* m, int var_id);
int  scm_add_edge(SCM_Model* m, int from, int to, double coef);
bool scm_has_edge(const SCM_Model* m, int from, int to);
void scm_build_adjacency(SCM_Model* m);
void scm_topological_sort(SCM_Model* m);
bool scm_is_dag(const SCM_Model* m);
void scm_get_parents(const SCM_Model* m, int var_id, SCM_VarSet* pa);
void scm_get_children(const SCM_Model* m, int var_id, SCM_VarSet* ch);
void scm_get_ancestors(const SCM_Model* m, int var_id, SCM_VarSet* an);
void scm_get_descendants(const SCM_Model* m, int var_id, SCM_VarSet* de);
void scm_markov_blanket(const SCM_Model* m, int var_id, SCM_VarSet* mb);
void scm_sample(SCM_Model* m);
void scm_sample_n(SCM_Model* m, double* data, int n_rows);
double scm_evaluate(SCM_Model* m, int var_id);
void scm_print(const SCM_Model* m);
void scm_varset_init(SCM_VarSet* vs);
void scm_varset_add(SCM_VarSet* vs, int id);
bool scm_varset_contains(const SCM_VarSet* vs, int id);
void scm_varset_union(SCM_VarSet* a, const SCM_VarSet* b, SCM_VarSet* out);
void scm_varset_intersect(SCM_VarSet* a, const SCM_VarSet* b, SCM_VarSet* out);
int  scm_varset_size(const SCM_VarSet* vs);
double scm_linear_eq(SCM_Model* m, int var_id);
#endif
