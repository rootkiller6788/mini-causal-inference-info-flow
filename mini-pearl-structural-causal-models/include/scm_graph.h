#ifndef SCM_GRAPH_H
#define SCM_GRAPH_H
#include "scm_core.h"
typedef struct { int parent1, parent2, child; } SCM_Collider;
bool scm_is_collider(const SCM_Model* m, int a, int b, int c);
void scm_find_colliders(const SCM_Model* m, SCM_Collider* colliders, int* n);
bool scm_d_separated(const SCM_Model* m, int x, int y, const SCM_VarSet* Z);
bool scm_d_separated_set(const SCM_Model* m, const SCM_VarSet* X, const SCM_VarSet* Y, const SCM_VarSet* Z);
void scm_moral_graph(const SCM_Model* m, bool moral[SCM_MAX_VARS][SCM_MAX_VARS]);
bool scm_find_directed_path(const SCM_Model* m, int from, int to, SCM_Path* p);
int  scm_all_directed_paths(const SCM_Model* m, int from, int to, SCM_Path* paths, int max_paths);
void scm_path_print(const SCM_Path* p);
bool scm_is_blocked(const SCM_Path* p, const SCM_VarSet* Z, const SCM_Model* m);
bool scm_path_has_collider(const SCM_Path* p, const SCM_Model* m);
bool scm_is_ancestor(const SCM_Model* m, int a, int b);
#endif
