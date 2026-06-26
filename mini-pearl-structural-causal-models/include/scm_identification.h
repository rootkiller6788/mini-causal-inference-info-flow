#ifndef SCM_IDENTIFICATION_H
#define SCM_IDENTIFICATION_H
#include "scm_core.h"
typedef struct { bool found; SCM_VarSet Z; bool sufficient; } SCM_BackDoorResult;
typedef struct { bool found; SCM_VarSet Z; SCM_VarSet M; bool sufficient; } SCM_FrontDoorResult;
SCM_BackDoorResult* scm_back_door(const SCM_Model* m, int treatment, int outcome);
void scm_back_door_free(SCM_BackDoorResult* br);
void scm_back_door_print(const SCM_BackDoorResult* br);
SCM_FrontDoorResult* scm_front_door(const SCM_Model* m, int treatment, int outcome, int mediator);
void scm_front_door_free(SCM_FrontDoorResult* fr);
void scm_front_door_print(const SCM_FrontDoorResult* fr);
bool scm_is_back_door_admissible(const SCM_Model* m, int treatment, int outcome, const SCM_VarSet* Z);
bool scm_is_front_door_admissible(const SCM_Model* m, int treatment, int outcome, int mediator);
bool scm_is_identifiable(const SCM_Model* m, int treatment, int outcome);
void scm_adjustment_sets(const SCM_Model* m, int treatment, int outcome, SCM_VarSet* sets, int* n_sets, int max_sets);
bool scm_is_minimal_adjustment(const SCM_Model* m, int treatment, int outcome, const SCM_VarSet* Z);
void scm_find_instrument(const SCM_Model* m, int treatment, int outcome, int* instrument);
bool scm_is_instrument(const SCM_Model* m, int z, int x, int y);
void scm_unobserved_confounder_test(const SCM_Model* m, int x, int y, bool* has_confounder);
#endif
