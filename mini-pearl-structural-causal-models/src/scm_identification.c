#include "scm_identification.h"
#include "scm_graph.h"
#include <stdlib.h>
#include <stdio.h>

SCM_BackDoorResult* scm_back_door(const SCM_Model* m, int treatment, int outcome) {
    SCM_BackDoorResult* br = calloc(1, sizeof(SCM_BackDoorResult));
    if (!br) return NULL;
    scm_varset_init(&br->Z);
    for (int i = 0; i < m->n_vars; i++) {
        if (i == treatment || i == outcome) continue;
        if (scm_is_ancestor(m, i, treatment) && scm_is_ancestor(m, i, outcome))
            scm_varset_add(&br->Z, i);
    }
    br->found = (br->Z.n > 0);
    br->sufficient = scm_is_back_door_admissible(m, treatment, outcome, &br->Z);
    return br;
}
void scm_back_door_free(SCM_BackDoorResult* br) { free(br); }
void scm_back_door_print(const SCM_BackDoorResult* br) {
    printf("BackDoor: found=%s sufficient=%s Z={", br->found?"YES":"NO", br->sufficient?"YES":"NO");
    for (int i = 0; i < br->Z.n; i++) printf("%d ", br->Z.nodes[i]);
    printf("}\n");
}
bool scm_is_back_door_admissible(const SCM_Model* m, int treatment, int outcome,
                                   const SCM_VarSet* Z) {
    for (int i = 0; i < Z->n; i++)
        if (scm_is_ancestor(m, treatment, Z->nodes[i])) return false;
    SCM_VarSet X, Y;
    scm_varset_init(&X); scm_varset_add(&X, treatment);
    scm_varset_init(&Y); scm_varset_add(&Y, outcome);
    return scm_d_separated_set(m, &X, &Y, Z);
}
SCM_FrontDoorResult* scm_front_door(const SCM_Model* m, int treatment,
                                      int outcome, int mediator) {
    SCM_FrontDoorResult* fr = calloc(1, sizeof(SCM_FrontDoorResult));
    if (!fr) return NULL;
    fr->found = scm_is_front_door_admissible(m, treatment, outcome, mediator);
    scm_varset_init(&fr->M); scm_varset_add(&fr->M, mediator);
    SCM_VarSet pa; scm_get_parents(m, treatment, &pa);
    scm_varset_union(&fr->Z, &pa, &fr->Z);
    fr->sufficient = fr->found;
    return fr;
}
void scm_front_door_free(SCM_FrontDoorResult* fr) { free(fr); }
void scm_front_door_print(const SCM_FrontDoorResult* fr) {
    printf("FrontDoor: found=%s M={", fr->found?"YES":"NO");
    for (int i = 0; i < fr->M.n; i++) printf("%d ", fr->M.nodes[i]);
    printf("} Z={");
    for (int i = 0; i < fr->Z.n; i++) printf("%d ", fr->Z.nodes[i]);
    printf("}\n");
}
bool scm_is_front_door_admissible(const SCM_Model* m, int treatment,
                                    int outcome, int mediator) {
    SCM_Path p1, p2;
    if (!scm_find_directed_path(m, treatment, mediator, &p1)) return false;
    if (!scm_find_directed_path(m, mediator, outcome, &p2)) return false;
    SCM_VarSet Z; scm_varset_init(&Z); scm_varset_add(&Z, treatment);
    SCM_VarSet O; scm_varset_init(&O); scm_varset_add(&O, outcome);
    SCM_VarSet mediator_set; scm_varset_init(&mediator_set);
    for (int i = 0; i < p1.len; i++) scm_varset_add(&mediator_set, p1.nodes[i]);
    return scm_d_separated_set(m, &Z, &O, &mediator_set);
}
bool scm_is_identifiable(const SCM_Model* m, int treatment, int outcome) {
    SCM_BackDoorResult* br = scm_back_door(m, treatment, outcome);
    bool id = br->found && br->sufficient;
    scm_back_door_free(br);
    return id;
}
void scm_adjustment_sets(const SCM_Model* m, int treatment, int outcome,
                           SCM_VarSet* sets, int* n_sets, int max_sets) {
    *n_sets = 0;
    SCM_BackDoorResult* br = scm_back_door(m, treatment, outcome);
    if (br->found && br->sufficient && *n_sets < max_sets)
        sets[(*n_sets)++] = br->Z;
    scm_back_door_free(br);
}
bool scm_is_minimal_adjustment(const SCM_Model* m, int treatment, int outcome,
                                 const SCM_VarSet* Z) {
    if (!scm_is_back_door_admissible(m, treatment, outcome, Z)) return false;
    for (int i = 0; i < Z->n; i++) {
        SCM_VarSet Zminus; scm_varset_init(&Zminus);
        for (int j = 0; j < Z->n; j++)
            if (j != i) scm_varset_add(&Zminus, Z->nodes[j]);
        if (scm_is_back_door_admissible(m, treatment, outcome, &Zminus)) return false;
    }
    return true;
}
void scm_find_instrument(const SCM_Model* m, int treatment, int outcome, int* instrument) {
    *instrument = -1;
    for (int z = 0; z < m->n_vars; z++)
        if (scm_is_instrument(m, z, treatment, outcome)) { *instrument = z; return; }
}
bool scm_is_instrument(const SCM_Model* m, int z, int x, int y) {
    SCM_VarSet empty; scm_varset_init(&empty);
    if (!scm_d_separated(m, z, y, &empty)) return false;
    if (!scm_has_edge(m, z, x)) return false;
    SCM_VarSet U; scm_varset_init(&U); scm_varset_add(&U, z);
    return scm_d_separated(m, x, y, &U);
}
void scm_unobserved_confounder_test(const SCM_Model* m, int x, int y,
                                      bool* has_confounder) {
    SCM_VarSet empty; scm_varset_init(&empty);
    *has_confounder = !scm_d_separated(m, x, y, &empty) && !scm_has_edge(m, x, y) && !scm_has_edge(m, y, x);
}
