#include "scm_core.h"
#include "scm_graph.h"
#include "scm_intervention.h"
#include "scm_identification.h"
#include "scm_counterfactual.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define EPS 1e-9
int main(void) {
    SCM_Model* m = scm_create(); assert(m);
    int z=scm_add_variable(m,"Z",SCM_EXOGENOUS);
    int x=scm_add_variable(m,"X",SCM_ENDOGENOUS);
    int y=scm_add_variable(m,"Y",SCM_ENDOGENOUS);
    assert(z==0&&x==1&&y==2);
    scm_add_edge(m,z,x,0.8); scm_add_edge(m,z,y,0.6); scm_add_edge(m,x,y,0.5);
    scm_set_equation(m,x,scm_linear_eq); scm_set_equation(m,y,scm_linear_eq);
    scm_build_adjacency(m); scm_topological_sort(m);
    assert(scm_is_dag(m)); assert(scm_has_edge(m,z,x));
    assert(!scm_is_exogenous(m,x)); assert(scm_is_exogenous(m,z));
    SCM_VarSet pa; scm_get_parents(m,y,&pa); assert(pa.n>=2);
    SCM_VarSet ch; scm_get_children(m,z,&ch); assert(ch.n>=2);
    SCM_VarSet an; scm_get_ancestors(m,y,&an); assert(scm_varset_contains(&an,z));
    SCM_VarSet de; scm_get_descendants(m,z,&de); assert(scm_varset_contains(&de,y));
    SCM_VarSet mb; scm_markov_blanket(m,x,&mb); assert(mb.n>0);
    m->vars[z].noise=0.5; scm_sample(m);
    assert(!isnan(m->vars[x].value)); assert(!isnan(m->vars[y].value));
    SCM_VarSet Zset; scm_varset_init(&Zset); scm_varset_add(&Zset,z);
    assert(scm_d_separated(m,x,y,&Zset));
    SCM_Collider colliders[16]; int nc; scm_find_colliders(m,colliders,&nc);
    bool moral[SCM_MAX_VARS][SCM_MAX_VARS]; scm_moral_graph(m,moral);
    SCM_Path p; assert(scm_find_directed_path(m,z,y,&p));
    SCM_Path paths[8]; int np=scm_all_directed_paths(m,z,y,paths,8); assert(np>0);
    SCM_Model* m_do=scm_do_intervention(m,x,1.0); assert(m_do);
    assert(m_do->intervened[x]); scm_free(m_do);
    SCM_CausalEffect* ce=scm_causal_effect_compute(m,x,y); assert(ce);
    scm_causal_effect_print(ce); scm_causal_effect_free(ce);
    double ace=scm_ace_linear(m,x,y); assert(fabs(ace-0.5)<0.01);
    SCM_BackDoorResult* br=scm_back_door(m,x,y); assert(br); assert(br->found);
    scm_back_door_print(br);
    SCM_FrontDoorResult* fr=scm_front_door(m,z,y,x); scm_front_door_print(fr);
    scm_front_door_free(fr);
    assert(scm_is_identifiable(m,x,y));
    SCM_VarSet sets[8]; int ns; scm_adjustment_sets(m,x,y,sets,&ns,8); assert(ns>0);
    int instr; scm_find_instrument(m,x,y,&instr);
    SCM_VarSet empty; scm_varset_init(&empty); bool app;
    scm_do_calculus_rule1(m,x,y,&empty,&empty,&app);
    scm_do_calculus_rule2(m,x,y,&empty,&app);
    scm_do_calculus_rule3(m,x,y,&empty,&empty,&app);
    int ev[]={x}; double evv[]={1.0};
    int iv[]={x}; double ivv[]={0.0};
    SCM_Counterfactual* cf=scm_counterfactual_compute(m,y,ev,evv,1,iv,ivv,1);
    assert(cf); scm_counterfactual_print(cf); scm_counterfactual_free(cf);
    SCM_ProbNecessitySuff* pns=scm_prob_necessity_sufficiency(m,x,y,1.0,1.0);
    assert(pns); scm_pns_print(pns); scm_pns_free(pns);
    double nde,nie; scm_mediation_analysis(m,x,z,y,&nde,&nie);
    assert(scm_monotonicity_check(m,x,y));
    bool has_uc; scm_unobserved_confounder_test(m,x,y,&has_uc);
    scm_back_door_free(br); scm_free(m);
    printf("All tests passed.\n");
    return 0;
}
