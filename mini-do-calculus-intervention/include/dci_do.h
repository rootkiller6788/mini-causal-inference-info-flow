#ifndef DCI_DO_H
#define DCI_DO_H

#include "dci_core.h"
#include "dci_graph.h"

/* ---- do-Calculus (Pearl, 1995) ---- */

/* Intervention structure */
typedef struct {
    int vars[DCI_MAX_VARS];     /* Variables to intervene on */
    double values[DCI_MAX_VARS]; /* Values to set */
    int n_vars;
    char name[128];
} DCI_Intervention;

/* Truncated SCM (after intervention) */
typedef struct {
    DCI_SCM* original;
    DCI_Intervention intervention;
    DCI_Graph* mutilated_graph;
    DCI_SCM* post_intervention_scm;
} DCI_TruncatedSCM;

/* do-Calculus Rule Set */
typedef enum {
    DCI_RULE1 = 1,  /* Insertion/deletion of observations */
    DCI_RULE2 = 2,  /* Action/observation exchange */
    DCI_RULE3 = 3   /* Insertion/deletion of actions */
} DCI_RuleType;

/* do-Calculus derivation step */
typedef struct {
    DCI_RuleType rule;
    char expression_before[256];
    char expression_after[256];
    char justification[256];
} DCI_DerivationStep;

/* do-Calculus derivation trace */
typedef struct {
    DCI_DerivationStep steps[32];
    int n_steps;
    bool success;
} DCI_Derivation;

/* ---- Intervention Operations ---- */
DCI_Intervention dci_intervention_create(const int* vars,
    const double* values, int n_vars);
DCI_TruncatedSCM* dci_truncated_scm_create(DCI_SCM* scm,
    DCI_Intervention* intervention);
void dci_truncated_scm_free(DCI_TruncatedSCM* tscm);
double dci_truncated_scm_probability(DCI_TruncatedSCM* tscm,
    int var_id, double value);

/* do-Calculus Rules */
DCI_Derivation dci_apply_rule1(DCI_SCM* scm, DCI_Graph* g,
    int y, const int* x, int n_x, const int* z, int n_z, const int* w, int n_w);
DCI_Derivation dci_apply_rule2(DCI_SCM* scm, DCI_Graph* g,
    int y, const int* x, int n_x, const int* z, int n_z, const int* w, int n_w);
DCI_Derivation dci_apply_rule3(DCI_SCM* scm, DCI_Graph* g,
    int y, const int* x, int n_x, const int* z, int n_z, const int* w, int n_w);

/* Full do-calculus derivation */
DCI_Derivation dci_do_calculus(DCI_SCM* scm,
    const char* target_expression, bool* identifiable);

/* ---- Truncated Factorization ---- */
double dci_truncated_factorization(DCI_SCM* scm,
    const int* intervention_vars, const double* intervention_vals,
    int n_intervention, const double* observed_vals);

/* ---- Causal Effect via do-Calculus ---- */
double dci_causal_effect_do(DCI_SCM* scm, int cause, int effect,
    double cause_value);

/* ---- Print derivation ---- */
void dci_derivation_print(const DCI_Derivation* der);

bool dci_identify_effect(DCI_SCM* scm, int x, int y, char* derivation, int max_len);
void dci_truncated_factorization_graph(DCI_SCM* scm, DCI_Graph* g, const int* x, int n_x, double* probs, int n_probs);
bool dci_identifiability_check(DCI_SCM* scm, int x, int y);
int dci_direct_causes(DCI_Graph* g, int node, int* direct_causes);
double dci_g_computation(DCI_SCM* scm, int treatment, int outcome, double treat_val, const int* time_vars, int n_time_vars, int n_samples);
double dci_sequential_backdoor(DCI_SCM* scm, const int* treatments, const double* treat_vals, int n_timepoints, int outcome, int n_samples);

#endif
