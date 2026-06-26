#ifndef DCI_COUNTERFACTUAL_H
#define DCI_COUNTERFACTUAL_H

#include "dci_core.h"
#include "dci_graph.h"
#include "dci_do.h"

/* ---- Counterfactual Reasoning (Pearl, 2000) ---- */

/* Counterfactual query */
typedef struct {
    int antecedent_var;
    double antecedent_value;
    int consequent_var;
    double consequent_value;
    char description[256];
} DCI_Counterfactual;

/* Counterfactual probability:
 * P(Y_{X=x} = y | E=e) where E is evidence
 */
typedef struct {
    double probability;
    int n_samples;
    bool is_identifiable;
} DCI_CounterfactualResult;

/* Twin network method */
typedef struct {
    DCI_SCM* factual;
    DCI_SCM* counterfactual;
    DCI_Graph* twin_graph;
} DCI_TwinNetwork;

/* Counterfactual computation */
DCI_CounterfactualResult dci_counterfactual_probability(
    DCI_SCM* scm, DCI_Counterfactual* query,
    const double* evidence, const int* evidence_vars,
    int n_evidence, int n_samples);

/* Twin network construction */
DCI_TwinNetwork* dci_twin_network_create(DCI_SCM* scm,
    int intervention_var, double intervention_value);
void dci_twin_network_free(DCI_TwinNetwork* tn);
double dci_twin_network_query(DCI_TwinNetwork* tn,
    int query_var, double query_value);

/* Probability of necessity and sufficiency
 * PN = P(Y_{X=0}=0 | X=1, Y=1) — probability that X was necessary
 * PS = P(Y_{X=1}=1 | X=0, Y=0) — probability that X was sufficient
 * PNS = P(Y_{X=1}=1, Y_{X=0}=0) — probability of necessity AND sufficiency
 */
double dci_probability_of_necessity(DCI_SCM* scm, int x, int y,
    int n_samples);
double dci_probability_of_sufficiency(DCI_SCM* scm, int x, int y,
    int n_samples);
double dci_probability_of_necessity_sufficiency(DCI_SCM* scm,
    int x, int y, int n_samples);

/* Individual-level counterfactual (specific factual case) */
double dci_individual_counterfactual(DCI_SCM* scm,
    const double* observed_values, const int* observed_vars,
    int n_observed, int intervention_var, double intervention_value,
    int query_var);

/* Mediation via counterfactuals */
double dci_causal_mediation(DCI_SCM* scm, int cause, int mediator,
    int effect, int n_samples);
bool dci_counterfactual_fairness(DCI_SCM* scm, int protected_attr, int outcome, double threshold, int n_samples);
double dci_path_specific_effect(DCI_SCM* scm, int cause, int effect, const int* path_nodes, int path_len, int n_samples);
void dci_counterfactual_explanation(DCI_SCM* scm, const double* observed, const int* obs_vars, int n_obs, int target_var, double* importance, int max_vars);
double dci_attributable_fraction(DCI_SCM* scm, int exposure, int outcome, int n_samples);

#endif
