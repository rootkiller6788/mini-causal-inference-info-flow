#ifndef IFN_DYNAMIC_H
#define IFN_DYNAMIC_H
#include "ifn_core.h"

/* Network information dynamics — how information flows and
 * evolves on network structures over time. */

IFN_NetworkDynamics ifn_network_dynamics_create(const IFN_CausalGraph* g,
    const double* initial_info, int n_steps);
void ifn_dynamics_free(IFN_NetworkDynamics* nd);
void ifn_dynamics_step(IFN_NetworkDynamics* nd, const IFN_CausalGraph* g);
void ifn_dynamics_run(IFN_NetworkDynamics* nd, const IFN_CausalGraph* g, int steps);

/* Network centrality for information flow */
double ifn_pagerank(const IFN_CausalGraph* g, double damping, int max_iter, double tol, double* pr);
double ifn_betweenness_centrality(const IFN_CausalGraph* g, int node);
double ifn_closeness_centrality(const IFN_CausalGraph* g, int node);
double ifn_eigenvector_centrality(const IFN_CausalGraph* g, double* centrality);

/* Network entropy measures */
double ifn_network_entropy(const IFN_CausalGraph* g);
double ifn_degree_entropy(const IFN_CausalGraph* g);
double ifn_flow_entropy(const IFN_NetworkDynamics* nd);

/* Information diffusion models */
void ifn_si_model(const IFN_CausalGraph* g, double beta, int n_steps,
    int* infected_counts);
void ifn_sir_model(const IFN_CausalGraph* g, double beta, double gamma,
    int n_steps, int* susceptible, int* infected, int* recovered);
double ifn_epidemic_threshold(const IFN_CausalGraph* g);

/* Effective information (integrated information theory concept) */
double ifn_effective_information(const IFN_CausalGraph* g,
    const double* state_probs, int mechanism_node);

#endif
