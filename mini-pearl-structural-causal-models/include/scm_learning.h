#ifndef SCM_LEARNING_H
#define SCM_LEARNING_H
#include "scm_core.h"

/* ============================================================================
 * Causal Structure Learning (L5 Algorithms)
 * Reference: Spirtes, Glymour, Scheines (2000) "Causation, Prediction, Search"
 *            Chickering (2002) "Optimal Structure Identification with Greedy Search"
 * ============================================================================ */

typedef struct {
    double correlation;
    double p_value;
    bool independent;
    double z_statistic;
} SCM_CITest;

typedef struct {
    double alpha;
    int max_conditioning;
    bool use_fisher_z;
    int max_iters;
} SCM_PCConfig;

typedef struct {
    bool adjacency[SCM_MAX_VARS][SCM_MAX_VARS];
    bool directed[SCM_MAX_VARS][SCM_MAX_VARS];
    int  orientation_rules_applied;
    double bic_score;
    int  n_edges;
} SCM_LearnResult;

SCM_CITest scm_partial_correlation_test(const double* data, int n_rows, int n_cols,
                                         int x, int y, const SCM_VarSet* Z, double alpha);
double scm_fisher_z(double r, int n, int k);
SCM_CITest scm_mutual_information_test(const double* data, int n_rows, int n_cols,
                                        int x, int y, const SCM_VarSet* Z, double alpha);
SCM_PCConfig scm_pc_config_default(void);
SCM_LearnResult* scm_pc_algorithm(const double* data, int n_rows, int n_cols,
                                   const SCM_PCConfig* config);
void scm_skeleton_discovery(const double* data, int n_rows, int n_cols,
                             SCM_LearnResult* result, const SCM_PCConfig* config);
void scm_orient_edges(SCM_LearnResult* result, const SCM_VarSet* separating_sets,
                       int n_vars);
int scm_orient_rule1(SCM_LearnResult* r, const SCM_VarSet* sepsets, int n_vars);
int scm_orient_rule2(SCM_LearnResult* r, int n_vars);
int scm_orient_rule3(SCM_LearnResult* r, int n_vars);
SCM_LearnResult* scm_ges_algorithm(const double* data, int n_rows, int n_cols,
                                    const SCM_PCConfig* config);
SCM_LearnResult* scm_hill_climbing_search(const double* data, int n_rows, int n_cols,
                                           const SCM_PCConfig* config);
SCM_LearnResult* scm_tabu_search_causal(const double* data, int n_rows, int n_cols,
                                         const SCM_PCConfig* config, int tabu_tenure);
double scm_bic_score(const double* data, int n_rows, int n_cols,
                      const SCM_LearnResult* graph);
double scm_aic_score(const double* data, int n_rows, int n_cols,
                      const SCM_LearnResult* graph);
SCM_LearnResult* scm_fci_algorithm(const double* data, int n_rows, int n_cols,
                                    const SCM_PCConfig* config);
bool scm_markov_equivalent(const SCM_LearnResult* a, const SCM_LearnResult* b, int n_vars);
void scm_bootstrap_stability(const double* data, int n_rows, int n_cols,
                              const SCM_PCConfig* config, int n_bootstrap,
                              double stability[SCM_MAX_VARS][SCM_MAX_VARS]);
void scm_learn_result_init(SCM_LearnResult* r);
void scm_learn_result_free(SCM_LearnResult* r);
void scm_learn_result_print(const SCM_LearnResult* r, int n_vars);

#endif
