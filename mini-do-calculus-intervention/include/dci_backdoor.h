#ifndef DCI_BACKDOOR_H
#define DCI_BACKDOOR_H

#include "dci_core.h"
#include "dci_graph.h"

/* ---- Back-Door & Front-Door Criteria (Pearl, 1993) ---- */

/* Back-door admissible set */
typedef struct {
    int variables[DCI_MAX_VARS];
    int n_vars;
    bool is_valid;
} DCI_BackdoorSet;

/* Front-door set */
typedef struct {
    int mediators[DCI_MAX_VARS];
    int n_mediators;
    bool is_valid;
} DCI_FrontdoorSet;

/* Back-door criterion */
DCI_BackdoorSet dci_backdoor_find(DCI_Graph* g, int cause, int effect);
bool dci_backdoor_check(DCI_Graph* g, int cause, int effect,
    const int* z_set, int n_z);
DCI_BackdoorSet dci_backdoor_minimal(DCI_Graph* g, int cause, int effect);

/* Back-door adjustment formula:
 * P(Y=y | do(X=x)) = Σ_z P(Y=y | X=x, Z=z) P(Z=z)
 */
double dci_backdoor_adjust(DCI_SCM* scm, DCI_Graph* g,
    int cause, int effect, double cause_val, double effect_val,
    DCI_BackdoorSet* z_set, int n_samples);

/* Front-door criterion */
DCI_FrontdoorSet dci_frontdoor_find(DCI_Graph* g, int cause, int effect);
bool dci_frontdoor_check(DCI_Graph* g, int cause, int effect,
    const int* mediators, int n_med);
double dci_frontdoor_adjust(DCI_SCM* scm, DCI_Graph* g,
    int cause, int effect, double cause_val, double effect_val,
    DCI_FrontdoorSet* med_set, int n_samples);

/* Covariate selection */
int dci_covariate_selection(DCI_Graph* g, int cause, int effect,
    int* sufficient_set, int max_size);

/* IV (Instrumental Variable) criterion */
bool dci_is_instrument(DCI_Graph* g, int instrument, int cause, int effect);

/* M-bias detection */
bool dci_detect_m_bias(DCI_Graph* g, int cause, int effect);

#endif
