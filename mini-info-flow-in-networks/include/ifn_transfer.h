#ifndef IFN_TRANSFER_H
#define IFN_TRANSFER_H
#include "ifn_core.h"

/* Transfer entropy: TE(X->Y) = I(Y_{t+1}; X_t | Y_t)
 * Quantifies directed information flow from X to Y. */

IFN_TransferResult ifn_transfer_entropy(const IFN_TimeSeries* source,
    const IFN_TimeSeries* target, int k_history, int l_history);

IFN_TransferResult ifn_transfer_entropy_discrete(const int* src, const int* tgt,
    int T, int n_src_states, int n_tgt_states, int k, int l);

double ifn_effective_transfer(const IFN_TimeSeries* source,
    const IFN_TimeSeries* target, int k, int l);

void ifn_transfer_significance(IFN_TransferResult* result,
    int n_surrogates, double alpha);

double ifn_transfer_entropy_normalized(const IFN_TimeSeries* source,
    const IFN_TimeSeries* target, int k, int l);

void ifn_transfer_matrix(const IFN_TimeSeries* multi_ts,
    int n_vars, int k, int l, double* te_matrix);

/* Symbolic transfer entropy (for noisy data) */
double ifn_symbolic_transfer_entropy(const double* src, const double* tgt,
    int T, int symbol_size, int k, int l);

#endif
