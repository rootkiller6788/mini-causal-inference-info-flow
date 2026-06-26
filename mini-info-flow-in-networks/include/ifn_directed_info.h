#ifndef IFN_DIRECTED_INFO_H
#define IFN_DIRECTED_INFO_H
#include "ifn_core.h"

/* Directed information theory (Massey 1990): I(X^n -> Y^n) = sum I(X^i; Y_i|Y^{i-1})
 * Generalizes mutual information to capture causal direction in time series.
 *
 * Key concepts:
 *   - Directed information: causal influence from past of X to present of Y
 *   - Causal conditioning: conditioning only on causal past
 *   - Kramer's causal conditioning: framing Granger causality in info-theory
 *   - Feedback capacity: maximum directed information through a channel
 *
 * Ref: Massey (1990) "Causality, feedback and directed information"
 *      Kramer (1998) "Directed Information for Channels with Feedback"
 *      Amblard & Michel (2013) Entropy 15(1), 113-141
 */

/* Directed information: I(X^T -> Y^T) = sum_{t=1}^{T} I(X^{t}; Y_t | Y^{t-1})
 * Measures cumulative causal influence of X on Y. */
double ifn_directed_information(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, int n_bins);

/* Transfer entropy as limiting form: TE_{X->Y} ≈ (1/T)*I(X^T -> Y^T) for large T */
double ifn_di_to_transfer_entropy(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, int n_bins);

/* Causally conditioned entropy: H(Y_t | Y^{t-1}, X^{t})
 * Key building block for directed information computation. */
double ifn_causally_conditioned_entropy(const double* history,
    const double* current, int T, int n_bins);

/* Feedback capacity: sup_{P(X|Y)} I(X -> Y) for a channel.
 * Upper bound on information transmission with feedback.
 * Approximated via Blahut-Arimoto with causal constraints. */
double ifn_feedback_capacity_bound(const double** channel, int n_inputs,
    int n_outputs, int n_steps, double tol, int max_iter);

#endif
