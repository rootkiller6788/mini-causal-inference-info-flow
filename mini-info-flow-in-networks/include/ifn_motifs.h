#ifndef IFN_MOTIFS_H
#define IFN_MOTIFS_H
#include "ifn_core.h"

/* Network motifs: recurring, statistically significant subgraph patterns.
 * Motifs reveal functional building blocks of information processing.
 *
 * Triadic census: classification of all 3-node subgraphs in directed networks.
 * 16 possible triad types (0-15 per Davis & Leinhardt 1972).
 *
 * Ref: Milo et al. (2002) Science 298, 824-827
 *      Alon (2007) Nature Reviews Genetics 8, 450-461
 *      Holland & Leinhardt (1976) Sociological Methodology 7, 1-45
 */

/* Triadic census: count each of the 16 directed triad types.
 * Returns array of 16 counts (types 003 through 300 per standard notation). */
void ifn_triadic_census(const IFN_CausalGraph* g, int* triad_counts);

/* Z-score for motif significance: (count - mean_random) / std_random.
 * n_random: number of random graphs with same degree sequence. */
void ifn_motif_zscore(const IFN_CausalGraph* g, int* triad_counts,
    int n_random, double* z_scores);

/* Motif-based information flow: quantify how each triad type contributes
 * to total information routing through the network. */
void ifn_motif_information_flow(const IFN_CausalGraph* g,
    const double* node_signals, int n_vars, double* motif_flow);

/* Feed-forward loop (FFL) detection: pattern X->Y, X->Z, Y->Z.
 * Two types: coherent FFL (same sign) and incoherent FFL (mixed signs).
 * Fundamental building block of gene regulatory networks (Alon 2007). */
int ifn_count_feedforward_loops(const IFN_CausalGraph* g,
    double* coherent_ratio);

/* Bi-fan motif: X->A, X->B, Y->A, Y->B. Dense signaling pattern.
 * Common in neuronal and transcriptional networks. */
int ifn_count_bifan_motifs(const IFN_CausalGraph* g);

#endif
