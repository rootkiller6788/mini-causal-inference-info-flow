#ifndef PARTIAL_GRANGER_H
#define PARTIAL_GRANGER_H
#include "ts_core.h"

/* Partial Granger Causality and Partial Directed Coherence
 *
 * Partial Granger causality (PGC) extends conditional Granger by removing
 * the influence of exogenous inputs and latent confounding variables via
 * orthogonalization before computing the Granger test.
 *
 * Partial Directed Coherence (PDC) measures direct influence in the
 * frequency domain, normalized so that 0 <= PDC <= 1 and column sums = 1.
 *
 * Reference: Guo, Seth, Kendrick, Zhou, Feng (2008), J. Neurosci. Methods, 172(1):79-93.
 *            Baccala & Sameshima (2001), Biol. Cybern., 84(6):463-474.
 */

/* Partial Directed Coherence result: PDC_{i<-j}(omega) for all i,j at each frequency */
typedef struct {
    int n_vars;
    int n_freqs;
    double* frequencies;
    double** pdc;    /* [n_freqs][n_vars*n_vars] */
    int max_lag;
} PDCResult;

/* Partial Granger causality test result */
typedef struct {
    double f_statistic;
    double p_value;
    double rss_restricted;
    double rss_unrestricted;
    int lag_order;
    int n_obs;
    bool is_significant;
    double significance_level;
} PartialGrangerResult;

/* EEG/MEG frequency bands for band-averaged PDC */
typedef struct {
    double delta[16];   /* up to 4x4 */
    double theta[16];
    double alpha[16];
    double beta[16];
    double gamma[16];
    int dim;
} PDCBands;

/* Create and free PDC result.
 * Complexity: allocates O(n_freqs * n_vars^2). */
PDCResult* pdc_create(int dim, int n_freqs);
void pdc_free(PDCResult* p);

/* Compute PDC from a fitted VAR model at n_freqs uniformly spaced frequencies.
 * PDC_{i<-j}(omega) = |Ā_{ij}(omega)| / sqrt(sum_k |Ā_{kj}(omega)|^2)
 * where Ā(omega) = I - sum_{k=1}^p A_k e^{-i k omega}.
 * Complexity: O(n_freqs * dim^3). */
int pdc_compute(PDCResult* pdc, const VARModel* var, int n_freqs);

/* Print PDC matrix at a specific frequency index. */
void pdc_print(const PDCResult* pdc, int freq_idx);

/* Average PDC over standard EEG frequency bands:
 * delta (0.5-4 Hz), theta (4-8 Hz), alpha (8-13 Hz), beta (13-30 Hz), gamma (30+ Hz).
 * Normalized frequencies map proportionally to [0, pi]. */
PDCBands* pdc_band_average(const PDCResult* pdc);
void pdc_bands_free(PDCBands* b);
void pdc_bands_print(const PDCBands* b);

/* Partial Granger causality test.
 * Orthogonalizes X and Y with respect to exogenous Z, then computes
 * standard Granger F-test on the residuals.
 *
 * PGC_{X->Y|Z} removes Z's linear influence before testing.
 * More robust to latent confounding than conditional Granger.
 *
 * Complexity: O(length * max_lag^2). */
PartialGrangerResult* pg_create(void);
void pg_free(PartialGrangerResult* pg);
int pg_test(PartialGrangerResult* pg, const double* x, const double* y,
    const double* z, int length, int max_lag, double alpha);
void pg_print(const PartialGrangerResult* pg);

#endif
