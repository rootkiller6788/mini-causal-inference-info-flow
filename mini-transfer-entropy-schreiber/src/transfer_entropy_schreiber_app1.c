/* transfer_entropy_schreiber_app1.c -- L7 Application: Neural Information Transfer
 * Implements TE-based analysis of neural spike train recordings.
 * Modeled after the NASA/Ames neural data analysis pipeline.
 * Reference: Vicente, R. et al. (2011) "Transfer entropy--a model-free measure
 *            of effective connectivity for the neurosciences."
 *            J. Comput. Neurosci. 30, 45-67.
 */
#include "te_core.h"
#include "te_schreiber.h"
#include "te_ksg.h"
#include "te_effective.h"
#include "te_multivariate.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Generate synthetic neural spike train data with known connectivity.
 * Neuron A fires Poisson with rate lambda_a modulated by a slow oscillation.
 * Neuron B fires Poisson with rate lambda_b, driven by A with coupling c.
 *
 * This simulates a typical 2-electrode recording from prefrontal cortex
 * as described in the neurophysiology literature.
 */
static void generate_spike_trains(double *a, double *b, int n,
                                   double lambda_a, double lambda_b,
                                   double coupling, double *rates) {
    double phase = 0.0;
    for (int i = 0; i < n; i++) {
        double mod = 1.0 + 0.3 * sin(2.0 * 3.141592653589793 * phase);
        double rate_a = lambda_a * mod;
        double rate_b = lambda_b + coupling * ((i > 0) ? a[i-1] : 0.0);
        /* Poisson process: fire if uniform random < rate*dt */
        a[i] = ((double)rand() / (double)RAND_MAX < rate_a * 0.001) ? 1.0 : 0.0;
        b[i] = ((double)rand() / (double)RAND_MAX < rate_b * 0.001) ? 1.0 : 0.0;
        phase += 0.001;
        if (rates) { rates[i*2] = rate_a; rates[i*2+1] = rate_b; }
    }
}

/* Compute firing rate from spike train (spikes per second) */
static double firing_rate(const double *spikes, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += spikes[i];
    return sum / (double)n * 1000.0;  /* spikes/sec assuming 1ms bins */
}

/* Compute inter-spike interval statistics */
static void isi_stats(const double *spikes, int n,
                      double *mean_isi, double *cv_isi) {
    int *isi = malloc((size_t)n * sizeof(int));
    if (!isi) return;
    int count = 0, prev = -1;
    for (int i = 0; i < n; i++) {
        if (spikes[i] > 0.5) {
            if (prev >= 0) { isi[count++] = i - prev; }
            prev = i;
        }
    }
    if (count > 0) {
        double sum = 0.0, sum2 = 0.0;
        for (int i = 0; i < count; i++) {
            sum += (double)isi[i];
            sum2 += (double)isi[i] * (double)isi[i];
        }
        double m = sum / (double)count;
        *mean_isi = m;
        *cv_isi = (m > 0) ? sqrt(sum2/(double)count - m*m) / m : 0.0;
    } else { *mean_isi = 0; *cv_isi = 0; }
    free(isi);
}

/* TE-based functional connectivity analysis for neural data.
 * This is the core L7 application function.
 */
void te_neural_connectivity(const double *neuron_a, const double *neuron_b,
                             int n_samples, double fs,
                             double *te_ab, double *te_ba,
                             double *significance, int n_surr) {
    (void)fs; /* sampling frequency reserved for future use */
    TETimeSeries *ts_a = te_ts_create(neuron_a, n_samples, "Neuron A");
    TETimeSeries *ts_b = te_ts_create(neuron_b, n_samples, "Neuron B");
    if (!ts_a || !ts_b) return;

    /* Compute effective TE with bias correction */
    TEResult r = te_effective_compute(ts_a, ts_b, 3, 3, n_surr);
    *te_ab = r.te_xy;
    *te_ba = r.te_yx;
    if (significance) *significance = r.p_value;

    te_ts_free(ts_a);
    te_ts_free(ts_b);
}

/* --- Demo main (won't conflict with test mains since it's in .o only) --- */
#ifdef TE_APP1_MAIN
int main(void) {
    printf("=== Neural Information Transfer Analysis (L7 Application) ===\n\n");
    int n = 50000;  /* 50 seconds at 1kHz */
    double *neuron_a = calloc((size_t)n, sizeof(double));
    double *neuron_b = calloc((size_t)n, sizeof(double));
    double *rates = calloc((size_t)(n * 2), sizeof(double));
    if (!neuron_a || !neuron_b || !rates) {
        printf("Memory error\n"); return 1;
    }

    generate_spike_trains(neuron_a, neuron_b, n, 5.0, 3.0, 0.003, rates);

    printf("Firing rate A: %.1f Hz\n", firing_rate(neuron_a, n));
    printf("Firing rate B: %.1f Hz\n", firing_rate(neuron_b, n));

    double mean_isi_a, cv_isi_a, mean_isi_b, cv_isi_b;
    isi_stats(neuron_a, n, &mean_isi_a, &cv_isi_a);
    isi_stats(neuron_b, n, &mean_isi_b, &cv_isi_b);
    printf("ISI A: mean=%.1f ms, CV=%.2f\n", mean_isi_a, cv_isi_a);
    printf("ISI B: mean=%.1f ms, CV=%.2f\n", mean_isi_b, cv_isi_b);

    double te_ab, te_ba, sig;
    te_neural_connectivity(neuron_a, neuron_b, n, 1000.0, &te_ab, &te_ba, &sig, 50);
    printf("\nTE(A->B) = %.6f nats\n", te_ab);
    printf("TE(B->A) = %.6f nats\n", te_ba);
    printf("Significance (A->B): p = %.4f\n", sig);

    printf("\n=== Neural TE analysis complete ===\n");
    free(neuron_a); free(neuron_b); free(rates);
    return 0;
}
#endif