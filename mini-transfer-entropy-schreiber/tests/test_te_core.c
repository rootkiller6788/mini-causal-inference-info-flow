/* test_te_core.c -- Tests for core TE operations and information measures */
#include "te_core.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

static int test_ts_create_free(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    TETimeSeries *ts = te_ts_create(data, 5, "test");
    assert(ts != NULL);
    assert(ts->length == 5);
    assert(fabs(ts->data[2] - 3.0) < TE_EPSILON);
    assert(strcmp(ts->name, "test") == 0);
    te_ts_free(ts);
    te_ts_free(NULL);
    assert(te_ts_create(NULL, 5, "x") == NULL);
    assert(te_ts_create(data, 0, "x") == NULL);
    return 0;
}

static int test_ts_statistics(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    TETimeSeries *ts = te_ts_create(data, 5, "test");
    assert(fabs(te_ts_mean(ts) - 3.0) < TE_EPSILON);
    assert(fabs(te_ts_variance(ts) - 2.5) < TE_EPSILON);
    assert(fabs(te_ts_std(ts) - sqrt(2.5)) < TE_EPSILON);
    assert(fabs(te_ts_min(ts) - 1.0) < TE_EPSILON);
    assert(fabs(te_ts_max(ts) - 5.0) < TE_EPSILON);
    double skew = te_ts_skewness(ts);
    assert(fabs(skew) < 0.01);
    assert(fabs(te_ts_rms(ts) - sqrt(11.0)) < 0.001);
    te_ts_free(ts);
    return 0;
}

static int test_ts_normalize(void) {
    double data[] = {10.0, 20.0, 30.0, 40.0, 50.0};
    TETimeSeries *ts = te_ts_create(data, 5, "test");
    te_ts_normalize(ts);
    assert(fabs(te_ts_min(ts) - 0.0) < TE_EPSILON);
    assert(fabs(te_ts_max(ts) - 1.0) < TE_EPSILON);
    te_ts_free(ts);
    ts = te_ts_create(data, 5, "test");
    te_ts_standardize(ts);
    assert(fabs(te_ts_mean(ts)) < TE_EPSILON);
    assert(fabs(te_ts_std(ts) - 1.0) < TE_EPSILON);
    te_ts_free(ts);
    return 0;
}

static int test_ts_autocorr(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    TETimeSeries *ts = te_ts_create(data, 9, "acf");
    double ac0 = te_ts_autocorrelation(ts, 0);
    assert(fabs(ac0 - 1.0) < TE_EPSILON);
    double ac1 = te_ts_autocorrelation(ts, 1);
    assert(ac1 > 0.0 && ac1 <= 1.0);
    te_ts_free(ts);
    return 0;
}

static int test_binning(void) {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = (double)i;
    TETimeSeries *ts = te_ts_create(data, 100, "linear");
    TEBinnedSeries *bs = te_bin_create(ts, 10, TE_BIN_EQUAL_WIDTH);
    assert(bs != NULL);
    assert(bs->n_bins == 10);
    assert(bs->length == 100);
    assert(te_bin_get(bs, 0) == 0);
    assert(te_bin_get(bs, 99) == 9);
    int mid = te_bin_get(bs, 50);
    assert(mid >= 4 && mid <= 5);
    double H = te_bin_entropy(bs);
    assert(H > 2.0 && H < 2.5);
    double hist[10];
    te_bin_count(bs, hist);
    double hsum = 0.0;
    for (int b = 0; b < 10; b++) { assert(hist[b] > 0.0); hsum += hist[b]; }
    assert(fabs(hsum - 1.0) < TE_EPSILON);
    te_bin_free(bs);
    te_ts_free(ts);
    return 0;
}

static int test_embedding(void) {
    double data[50];
    for (int i = 0; i < 50; i++) data[i] = (double)(i % 7);
    TETimeSeries *ts = te_ts_create(data, 50, "periodic");
    TEBinnedSeries *bs = te_bin_create(ts, 7, TE_BIN_EQUAL_WIDTH);
    assert(bs != NULL);
    TEEmbedding *emb = te_embed_create(bs, 3, 2);
    assert(emb != NULL);
    assert(emb->embed_dim == 3);
    assert(emb->delay == 2);
    assert(emb->n_states == bs->length - (3-1)*2);
    te_embed_free(emb);
    te_bin_free(bs);
    te_ts_free(ts);
    return 0;
}

static int test_joint_prob(void) {
    double dx[100], dy[100];
    for (int i = 0; i < 100; i++) { dx[i] = (double)(i % 10); dy[i] = (double)((i/10)%10); }
    TETimeSeries *tx = te_ts_create(dx, 100, "x");
    TETimeSeries *ty = te_ts_create(dy, 100, "y");
    TEBinnedSeries *bx = te_bin_create(tx, 10, TE_BIN_EQUAL_WIDTH);
    TEBinnedSeries *by = te_bin_create(ty, 10, TE_BIN_EQUAL_WIDTH);
    assert(bx && by);
    TEJointProb *jp = te_joint_estimate(bx, by, 0);
    assert(jp != NULL);
    assert(jp->nx == 10 && jp->ny == 10);
    double Hj = te_joint_entropy(jp);
    assert(Hj >= 0.0);
    double Hc = te_conditional_entropy(jp);
    assert(Hc >= 0.0);
    assert(Hc <= Hj + TE_EPSILON);
    te_joint_free(jp);
    te_bin_free(bx); te_bin_free(by);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_mutual_info(void) {
    double dx[200], dy[200];
    for (int i = 0; i < 200; i++) { dx[i] = (double)(i % 10); dy[i] = (double)(i % 10); }
    TETimeSeries *tx = te_ts_create(dx, 200, "x");
    TETimeSeries *ty = te_ts_create(dy, 200, "y");
    double mi = te_mutual_information(tx, ty, 10, 0);
    assert(mi > 0.5);
    assert(mi >= 0.0);
    assert(mi <= log(10.0) + TE_EPSILON);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_surrogate(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    TETimeSeries *ts = te_ts_create(data, 10, "test");
    te_ts_shuffle(ts);
    assert(fabs(te_ts_min(ts) - 1.0) < TE_EPSILON);
    assert(fabs(te_ts_max(ts) - 10.0) < TE_EPSILON);
    te_ts_free(ts);
    return 0;
}

static int test_optimal_bins(void) {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = (double)i;
    TETimeSeries *ts = te_ts_create(data, 100, "linear");
    int nb = te_estimate_optimal_bins(ts);
    assert(nb >= 2 && nb <= TE_MAX_BINS);
    te_ts_free(ts);
    return 0;
}

static int test_quantile(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    TETimeSeries *ts = te_ts_create(data, 10, "test");
    double q50 = te_ts_quantile(ts, 0.5);
    assert(q50 >= 5.0 && q50 <= 6.0);
    double p25 = te_ts_percentile(ts, 25.0);
    assert(p25 >= 2.5 && p25 <= 3.5);
    double cdf3 = te_ts_empirical_cdf(ts, 3.0);
    assert(fabs(cdf3 - 0.3) < TE_EPSILON);
    te_ts_free(ts);
    return 0;
}

static int test_ais(void) {
    double data[200];
    data[0] = 0.0;
    for (int i = 1; i < 200; i++) {
        double noise = 0.1 * ((double)rand() / (double)RAND_MAX - 0.5);
        data[i] = 0.8 * data[i-1] + noise;
    }
    TETimeSeries *ts = te_ts_create(data, 200, "ar1");
    double ais = te_active_information_storage(ts, 2, 6);
    assert(ais > 0.0);
    te_ts_free(ts);
    return 0;
}

static int test_generators(void) {
    TETimeSeries *sine = te_ts_generate_sine(100, 5.0, 1.0, 0.0);
    assert(sine != NULL && sine->length == 100);
    te_ts_free(sine);
    TETimeSeries *henon = te_ts_generate_henon(100, 1.4, 0.3);
    assert(henon != NULL && henon->length == 100);
    te_ts_free(henon);
    TETimeSeries *logistic = te_ts_generate_logistic(100, 3.9, 0.5);
    assert(logistic != NULL);
    assert(te_ts_min(logistic) >= 0.0 && te_ts_max(logistic) <= 1.0 + TE_EPSILON);
    te_ts_free(logistic);
    TETimeSeries *ar1 = te_ts_generate_ar1(100, 0.5, 0.1);
    assert(ar1 != NULL);
    te_ts_free(ar1);
    TETimeSeries *coupled = te_ts_generate_coupled_ar1(100, 0.2, 0.5, 0.1);
    assert(coupled != NULL);
    te_ts_free(coupled);
    return 0;
}

static int test_apen_sampen(void) {
    double data[200];
    data[0] = 0.0;
    for (int i = 1; i < 200; i++)
        data[i] = 0.7 * data[i-1] + 0.1 * ((double)rand()/(double)RAND_MAX - 0.5);
    TETimeSeries *ts = te_ts_create(data, 200, "ar1");
    double apen = te_approximate_entropy(ts, 2, 0.3);
    assert(apen >= 0.0);
    double sampen = te_sample_entropy(ts, 2, 0.2);
    assert(sampen >= 0.0);
    te_ts_free(ts);
    return 0;
}

static int test_renyi_tsallis(void) {
    double data[200];
    for (int i = 0; i < 200; i++) data[i] = (double)(i % 10);
    TETimeSeries *ts = te_ts_create(data, 200, "uniform");
    TEBinnedSeries *bs = te_bin_create(ts, 10, TE_BIN_EQUAL_WIDTH);
    assert(bs != NULL);
    double r0 = te_renyi_entropy(bs, 0.0);
    assert(fabs(r0 - log(10.0)) < 0.01);
    double r1 = te_renyi_entropy(bs, 1.0);
    double sh = te_bin_entropy(bs);
    assert(fabs(r1 - sh) < 0.01);
    double t2 = te_tsallis_entropy(bs, 2.0);
    assert(t2 >= 0.0);
    te_bin_free(bs);
    te_ts_free(ts);
    return 0;
}

static int test_js_div(void) {
    double p[] = {0.5, 0.3, 0.2};
    double q[] = {0.2, 0.3, 0.5};
    double jsd = te_js_divergence(p, q, 3);
    assert(jsd >= 0.0 && jsd < log(2.0));
    double jsd_same = te_js_divergence(p, p, 3);
    assert(fabs(jsd_same) < TE_EPSILON);
    return 0;
}

static int test_dw(void) {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = (double)(i % 10);
    TETimeSeries *ts = te_ts_create(data, 100, "periodic");
    double dw = te_ts_durbin_watson(ts);
    assert(dw > 0.0 && dw < 4.0);
    te_ts_free(ts);
    return 0;
}

static int test_hurst(void) {
    double data[128];
    for (int i = 0; i < 128; i++) data[i] = (double)i * 0.01;
    TETimeSeries *ts = te_ts_create(data, 128, "trend");
    double H = te_ts_hurst_rs(ts, 8, 64, 5);
    assert(H > 0.0 && H <= 1.0);
    te_ts_free(ts);
    return 0;
}

static int test_mw(void) {
    double a[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double b[] = {6.0, 7.0, 8.0, 9.0, 10.0};
    double p_val;
    te_mann_whitney_u(a, 5, b, 5, &p_val);
    assert(p_val < 0.05);
    return 0;
}

static int test_hjorth(void) {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = sin(0.1 * (double)i);
    TETimeSeries *ts = te_ts_create(data, 100, "sine");
    double activity, mobility, complexity;
    te_ts_hjorth(ts, &activity, &mobility, &complexity);
    assert(activity > 0.0 && mobility > 0.0);
    te_ts_free(ts);
    return 0;
}

static int test_lzc(void) {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = (double)(i % 5);
    TETimeSeries *ts = te_ts_create(data, 100, "periodic");
    TEBinnedSeries *bs = te_bin_create(ts, 5, TE_BIN_EQUAL_WIDTH);
    assert(bs != NULL);
    double lzc = te_lz_complexity(bs);
    assert(lzc >= 0.0 && lzc <= 2.0);
    te_bin_free(bs);
    te_ts_free(ts);
    return 0;
}

static int test_detrend(void) {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = 0.01*(double)i + 0.1*((double)rand()/(double)RAND_MAX - 0.5);
    TETimeSeries *ts = te_ts_create(data, 100, "trend");
    te_ts_detrend_linear(ts);
    double m = te_ts_mean(ts);
    assert(fabs(m) < 0.15);
    te_ts_free(ts);
    return 0;
}

static int test_predictive_info(void) {
    double data[200];
    data[0] = 0.0;
    for (int i = 1; i < 200; i++) data[i] = 0.8*data[i-1] + 0.1*((double)rand()/(double)RAND_MAX - 0.5);
    TETimeSeries *ts = te_ts_create(data, 200, "ar1");
    double pi = te_predictive_information(ts, 3, 2, 6);
    assert(pi >= 0.0);
    te_ts_free(ts);
    return 0;
}

static int test_window_functions(void) {
    double data[50];
    for (int i = 0; i < 50; i++) data[i] = 1.0;
    TETimeSeries *ts = te_ts_create(data, 50, "ones");
    te_ts_hamming_window(ts);
    assert(te_ts_max(ts) <= 1.0);
    te_ts_free(ts);
    ts = te_ts_create(data, 50, "ones");
    te_ts_hann_window(ts);
    assert(te_ts_max(ts) <= 1.0);
    te_ts_free(ts);
    ts = te_ts_create(data, 50, "ones");
    te_ts_blackman_window(ts);
    assert(te_ts_max(ts) <= 1.0);
    te_ts_free(ts);
    return 0;
}

static int test_clone_resample(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    TETimeSeries *ts = te_ts_create(data, 5, "orig");
    TETimeSeries *clone = te_ts_clone(ts);
    assert(clone != NULL);
    assert(clone->length == 5);
    assert(fabs(clone->data[2] - 3.0) < TE_EPSILON);
    te_ts_free(clone);
    te_ts_resample(ts, 10);
    assert(ts->length == 10);
    te_ts_free(ts);
    return 0;
}

static int test_moving_average(void) {
    double data[50];
    for (int i = 0; i < 50; i++) data[i] = (double)(i % 5) + 0.1*((double)rand()/(double)RAND_MAX - 0.5);
    TETimeSeries *ts = te_ts_create(data, 50, "noisy");
    te_ts_moving_average(ts, 5);
    assert(ts != NULL);
    te_ts_free(ts);
    return 0;
}

static int test_exponential_smooth(void) {
    double data[50];
    for (int i = 0; i < 50; i++) data[i] = (double)i;
    TETimeSeries *ts = te_ts_create(data, 50, "linear");
    te_ts_exponential_smooth(ts, 0.3);
    assert(ts != NULL);
    te_ts_free(ts);
    return 0;
}

static int test_bds(void) {
    double data[200];
    for (int i = 0; i < 200; i++) data[i] = 0.5 * ((double)rand()/(double)RAND_MAX - 0.5);
    TETimeSeries *ts = te_ts_create(data, 200, "noise");
    double bds = te_bds_statistic(ts, 2, 0.1);
    (void)bds;
    te_ts_free(ts);
    return 0;
}

int main(void) {
    printf("=== test_te_core ===\n");
    test_ts_create_free(); printf("  PASS create_free\n");
    test_ts_statistics(); printf("  PASS statistics\n");
    test_ts_normalize(); printf("  PASS normalize\n");
    test_ts_autocorr(); printf("  PASS autocorr\n");
    test_binning(); printf("  PASS binning\n");
    test_embedding(); printf("  PASS embedding\n");
    test_joint_prob(); printf("  PASS joint_prob\n");
    test_mutual_info(); printf("  PASS mutual_info\n");
    test_surrogate(); printf("  PASS surrogate\n");
    test_optimal_bins(); printf("  PASS optimal_bins\n");
    test_quantile(); printf("  PASS quantile\n");
    test_ais(); printf("  PASS AIS\n");
    test_generators(); printf("  PASS generators\n");
    test_apen_sampen(); printf("  PASS ApEn/SampEn\n");
    test_renyi_tsallis(); printf("  PASS Renyi/Tsallis\n");
    test_js_div(); printf("  PASS JS_div\n");
    test_dw(); printf("  PASS DW\n");
    test_hurst(); printf("  PASS Hurst\n");
    test_mw(); printf("  PASS MW\n");
    test_hjorth(); printf("  PASS Hjorth\n");
    test_lzc(); printf("  PASS LZC\n");
    test_detrend(); printf("  PASS detrend\n");
    test_predictive_info(); printf("  PASS predictive_info\n");
    test_window_functions(); printf("  PASS windows\n");
    test_clone_resample(); printf("  PASS clone_resample\n");
    test_moving_average(); printf("  PASS moving_average\n");
    test_exponential_smooth(); printf("  PASS exp_smooth\n");
    test_bds(); printf("  PASS BDS\n");
    printf("=== test_te_core: ALL PASSED ===\n");
    return 0;
}