#include "ts_core.h"
#include "granger_test.h"
#include "transfer_entropy.h"
#include "spectral_granger.h"
#include "conditional_granger.h"
#include "fevd.h"
#include "nonlinear_granger.h"
#include "causality_graph.h"
#include "time_varying_granger.h"
#include "granger_bootstrap.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* L1-L2: Vector operations */
void t1(void) {
    Vector* v = vec_create(3);
    assert(v); v->data[0] = 3; v->data[1] = 4;
    assert(fabs(vec_norm(v) - 5.0) < 1e-9);
    vec_free(v);
    printf("  vec_norm: PASS\n");
}

void t2(void) {
    Matrix* A = mat_create(2, 2);
    mat_set(A, 0, 0, 1); mat_set(A, 1, 1, 4);
    assert(fabs(mat_trace(A) - 5.0) < 1e-9);
    mat_free(A);
    printf("  mat_trace: PASS\n");
}

void t3(void) {
    Matrix* I = mat_identity(3);
    assert(I); assert(fabs(mat_get(I, 1, 1) - 1.0) < 1e-9);
    assert(fabs(mat_get(I, 0, 1)) < 1e-9);
    mat_free(I);
    printf("  mat_identity: PASS\n");
}

/* L1: Time series creation and simulation */
void t4(void) {
    TimeSeries* ts = ts_simulate_ar1(0.5, 0.2, 100);
    assert(ts); assert(ts->length == 100); assert(ts->dim == 1);
    ts_free(ts);
    printf("  ar1_sim: PASS\n");
}

void t5(void) {
    VARModel* v = var_create(1, 2);
    assert(v); assert(v->dim == 1); assert(v->p == 2);
    var_free(v);
    printf("  var_create: PASS\n");
}

/* L3-L5: VAR fitting */
void t6(void) {
    TimeSeries* ts = ts_simulate_ar1(0.7, 0.1, 200);
    VARModel* v = var_create(1, 2);
    var_fit(v, ts);
    assert(v->sigma2 > 0);
    assert(fabs(v->A[0]->data[0] - 0.7) < 0.2);
    var_free(v); ts_free(ts);
    printf("  var_fit_ar1: PASS\n");
}

/* L5: Granger test creation */
void t7(void) {
    GrangerTestResult* gt = gt_create(); assert(gt);
    gt_free(gt);
    printf("  gt_create: PASS\n");
}

/* L2: Transfer entropy creation */
void t8(void) {
    TransferEntropy* te = te_create(); assert(te);
    te_free(te);
    printf("  te_create: PASS\n");
}

/* L2: Spectral Granger creation */
void t9(void) {
    SpectralGranger* sg = sg_create(10); assert(sg);
    assert(sg->n_freqs == 10);
    sg_free(sg);
    printf("  sg_create: PASS\n");
}

/* L2: Conditional Granger creation */
void t10(void) {
    ConditionalGranger* cg = cg_create(); assert(cg);
    cg_free(cg);
    printf("  cg_create: PASS\n");
}

/* L4-L5: Granger F-test on coupled VAR(1) */
void t11(void) {
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.3, 0.0, 0.5, 200);
    assert(ts);
    double* x = ts_extract_dim(ts, 0);
    double* y = ts_extract_dim(ts, 1);
    assert(x && y);
    GrangerTestResult gt;
    gt_test(&gt, x, y, 200, 5, 0.05);
    assert(gt.f_statistic > 0);
    free(x); free(y); ts_free(ts);
    printf("  granger_test: PASS (F=%.2f)\n", gt.f_statistic);
}

/* L5: Entropy estimation positivity */
void t12(void) {
    TimeSeries* ts = ts_simulate_ar1(0.5, 0.2, 100);
    double* x = ts_extract_dim(ts, 0);
    double h = entropy_histogram(x, 100, 10);
    assert(h > 0);
    free(x); ts_free(ts);
    printf("  entropy: PASS (%.4f)\n", h);
}

/* L5: Transfer entropy computation */
void t13(void) {
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.0, 0.3, 0.5, 100);
    double* x = ts_extract_dim(ts, 0);
    double* y = ts_extract_dim(ts, 1);
    TransferEntropy te;
    te_compute(&te, x, y, 100, 5, 3);
    assert(te.te_xy >= 0);
    free(x); free(y); ts_free(ts);
    printf("  te_compute: PASS (%.4f)\n", te.te_xy);
}

/* L5: VAR order selection (fixed: non-NULL ts) */
void t14(void) {
    TimeSeries* ts = ts_simulate_ar1(0.5, 0.2, 100);
    int best_p; double aic;
    int result = var_order_select(ts, 4, 1, &best_p, &aic);
    assert(result >= 1 && result <= 4);
    ts_free(ts);
    printf("  order_select: PASS (p=%d)\n", result);
}

/* L3: Multivariate VAR fitting */
void t15(void) {
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.3, 0.0, 0.5, 100);
    VARModel* v = var_create(2, 2);
    var_fit(v, ts);
    assert(v->sigma2 > 0);
    var_free(v); ts_free(ts);
    printf("  var_fit_2d: PASS\n");
}

/* L1: Time series slicing */
void t16(void) {
    TimeSeries* ts = ts_simulate_ar1(0.5, 0.2, 50);
    TimeSeries* s = ts_slice(ts, 10, 20);
    assert(s); assert(s->length == 20);
    ts_free(s); ts_free(ts);
    printf("  ts_slice: PASS\n");
}

/* L5: Spectral Granger computation */
void t17(void) {
    SpectralGranger* sg = sg_create(20);
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.3, 0.0, 0.5, 50);
    double* x = ts_extract_dim(ts, 0);
    double* y = ts_extract_dim(ts, 1);
    sg_compute(sg, x, y, 50, 20, 3);
    assert(sg->n_freqs == 20);
    free(x); free(y); ts_free(ts); sg_free(sg);
    printf("  sg_compute: PASS\n");
}

/* L3: Cholesky decomposition */
void t18(void) {
    Matrix* A = mat_create(3, 3);
    mat_set(A, 0, 0, 4); mat_set(A, 0, 1, 2); mat_set(A, 0, 2, 1);
    mat_set(A, 1, 0, 2); mat_set(A, 1, 1, 5); mat_set(A, 1, 2, 3);
    mat_set(A, 2, 0, 1); mat_set(A, 2, 1, 3); mat_set(A, 2, 2, 6);
    Matrix* L = cholesky_decompose(A);
    assert(L != NULL);
    /* Verify L * L^T recovers A */
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double sum = 0;
            for (int k = 0; k < 3; k++)
                sum += L->data[i*3+k] * L->data[j*3+k];
            assert(fabs(sum - A->data[i*3+j]) < 1e-6);
        }
    mat_free(A); mat_free(L);
    printf("  cholesky: PASS\n");
}

/* L3: IRF computation */
void t19(void) {
    VARModel* var = var_create(2, 1);
    mat_set(var->A[0], 0, 0, 0.5); mat_set(var->A[0], 0, 1, 0.1);
    mat_set(var->A[0], 1, 0, 0.0); mat_set(var->A[0], 1, 1, 0.5);
    Matrix** irf = compute_irf(var, 3);
    assert(irf != NULL);
    /* Phi_0 = I */
    assert(fabs(irf[0]->data[0] - 1.0) < 1e-9);
    assert(fabs(irf[0]->data[3] - 1.0) < 1e-9);
    for (int h = 0; h < 3; h++) mat_free(irf[h]);
    free(irf); var_free(var);
    printf("  compute_irf: PASS\n");
}

/* L5: FEVD computation */
void t20(void) {
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.2, 0.0, 0.5, 200);
    VARModel* var = var_create(2, 2);
    var_fit(var, ts);
    FEVDResult* fevd = fevd_compute(var, 5);
    assert(fevd != NULL);
    assert(fevd->n_vars == 2);
    assert(fevd->horizon == 5);
    /* Each row should sum to approximately 1 */
    for (int h = 0; h < 5; h++) {
        double row0 = fevd->matrix[h][0] + fevd->matrix[h][1];
        double row1 = fevd->matrix[h][2] + fevd->matrix[h][3];
        assert(fabs(row0 - 1.0) < 1e-6);
        assert(fabs(row1 - 1.0) < 1e-6);
    }
    fevd_free(fevd); var_free(var); ts_free(ts);
    printf("  fevd_compute: PASS\n");
}

/* L8: Nonlinear Granger */
void t21(void) {
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.3, 0.0, 0.5, 100);
    double* x = ts_extract_dim(ts, 0);
    double* y = ts_extract_dim(ts, 1);
    NonlinearGrangerResult* nlg = nlg_create();
    assert(nlg != NULL);
    double sigma = kernel_median_heuristic(x, 100);
    assert(sigma > 0);
    nlg_test(nlg, x, y, 100, 3, sigma, 0.05);
    assert(nlg->n_obs > 0);
    nlg_free(nlg); free(x); free(y); ts_free(ts);
    printf("  nlg_test: PASS\n");
}

/* L5: Causality graph construction */
void t22(void) {
    TimeSeries* ts = ts_create(50, 3);
    for (int t = 0; t < 50; t++) {
        ts->values[t*3+0] = 0.5*((t>0?ts->values[(t-1)*3+0]:0.0)) + 0.1*((double)rand()/RAND_MAX-0.5);
        ts->values[t*3+1] = 0.3*((t>0?ts->values[(t-1)*3+0]:0.0)) + 0.4*((t>0?ts->values[(t-1)*3+1]:0.0)) + 0.1*((double)rand()/RAND_MAX-0.5);
        ts->values[t*3+2] = 0.5*((t>0?ts->values[(t-1)*3+2]:0.0)) + 0.1*((double)rand()/RAND_MAX-0.5);
    }
    CausalityGraph* cg = graph_build(ts, 2, 0.05);
    assert(cg != NULL);
    assert(cg->n_vars == 3);
    graph_free(cg); ts_free(ts);
    printf("  graph_build: PASS\n");
}

/* L8: Bootstrap Granger */
void t23(void) {
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.3, 0.0, 0.5, 100);
    double* x = ts_extract_dim(ts, 0);
    double* y = ts_extract_dim(ts, 1);
    BootstrapGrangerResult* boot = boot_create(50);
    assert(boot != NULL);
    boot_parametric(boot, x, y, 100, 2, 0.05);
    assert(boot->observed_F >= 0);
    assert(boot->bootstrap_p_value >= 0 && boot->bootstrap_p_value <= 1);
    boot_free(boot); free(x); free(y); ts_free(ts);
    printf("  boot_parametric: PASS\n");
}

/* L8: Time-varying Granger */
void t24(void) {
    double* x = malloc(200 * sizeof(double));
    double* y = malloc(200 * sizeof(double));
    x[0] = 0; y[0] = 0;
    for (int t = 1; t < 200; t++) {
        x[t] = 0.5*x[t-1] + 0.1*((double)rand()/RAND_MAX-0.5);
        y[t] = 0.3*x[t-1] + 0.5*y[t-1] + 0.1*((double)rand()/RAND_MAX-0.5);
    }
    TVGrangerResult* tvg = tvg_create();
    assert(tvg != NULL);
    tvg_compute(tvg, x, y, 200, 50, 10, 3, 0.05, WINDOW_SLIDING);
    assert(tvg->n_windows > 0);
    tvg_free(tvg); free(x); free(y);
    printf("  tvg_compute: PASS\n");
}

/* L5: Matrix multiplication */
void t25(void) {
    Matrix* A = mat_create(2, 3);
    Matrix* B = mat_create(3, 2);
    for (int i = 0; i < 6; i++) A->data[i] = i + 1;
    for (int i = 0; i < 6; i++) B->data[i] = (i % 2 == 0) ? 1 : 0;
    Matrix* C = mat_mul(A, B);
    assert(C != NULL);
    assert(C->rows == 2 && C->cols == 2);
    mat_free(A); mat_free(B); mat_free(C);
    printf("  mat_mul: PASS\n");
}

/* L5: Matrix determinant and inverse consistency */
void t26(void) {
    Matrix* A = mat_create(2, 2);
    mat_set(A, 0, 0, 4); mat_set(A, 0, 1, 7);
    mat_set(A, 1, 0, 2); mat_set(A, 1, 1, 6);
    double det = mat_det(A);
    assert(fabs(det - 10.0) < 1e-6);
    Matrix* inv = mat_inverse(A);
    assert(inv != NULL);
    Matrix* prod = mat_mul(A, inv);
    assert(fabs(prod->data[0] - 1.0) < 1e-6);
    assert(fabs(prod->data[3] - 1.0) < 1e-6);
    mat_free(A); mat_free(inv); mat_free(prod);
    printf("  mat_inverse: PASS\n");
}

int main(void) {
    printf("=== Granger Causality Test Suite ===\n"); fflush(stdout);
    t1(); t2(); t3(); t4(); t5(); t6(); t7(); t8(); t9(); t10();
    t11(); t12(); t13(); t14(); t15(); t16(); t17(); t18(); t19();
    t20(); t21(); t22(); t23(); t24(); t25(); t26();
    printf("=== All 26 tests passed ===\n");
    return 0;
}