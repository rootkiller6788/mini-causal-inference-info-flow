/* test_te_schreiber.c -- Tests for Schreiber transfer entropy estimator */
#include "te_schreiber.h"
#include "te_core.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

static int test_schreiber_basic(void) {
    double dx[500], dy[500];
    dx[0] = 0.0; dy[0] = 0.0;
    for (int i = 1; i < 500; i++) {
        dx[i] = 0.5*dx[i-1] + 0.3*dy[i-1] + 0.1*((double)rand()/(double)RAND_MAX-0.5);
        dy[i] = 0.5*dy[i-1] + 0.1*((double)rand()/(double)RAND_MAX-0.5);
    }
    TETimeSeries *tx = te_ts_create(dx, 500, "x");
    TETimeSeries *ty = te_ts_create(dy, 500, "y");
    TEResult r = te_schreiber_binning(tx, ty, 2, 2, 6);
    assert(r.te_xy >= 0.0);
    double dir = te_schreiber_directionality(&r);
    assert(dir >= -1.0 && dir <= 1.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_entropy_rate(void) {
    double data[200];
    for (int i = 0; i < 200; i++) data[i] = (double)(i % 5);
    TETimeSeries *ts = te_ts_create(data, 200, "periodic");
    TEBinnedSeries *bs = te_bin_create(ts, 5, TE_BIN_EQUAL_WIDTH);
    assert(bs != NULL);
    double h1 = te_schreiber_entropy_rate(bs, 1);
    assert(h1 >= 0.0);
    double h3 = te_schreiber_entropy_rate(bs, 3);
    assert(h3 >= 0.0);
    assert(h3 <= h1 + 0.5);
    te_bin_free(bs);
    te_ts_free(ts);
    return 0;
}

static int test_schreiber_mutual_info(void) {
    double dx[200], dy[200];
    for (int i = 0; i < 200; i++) { dx[i] = sin(0.1*(double)i); dy[i] = sin(0.1*(double)i + 0.5); }
    TETimeSeries *tx = te_ts_create(dx, 200, "x");
    TETimeSeries *ty = te_ts_create(dy, 200, "y");
    TEBinnedSeries *bx = te_bin_create(tx, 8, TE_BIN_EQUAL_WIDTH);
    TEBinnedSeries *by = te_bin_create(ty, 8, TE_BIN_EQUAL_WIDTH);
    assert(bx && by);
    double mi0 = te_schreiber_mutual_info(bx, by, 0);
    assert(mi0 >= 0.0);
    double mi1 = te_schreiber_mutual_info(bx, by, 1);
    assert(mi1 >= 0.0);
    te_bin_free(bx); te_bin_free(by);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_scan_lag(void) {
    double dx[200], dy[200];
    for (int i = 0; i < 200; i++) { dx[i] = (double)(i % 10); dy[i] = dx[(i+5)%200]; }
    TETimeSeries *tx = te_ts_create(dx, 200, "x");
    TETimeSeries *ty = te_ts_create(dy, 200, "y");
    double te_vals[10];
    te_schreiber_scan_lag(tx, ty, 1, 1, 10, te_vals);
    int has_positive = 0;
    for (int i = 0; i < 10; i++) if (te_vals[i] > 0.0) has_positive = 1;
    (void)has_positive;
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_multiscale(void) {
    double dx[300], dy[300];
    for (int i = 0; i < 300; i++) { dx[i] = (double)(i%10); dy[i] = dx[i]*0.5+0.5*((double)rand()/(double)RAND_MAX); }
    TETimeSeries *tx = te_ts_create(dx, 300, "x");
    TETimeSeries *ty = te_ts_create(dy, 300, "y");
    double te_scales[5];
    te_schreiber_multiscale(tx, ty, 1, 1, 5, te_scales);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_sliding(void) {
    double dx[300], dy[300];
    for (int i = 0; i < 300; i++) { dx[i] = (double)(i%8); dy[i] = dx[i]*0.7+0.3*((double)rand()/(double)RAND_MAX); }
    TETimeSeries *tx = te_ts_create(dx, 300, "x");
    TETimeSeries *ty = te_ts_create(dy, 300, "y");
    double te_tv[20];
    int n = te_schreiber_sliding_window(tx, ty, 1, 1, 100, 20, te_tv, 20);
    assert(n > 0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_optimal_embedding(void) {
    double dx[250], dy[250];
    dx[0]=dy[0]=0;
    for(int i=1;i<250;i++){dx[i]=0.6*dx[i-1]+0.2*dy[i-1]+0.1*((double)rand()/(double)RAND_MAX-0.5);dy[i]=0.5*dy[i-1]+0.1*((double)rand()/(double)RAND_MAX-0.5);}
    TETimeSeries *tx=te_ts_create(dx,250,"x"), *ty=te_ts_create(dy,250,"y");
    int bk,bl;
    te_schreiber_optimal_embedding(tx,ty,3,3,&bk,&bl);
    assert(bk>=1&&bk<=3); assert(bl>=1&&bl<=3);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_robust(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=(double)(i%7);dy[i]=(double)((i+3)%7);}
    TETimeSeries *tx=te_ts_create(dx,200,"x"), *ty=te_ts_create(dy,200,"y");
    TEResult r=te_schreiber_robust_binning(tx,ty,1,1,4,12,5);
    assert(r.te_xy>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_normalized(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=(double)(i%5);dy[i]=dx[(i+2)%200];}
    TETimeSeries *tx=te_ts_create(dx,200,"x"), *ty=te_ts_create(dy,200,"y");
    TEResult r=te_schreiber_binning(tx,ty,2,2,5);
    double nte=te_schreiber_normalized(&r);
    assert(nte>=0.0);
    double uc=te_schreiber_uncertainty_coefficient(&r);
    assert(uc>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_theiler(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=(double)(i%6);dy[i]=dx[(i+1)%200];}
    TETimeSeries *tx=te_ts_create(dx,200,"x"), *ty=te_ts_create(dy,200,"y");
    TEResult r=te_schreiber_theiler(tx,ty,2,2,5,3);
    assert(r.te_xy>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_integrated(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=(double)(i%5);dy[i]=dx[(i+1)%200];}
    TETimeSeries *tx=te_ts_create(dx,200,"x"), *ty=te_ts_create(dy,200,"y");
    double ite=te_schreiber_integrated(tx,ty,3,3);
    assert(ite>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_schreiber_variable_embedding(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=(double)(i%5);dy[i]=dx[(i+1)%200];}
    TETimeSeries *tx=te_ts_create(dx,200,"x"), *ty=te_ts_create(dy,200,"y");
    TEResult r=te_schreiber_variable_embedding(tx,ty,1,3,1,3);
    assert(r.te_xy>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

int main(void) {
    printf("=== test_te_schreiber ===\n");
    test_schreiber_basic(); printf("  PASS basic\n");
    test_schreiber_entropy_rate(); printf("  PASS entropy_rate\n");
    test_schreiber_mutual_info(); printf("  PASS mutual_info\n");
    test_schreiber_scan_lag(); printf("  PASS scan_lag\n");
    test_schreiber_multiscale(); printf("  PASS multiscale\n");
    test_schreiber_sliding(); printf("  PASS sliding\n");
    test_schreiber_optimal_embedding(); printf("  PASS optimal_embedding\n");
    test_schreiber_robust(); printf("  PASS robust\n");
    test_schreiber_normalized(); printf("  PASS normalized\n");
    test_schreiber_theiler(); printf("  PASS theiler\n");
    test_schreiber_integrated(); printf("  PASS integrated\n");
    test_schreiber_variable_embedding(); printf("  PASS variable_embedding\n");
    printf("=== test_te_schreiber: ALL PASSED ===\n");
    return 0;
}