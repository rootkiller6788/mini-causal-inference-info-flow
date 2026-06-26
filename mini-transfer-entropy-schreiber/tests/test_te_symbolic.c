/* test_te_symbolic.c -- Tests for symbolic/permutation transfer entropy */
#include "te_symbolic.h"
#include "te_effective.h"
#include "te_core.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

static int test_permutation_encode(void) {
    double w[]={3.0,1.0,4.0,2.0};
    int pat[4];
    int idx=te_symbolic_permutation_encode(w,4,pat);
    assert(idx>=0);
    assert(idx<24);
    assert(pat[0]==2);
    assert(pat[1]==0);
    assert(pat[2]==3);
    assert(pat[3]==1);
    return 0;
}

static int test_pattern_count(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.1*(double)i)+0.05*((double)rand()/(double)RAND_MAX-0.5);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    int np=1; for(int i=2;i<=4;i++) np*=i;
    int *counts=calloc(np,sizeof(int));
    te_symbolic_pattern_count(ts,4,2,counts,np);
    int total=0;
    for(int i=0;i<np;i++) total+=counts[i];
    assert(total>0);
    free(counts);
    te_ts_free(ts);
    return 0;
}

static int test_permutation_entropy(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.1*(double)i);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    double pe=te_symbolic_permutation_entropy(ts,4,2);
    assert(pe>0.0);
    double pe_norm=pe/log(24.0);
    assert(pe_norm>=0.0&&pe_norm<=1.0+TE_EPSILON);
    te_ts_free(ts);
    return 0;
}

static int test_symbolic_te(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=sin(0.1*(double)i)+0.1*((double)rand()/(double)RAND_MAX-0.5);dy[i]=sin(0.1*(double)i+0.5)+0.1*((double)rand()/(double)RAND_MAX-0.5);}
    TETimeSeries *tx=te_ts_create(dx,200,"x"),*ty=te_ts_create(dy,200,"y");
    TEResult r=te_symbolic_compute(tx,ty,4,2);
    assert(r.te_xy>=0.0);
    assert(r.te_yx>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_weighted_pe(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.1*(double)i)+0.05*((double)rand()/(double)RAND_MAX-0.5);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    double wpe=te_symbolic_weighted_permutation_entropy(ts,4,2);
    assert(wpe>=0.0);
    te_ts_free(ts);
    return 0;
}

static int test_multiscale_pe(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.1*(double)i);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    double pe_scales[5];
    te_symbolic_multiscale_pe(ts,3,2,5,pe_scales);
    te_ts_free(ts);
    return 0;
}

static int test_transition_matrix(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.1*(double)i);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    int np=1; for(int i=2;i<=3;i++) np*=i;
    double *tm=calloc(np*np,sizeof(double));
    int rc=te_symbolic_transition_matrix(ts,3,2,tm);
    assert(rc==0);
    free(tm);
    te_ts_free(ts);
    return 0;
}

static int test_aape(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.1*(double)i);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    double aape=te_symbolic_aape(ts,3,2,0.5);
    assert(aape>=0.0);
    te_ts_free(ts);
    return 0;
}

static int test_conditional_pe(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.1*(double)i)+0.05*((double)rand()/(double)RAND_MAX-0.5);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    double cpe=te_symbolic_conditional_pe(ts,3,2);
    assert(cpe>=0.0);
    te_ts_free(ts);
    return 0;
}

static int test_symbolic_bias_corrected(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=sin(0.1*(double)i);dy[i]=sin(0.1*(double)i+1.0);}
    TETimeSeries *tx=te_ts_create(dx,200,"x"),*ty=te_ts_create(dy,200,"y");
    double tec=te_symbolic_bias_corrected(tx,ty,4,2,20);
    assert(tec>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

int main(void) {
    printf("=== test_te_symbolic ===\n");
    test_permutation_encode(); printf("  PASS permutation_encode\n");
    test_pattern_count(); printf("  PASS pattern_count\n");
    test_permutation_entropy(); printf("  PASS permutation_entropy\n");
    test_symbolic_te(); printf("  PASS symbolic_te\n");
    test_weighted_pe(); printf("  PASS weighted_pe\n");
    test_multiscale_pe(); printf("  PASS multiscale_pe\n");
    test_transition_matrix(); printf("  PASS transition_matrix\n");
    test_aape(); printf("  PASS AAPE\n");
    test_conditional_pe(); printf("  PASS conditional_pe\n");
    test_symbolic_bias_corrected(); printf("  PASS bias_corrected\n");
    printf("=== test_te_symbolic: ALL PASSED ===\n");
    return 0;
}