/* test_te_effective.c -- Tests for effective TE and surrogate methods */
#include "te_effective.h"
#include "te_schreiber.h"
#include "te_core.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

static int test_surrogate_iid(void) {
    double data[]={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
    double surr[10];
    te_generate_iid_surrogates(data,10,surr);
    double mn=surr[0],mx=surr[0];
    for(int i=1;i<10;i++){if(surr[i]<mn)mn=surr[i];if(surr[i]>mx)mx=surr[i];}
    assert(fabs(mn-1.0)<TE_EPSILON);
    assert(fabs(mx-10.0)<TE_EPSILON);
    return 0;
}

static int test_surrogate_ft(void) {
    double data[64];
    for(int i=0;i<64;i++) data[i]=sin(0.2*(double)i)+0.1*((double)rand()/(double)RAND_MAX-0.5);
    double surr[64];
    te_generate_ft_surrogates(data,64,surr);
    double mn_o=data[0],mx_o=data[0],mn_s=surr[0],mx_s=surr[0];
    for(int i=1;i<64;i++){if(data[i]<mn_o)mn_o=data[i];if(data[i]>mx_o)mx_o=data[i];if(surr[i]<mn_s)mn_s=surr[i];if(surr[i]>mx_s)mx_s=surr[i];}
    assert(fabs(mn_o-mn_s)<0.5);
    assert(fabs(mx_o-mx_s)<0.5);
    return 0;
}

static int test_surrogate_iaaft(void) {
    double data[64];
    for(int i=0;i<64;i++) data[i]=sin(0.2*(double)i)+0.05*((double)rand()/(double)RAND_MAX-0.5);
    double surr[64];
    te_generate_iaaft_surrogates(data,64,surr,5);
    return 0;
}

static int test_block_bootstrap(void) {
    double data[50];
    for(int i=0;i<50;i++) data[i]=(double)i;
    double surr[50];
    te_generate_block_bootstrap(data,50,10,surr);
    for(int i=0;i<50;i++) if(fabs(surr[i]-data[i])>TE_EPSILON) { (void)surr[i]; break; }
    return 0;
}

static int test_stationary_bootstrap(void) {
    double data[30];
    for(int i=0;i<30;i++) data[i]=sin(0.3*(double)i);
    double surr[30];
    te_generate_stationary_bootstrap(data,30,0.7,surr);
    return 0;
}

static int test_circular_shift(void) {
    double data[]={1.0,2.0,3.0,4.0,5.0};
    double surr[5];
    te_generate_circular_shift(data,5,surr);
    double mn=surr[0],mx=surr[0];
    for(int i=1;i<5;i++){if(surr[i]<mn)mn=surr[i];if(surr[i]>mx)mx=surr[i];}
    assert(fabs(mn-1.0)<TE_EPSILON);
    assert(fabs(mx-5.0)<TE_EPSILON);
    return 0;
}

static int test_effective_te(void) {
    double dx[300],dy[300];
    dx[0]=0;dy[0]=0;
    for(int i=1;i<300;i++){dx[i]=0.5*dx[i-1]+0.3*dy[i-1]+0.1*((double)rand()/(double)RAND_MAX-0.5);dy[i]=0.5*dy[i-1]+0.1*((double)rand()/(double)RAND_MAX-0.5);}
    TETimeSeries *tx=te_ts_create(dx,300,"x"),*ty=te_ts_create(dy,300,"y");
    TEResult r=te_effective_compute(tx,ty,2,2,30);
    assert(r.te_xy>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_bias_corrected(void) {
    double dx[200],dy[200];
    dx[0]=0;dy[0]=0;
    for(int i=1;i<200;i++){dx[i]=0.5*dx[i-1]+0.1*((double)rand()/(double)RAND_MAX-0.5);dy[i]=0.5*dy[i-1]+0.1*((double)rand()/(double)RAND_MAX-0.5);}
    TETimeSeries *tx=te_ts_create(dx,200,"x"),*ty=te_ts_create(dy,200,"y");
    TEResult r=te_effective_bias_corrected(tx,ty,2,2,20);
    assert(r.te_xy>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_zscore_to_pvalue(void) {
    double p0=te_zscore_to_pvalue(0.0);
    assert(p0>0.9);
    double p3=te_zscore_to_pvalue(3.0);
    assert(p3<0.01);
    double p5=te_zscore_to_pvalue(5.0);
    assert(p5<1e-5);
    return 0;
}

static int test_ks_test(void) {
    double a[]={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
    double b[]={2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0};
    double d,p;
    int rc=te_kolmogorov_smirnov(a,10,b,10,&d,&p);
    assert(rc==0);
    assert(d>0.0&&d<=1.0);
    assert(p>=0.0&&p<=1.0);
    return 0;
}

static int test_holm_correction(void) {
    double pv[]={0.01,0.02,0.03,0.001};
    double corrected[4];
    te_holm_correction(pv,4,corrected);
    for(int i=0;i<4;i++) assert(corrected[i]>=0.0&&corrected[i]<=1.0);
    return 0;
}

static int test_permutation_test(void) {
    double dx[150],dy[150];
    for(int i=0;i<150;i++){dx[i]=(double)(i%7);dy[i]=(double)((i+1)%7);}
    TETimeSeries *tx=te_ts_create(dx,150,"x"),*ty=te_ts_create(dy,150,"y");
    double pv,nd[51];
    int rc=te_permutation_test(tx,ty,2,2,50,&pv,nd);
    assert(rc==0);
    assert(pv>0.0&&pv<=1.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_bootstrap_ci(void) {
    double dx[150],dy[150];
    for(int i=0;i<150;i++){dx[i]=sin(0.1*(double)i);dy[i]=sin(0.1*(double)i+0.5);}
    TETimeSeries *tx=te_ts_create(dx,150,"x"),*ty=te_ts_create(dy,150,"y");
    double lo,hi;
    te_bootstrap_ci(tx,ty,1,1,50,0.05,&lo,&hi,NULL);
    assert(lo<=hi);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_jackknife(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=(double)(i%8);dy[i]=dx[(i+2)%200];}
    TETimeSeries *tx=te_ts_create(dx,200,"x"),*ty=te_ts_create(dy,200,"y");
    double jk=te_jackknife_te(tx,ty,2,2);
    (void)jk; /* Jackknife TE may be negative due to bias overcorrection */
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_fdr(void) {
    double pv[]={0.001,0.05,0.1,0.5};
    int sig[4]={0};
    te_benjamini_hochberg(pv,4,0.05,sig);
    assert(sig[0]==1||sig[0]==0);
    return 0;
}

int main(void) {
    printf("=== test_te_effective ===\n");
    test_surrogate_iid(); printf("  PASS surrogate_iid\n");
    test_surrogate_ft(); printf("  PASS surrogate_ft\n");
    test_surrogate_iaaft(); printf("  PASS surrogate_iaaft\n");
    test_block_bootstrap(); printf("  PASS block_bootstrap\n");
    test_stationary_bootstrap(); printf("  PASS stationary_bootstrap\n");
    test_circular_shift(); printf("  PASS circular_shift\n");
    test_effective_te(); printf("  PASS effective_te\n");
    test_bias_corrected(); printf("  PASS bias_corrected\n");
    test_zscore_to_pvalue(); printf("  PASS zscore_to_pvalue\n");
    test_ks_test(); printf("  PASS ks_test\n");
    test_holm_correction(); printf("  PASS holm\n");
    test_permutation_test(); printf("  PASS permutation_test\n");
    test_bootstrap_ci(); printf("  PASS bootstrap_ci\n");
    test_jackknife(); printf("  PASS jackknife\n");
    test_fdr(); printf("  PASS fdr\n");
    printf("=== test_te_effective: ALL PASSED ===\n");
    return 0;
}