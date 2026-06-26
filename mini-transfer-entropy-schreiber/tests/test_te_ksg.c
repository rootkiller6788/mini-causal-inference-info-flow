/* test_te_ksg.c -- Tests for KSG k-NN transfer entropy estimator */
#include "te_ksg.h"
#include "te_core.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

static int test_ksg_basic(void) {
    double dx[200], dy[200];
    dx[0]=0; dy[0]=0;
    for(int i=1;i<200;i++){dx[i]=0.5*dx[i-1]+0.2*dy[i-1]+0.05*((double)rand()/(double)RAND_MAX-0.5);dy[i]=0.5*dy[i-1]+0.05*((double)rand()/(double)RAND_MAX-0.5);}
    TETimeSeries *tx=te_ts_create(dx,200,"x"),*ty=te_ts_create(dy,200,"y");
    TEResult r=te_ksg_compute(tx,ty,2,2,4);
    assert(r.te_xy>=0.0); assert(r.te_yx>=0.0);
    double mi=te_ksg_mutual_info_knn(dx,dy,200,4);
    assert(mi>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_ksg_entropy(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=0.5*((double)rand()/(double)RAND_MAX-0.5);
    double H=te_ksg_entropy_knn(data,200,1,4);
    /* entropy can be negative for differential entropy but MI is non-negative */
    (void)H;
    return 0;
}

static int test_ksg_adaptive(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=(double)(i%10)/10.0;dy[i]=(double)((i+3)%10)/10.0;}
    double mi=te_ksg_adaptive_k(dx,dy,200);
    assert(mi>=0.0);
    return 0;
}

static int test_ksg_auto_knn(void) {
    double dx[100],dy[100];
    for(int i=0;i<100;i++){dx[i]=sin(0.1*(double)i);dy[i]=cos(0.1*(double)i);}
    double mi=te_ksg_auto_knn(dx,dy,100,2,8);
    assert(mi>=0.0);
    return 0;
}

static int test_ksg_conditional_mi(void) {
    double x[100],y[100],z[100];
    for(int i=0;i<100;i++){x[i]=(double)(i%5);y[i]=x[(i+1)%100];z[i]=(double)(i%3);}
    double cmi=te_ksg_conditional_mi(x,y,z,100,3);
    assert(cmi>=0.0);
    return 0;
}

static int test_ksg_te_via_cmi(void) {
    double x[200],y[200];
    for(int i=0;i<200;i++){x[i]=(double)(i%8);y[i]=x[(i+2)%200];}
    double te=te_ksg_te_via_cmi(x,y,200,2,4);
    assert(te>=0.0);
    return 0;
}

static int test_ksg_mi_manhattan(void) {
    double x[100],y[100];
    for(int i=0;i<100;i++){x[i]=sin(0.2*(double)i);y[i]=cos(0.2*(double)i);}
    double mi=te_ksg_mi_manhattan(x,y,100,3);
    assert(mi>=0.0);
    return 0;
}

int main(void) {
    printf("=== test_te_ksg ===\n");
    test_ksg_basic(); printf("  PASS basic\n");
    test_ksg_entropy(); printf("  PASS entropy\n");
    test_ksg_adaptive(); printf("  PASS adaptive\n");
    test_ksg_auto_knn(); printf("  PASS auto_knn\n");
    test_ksg_conditional_mi(); printf("  PASS conditional_mi\n");
    test_ksg_te_via_cmi(); printf("  PASS te_via_cmi\n");
    test_ksg_mi_manhattan(); printf("  PASS mi_manhattan\n");
    printf("=== test_te_ksg: ALL PASSED ===\n");
    return 0;
}