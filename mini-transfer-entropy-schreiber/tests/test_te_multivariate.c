/* test_te_multivariate.c -- Tests for multivariate/conditional TE */
#include "te_multivariate.h"
#include "te_schreiber.h"
#include "te_effective.h"
#include "te_core.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

static int test_conditional_te(void) {
    double dx[300],dy[300],dz[300];
    for(int i=0;i<300;i++){dx[i]=(double)(i%8);dy[i]=dx[(i+2)%300];dz[i]=(double)(i%4);}
    TETimeSeries *tx=te_ts_create(dx,300,"x"),*ty=te_ts_create(dy,300,"y"),*tz=te_ts_create(dz,300,"z");
    TEResult r=te_conditional_transfer_entropy(tx,ty,tz,2,2,2);
    assert(r.te_xy>=0.0);
    double pte=te_partial_transfer_entropy(tx,ty,tz,2,2,2);
    assert(pte>=0.0);
    te_ts_free(tx); te_ts_free(ty); te_ts_free(tz);
    return 0;
}

static int test_causal_network(void) {
    double d0[200],d1[200],d2[200];
    for(int i=0;i<200;i++){d0[i]=sin(0.1*(double)i);d1[i]=d0[i]*0.5+0.5*((double)rand()/(double)RAND_MAX);d2[i]=cos(0.1*(double)i);}
    TETimeSeries *ts[3];
    ts[0]=te_ts_create(d0,200,"v0"); ts[1]=te_ts_create(d1,200,"v1"); ts[2]=te_ts_create(d2,200,"v2");
    double *mat[3];
    for(int i=0;i<3;i++){mat[i]=calloc(3,sizeof(double));}
    te_causal_network(3,ts,2,2,mat);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) assert(mat[i][j]>=0.0);
    double density,reciprocity,transitivity;
    double *flat=calloc(9,sizeof(double));
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)flat[i*3+j]=mat[i][j];
    te_network_density(flat,3,&density,&reciprocity,&transitivity);
    assert(density>=0.0&&density<=1.0);
    free(flat);
    for(int i=0;i<3;i++) free(mat[i]);
    te_ts_free(ts[0]); te_ts_free(ts[1]); te_ts_free(ts[2]);
    return 0;
}

static int test_pid(void) {
    double dt[200],s1[200],s2[200];
    for(int i=0;i<200;i++){dt[i]=sin(0.1*(double)i);s1[i]=dt[i]*0.6+0.4*((double)rand()/(double)RAND_MAX);s2[i]=dt[i]*0.4+0.6*((double)rand()/(double)RAND_MAX);}
    TETimeSeries *tt=te_ts_create(dt,200,"target"),*ts1=te_ts_create(s1,200,"s1"),*ts2=te_ts_create(s2,200,"s2");
    double red,u1,u2,syn;
    te_partial_information_decomposition(tt,ts1,ts2,2,&red,&u1,&u2,&syn);
    assert(red>=0.0); assert(u1>=0.0); assert(u2>=0.0); assert(syn>=0.0);
    te_ts_free(tt); te_ts_free(ts1); te_ts_free(ts2);
    return 0;
}

static int test_pagerank(void) {
    double mat[9]={0,0.3,0.1,0.1,0,0.2,0.1,0.2,0};
    double ranks[3];
    te_pagerank(mat,3,0.85,ranks,50);
    double sum=0; for(int i=0;i<3;i++){assert(ranks[i]>=0.0);sum+=ranks[i];}
    assert(fabs(sum-1.0)<0.01);
    return 0;
}

static int test_label_propagation(void) {
    double mat[16]={0};
    for(int i=0;i<4;i++) for(int j=0;j<4;j++) if(i!=j) mat[i*4+j]=0.1;
    int labels[4];
    int nc=te_label_propagation(mat,4,labels,20);
    assert(nc>=1);
    return 0;
}

static int test_granger_f(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=sin(0.1*(double)i);dy[i]=cos(0.1*(double)i);}
    TETimeSeries *tx=te_ts_create(dx,200,"x"),*ty=te_ts_create(dy,200,"y");
    double F=te_granger_f_test(tx,ty,3);
    assert(F>=0.0);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

static int test_network_degrees(void) {
    int adj[9]={0,1,0,0,0,1,1,0,0};
    int in[3],out[3];
    te_network_all_degrees(adj,3,in,out);
    assert(in[0]==1); assert(out[0]==1);
    return 0;
}

static int test_optimal_delay(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.2*(double)i);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    int d=te_estimate_optimal_delay(ts,20,8);
    assert(d>=1&&d<=20);
    te_ts_free(ts);
    return 0;
}

static int test_fnn(void) {
    double data[200];
    for(int i=0;i<200;i++) data[i]=sin(0.1*(double)i)+0.05*((double)rand()/(double)RAND_MAX-0.5);
    TETimeSeries *ts=te_ts_create(data,200,"sin");
    int dim=te_estimate_embedding_dimension(ts,6,5,0.5,0.1);
    assert(dim>=1&&dim<=6);
    te_ts_free(ts);
    return 0;
}

static int test_embedding_spectrum(void) {
    double dx[200],dy[200];
    for(int i=0;i<200;i++){dx[i]=(double)(i%8);dy[i]=dx[(i+1)%200];}
    TETimeSeries *tx=te_ts_create(dx,200,"x"),*ty=te_ts_create(dy,200,"y");
    double tev[5];
    te_embedding_spectrum(tx,ty,5,tev);
    te_ts_free(tx); te_ts_free(ty);
    return 0;
}

int main(void) {
    printf("=== test_te_multivariate ===\n");
    test_conditional_te(); printf("  PASS conditional_te\n");
    test_causal_network(); printf("  PASS causal_network\n");
    test_pid(); printf("  PASS PID\n");
    test_pagerank(); printf("  PASS PageRank\n");
    test_label_propagation(); printf("  PASS label_propagation\n");
    test_granger_f(); printf("  PASS Granger_F\n");
    test_network_degrees(); printf("  PASS network_degrees\n");
    test_optimal_delay(); printf("  PASS optimal_delay\n");
    test_fnn(); printf("  PASS FNN\n");
    test_embedding_spectrum(); printf("  PASS embedding_spectrum\n");
    printf("=== test_te_multivariate: ALL PASSED ===\n");
    return 0;
}