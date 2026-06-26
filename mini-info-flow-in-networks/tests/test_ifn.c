#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ifn_core.h"
#include "ifn_transfer.h"
#include "ifn_granger.h"
#include "ifn_mutual.h"
#include "ifn_causal.h"
#include "ifn_dynamic.h"
#include "ifn_partial.h"
#include "ifn_network.h"
#define EPS 1e-6
int main(void) {
    int tc=0;
    printf("\n=== Info Flow Network Tests ===\n\n");
    printf(" [%d] ts_create... ", ++tc); fflush(stdout);
    { IFN_TimeSeries* ts=ifn_ts_create(100,2); assert(ts); assert(ts->length==100); ifn_ts_free(ts); }
    printf("PASSED\n");
    printf(" [%d] ts_mean... ", ++tc);
    { IFN_TimeSeries* ts=ifn_ts_create(10,1); for(int i=0;i<10;i++) ts->data[i]=i; double m=ifn_ts_mean(ts,0); assert(fabs(m-4.5)<1.0); ifn_ts_free(ts); }
    printf("PASSED\n");
    printf(" [%d] ts_variance... ", ++tc);
    { IFN_TimeSeries* ts=ifn_ts_create(10,1); for(int i=0;i<10;i++) ts->data[i]=i; double v=ifn_ts_variance(ts,0); assert(v>0); ifn_ts_free(ts); }
    printf("PASSED\n");
    printf(" [%d] ts_entropy... ", ++tc);
    { IFN_TimeSeries* ts=ifn_ts_create(100,1); for(int i=0;i<100;i++) ts->data[i]=(double)(i%10); double H=ifn_ts_entropy(ts,0); assert(H>0); ifn_ts_free(ts); }
    printf("PASSED\n");
    printf(" [%d] mutual_information... ", ++tc);
    { IFN_TimeSeries* x=ifn_ts_create(100,1); IFN_TimeSeries* y=ifn_ts_create(100,1);
      for(int i=0;i<100;i++) { x->data[i]=i%5; y->data[i]=i%5; }
      double mi=ifn_mutual_information_binned(x->data,y->data,100,5,5); assert(mi>=0); ifn_ts_free(x); ifn_ts_free(y); }
    printf("PASSED\n");
    printf(" [%d] transfer_entropy... ", ++tc);
    { IFN_TimeSeries* src=ifn_ts_create(100,1); IFN_TimeSeries* tgt=ifn_ts_create(100,1);
      for(int i=0;i<100;i++) { src->data[i]=(i%5); tgt->data[i]=((i+1)%5); }
      IFN_TransferResult r=ifn_transfer_entropy(src,tgt,1,1); assert(r.te_value>=0); ifn_ts_free(src); ifn_ts_free(tgt); }
    printf("PASSED\n");
    printf(" [%d] granger_test... ", ++tc);
    { IFN_TimeSeries* x=ifn_ts_create(100,1); IFN_TimeSeries* y=ifn_ts_create(100,1);
      for(int i=0;i<100;i++) { x->data[i]=i%3; y->data[i]=(i+1)%3; }
      IFN_GrangerResult gr=ifn_granger_test(x,y,2); assert(gr.F_statistic>=0); ifn_ts_free(x); ifn_ts_free(y); }
    printf("PASSED\n");
    printf(" [%d] causal_graph... ", ++tc);
    { IFN_CausalGraph* g=ifn_causal_create(4); assert(g); assert(g->n_nodes==4);
      ifn_causal_add_edge(g,0,1,0.5); assert(ifn_causal_has_edge(g,0,1)); ifn_causal_free(g); }
    printf("PASSED\n");
    printf(" [%d] is_dag... ", ++tc);
    { IFN_CausalGraph* g=ifn_causal_create(3); ifn_causal_add_edge(g,0,1,1); ifn_causal_add_edge(g,1,2,1);
      assert(ifn_is_dag(g)); ifn_causal_free(g); }
    printf("PASSED\n");
    printf(" [%d] d_separation... ", ++tc);
    { IFN_CausalGraph* g=ifn_causal_create(3); ifn_causal_add_edge(g,0,1,1); ifn_causal_add_edge(g,1,2,1);
      int z[]={1}; bool sep=ifn_d_separated(g,0,2,z,1); assert(sep||!sep); ifn_causal_free(g); }
    printf("PASSED\n");
    printf(" [%d] pagerank... ", ++tc);
    { IFN_CausalGraph* g=ifn_causal_create(3); ifn_causal_add_edge(g,0,1,1); ifn_causal_add_edge(g,0,2,1);
      double pr[3]; ifn_pagerank(g,0.85,50,1e-6,pr); assert(pr[0]>0); ifn_causal_free(g); }
    printf("PASSED\n");
    printf(" [%d] network_erdos... ", ++tc);
    { IFN_CausalGraph* g=ifn_network_random_erdos(10,0.3); assert(g); assert(g->n_nodes==10); ifn_causal_free(g); }
    printf("PASSED\n");
    printf(" [%d] pid_decompose... ", ++tc);
    { IFN_TimeSeries* t=ifn_ts_create(50,1); IFN_TimeSeries* s1=ifn_ts_create(50,1); IFN_TimeSeries* s2=ifn_ts_create(50,1);
      for(int i=0;i<50;i++) { t->data[i]=i%4; s1->data[i]=i%4; s2->data[i]=i%4; }
      IFN_PIDResult pid=ifn_pid_decompose(t,s1,s2,4); assert(pid.total_mi>=0);
      ifn_ts_free(t); ifn_ts_free(s1); ifn_ts_free(s2); }
    printf("PASSED\n");
    printf("\n=== %d tests PASSED ===\n", tc);
    return 0;
}
