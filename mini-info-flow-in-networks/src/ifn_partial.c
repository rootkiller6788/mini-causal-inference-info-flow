#include "ifn_partial.h"
#include "ifn_transfer.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

IFN_PIDResult ifn_pid_decompose(const IFN_TimeSeries* target,
    const IFN_TimeSeries* source1, const IFN_TimeSeries* source2, int n_bins) {
    IFN_PIDResult pid; memset(&pid,0,sizeof(pid));
    if (!target||!source1||!source2) return pid;
    int n=target->length;
    double I_s1=ifn_mutual_information_binned(target->data,source1->data,n,n_bins,n_bins);
    double I_s2=ifn_mutual_information_binned(target->data,source2->data,n,n_bins,n_bins);
    double* s12=malloc(n*2*sizeof(double));
    if (!s12) return pid;
    for (int i=0;i<n;i++) { s12[i*2]=source1->data[i]; s12[i*2+1]=source2->data[i]; }
    double I_s12=ifn_mutual_information_binned(target->data,s12,n,n_bins,2*n_bins);
    free(s12);
    pid.total_mi=I_s12;
    double I_s1s2=ifn_mutual_information_binned(source1->data,source2->data,n,n_bins,n_bins);
    double r_min=I_s1<I_s2?I_s1:I_s2;
    pid.redundant=r_min<I_s1s2?r_min:I_s1s2;
    if (pid.redundant<0) pid.redundant=0;
    pid.unique_y=I_s1-pid.redundant; if (pid.unique_y<0) pid.unique_y=0;
    pid.unique_z=I_s2-pid.redundant; if (pid.unique_z<0) pid.unique_z=0;
    pid.synergistic=pid.total_mi-pid.unique_y-pid.unique_z-pid.redundant;
    if (pid.synergistic<0) pid.synergistic=0;
    return pid;
}

/* Discrete PID decomposition.
 * Computes joint probability distributions from discrete state sequences
 * then decomposes I(T;S1,S2) into unique, redundant, synergistic components.
 * Uses Williams & Beer (2010) redundancy lattice framework. */
IFN_PIDResult ifn_pid_decompose_discrete(const int* target, const int* src1,
    const int* src2, int T, int n_st, int n_s1, int n_s2) {
    IFN_PIDResult pid; memset(&pid,0,sizeof(pid));
    if (!target||!src1||!src2||T<1||n_st<1||n_s1<1||n_s2<1) return pid;

    /* Joint probability P(T, S1, S2) */
    int js=n_st*n_s1*n_s2;
    double* joint=calloc(js,sizeof(double));
    if (!joint) return pid;
    for (int t=0;t<T;t++) {
        if (target[t]<n_st&&src1[t]<n_s1&&src2[t]<n_s2)
            joint[(target[t]*n_s1+src1[t])*n_s2+src2[t]]+=1.0/T;
    }
    /* Marginals */
    double* p_t=calloc(n_st,sizeof(double));
    double* p_s1=calloc(n_s1,sizeof(double));
    double* p_s2=calloc(n_s2,sizeof(double));
    for (int t=0;t<n_st;t++) for (int a=0;a<n_s1;a++) for (int b=0;b<n_s2;b++) {
        double p=joint[(t*n_s1+a)*n_s2+b];
        p_t[t]+=p; p_s1[a]+=p; p_s2[b]+=p;
    }
    /* Conditional MI: I(T;S1) and I(T;S2) */
    double I_t_s1=0.0, I_t_s2=0.0;
    for (int t=0;t<n_st;t++) for (int a=0;a<n_s1;a++) {
        double p_ta=0.0; for (int b=0;b<n_s2;b++) p_ta+=joint[(t*n_s1+a)*n_s2+b];
        if (p_ta>1e-12&&p_t[t]>1e-12&&p_s1[a]>1e-12)
            I_t_s1+=p_ta*log(p_ta/(p_t[t]*p_s1[a]));
    }
    for (int t=0;t<n_st;t++) for (int b=0;b<n_s2;b++) {
        double p_tb=0.0; for (int a=0;a<n_s1;a++) p_tb+=joint[(t*n_s1+a)*n_s2+b];
        if (p_tb>1e-12&&p_t[t]>1e-12&&p_s2[b]>1e-12)
            I_t_s2+=p_tb*log(p_tb/(p_t[t]*p_s2[b]));
    }
    /* Joint MI: I(T; S1,S2) */
    double I_joint=0.0;
    for (int t=0;t<n_st;t++) for (int a=0;a<n_s1;a++) for (int b=0;b<n_s2;b++) {
        double p_all=joint[(t*n_s1+a)*n_s2+b];
        double p_s=0.0; for (int tt=0;tt<n_st;tt++) p_s+=joint[(tt*n_s1+a)*n_s2+b];
        if (p_all>1e-12&&p_t[t]>1e-12&&p_s>1e-12)
            I_joint+=p_all*log(p_all/(p_t[t]*p_s));
    }
    pid.total_mi=I_joint>0?I_joint:0.0;
    /* Redundancy: min of I(T;S1) and I(T;S2) via MMI (Minimum Mutual Information) */
    double r_min=I_t_s1<I_t_s2?I_t_s1:I_t_s2;
    pid.redundant=r_min>0?r_min:0.0;
    pid.unique_y=I_t_s1-pid.redundant; if (pid.unique_y<0) pid.unique_y=0;
    pid.unique_z=I_t_s2-pid.redundant; if (pid.unique_z<0) pid.unique_z=0;
    pid.synergistic=pid.total_mi-pid.unique_y-pid.unique_z-pid.redundant;
    if (pid.synergistic<0) pid.synergistic=0;

    free(joint); free(p_t); free(p_s1); free(p_s2);
    return pid;
}

double ifn_unique_information(const IFN_TimeSeries* target,
    const IFN_TimeSeries* source, const IFN_TimeSeries* other, int n_bins) {
    IFN_PIDResult pid=ifn_pid_decompose(target,source,other,n_bins);
    return pid.unique_y;
}

double ifn_redundant_information(const IFN_TimeSeries* target,
    const IFN_TimeSeries* src1, const IFN_TimeSeries* src2, int n_bins) {
    IFN_PIDResult pid=ifn_pid_decompose(target,src1,src2,n_bins);
    return pid.redundant;
}

double ifn_synergistic_information(const IFN_TimeSeries* target,
    const IFN_TimeSeries* src1, const IFN_TimeSeries* src2, int n_bins) {
    IFN_PIDResult pid=ifn_pid_decompose(target,src1,src2,n_bins);
    return pid.synergistic;
}

double ifn_transfer_synergy(const IFN_TimeSeries* source1,
    const IFN_TimeSeries* source2, const IFN_TimeSeries* target, int k, int l) {
    IFN_TransferResult te1=ifn_transfer_entropy(source1,target,k,l);
    IFN_TransferResult te2=ifn_transfer_entropy(source2,target,k,l);
    int n=target->length;
    double* s12=malloc(n*2*sizeof(double));
    if (!s12) return 0.0;
    for (int i=0;i<n;i++) { s12[i*2]=source1->data[i]; s12[i*2+1]=source2->data[i]; }
    IFN_TimeSeries joint={s12,n,2,0,0,NULL,0};
    IFN_TransferResult tej=ifn_transfer_entropy(&joint,target,k,l);
    free(s12);
    double synergy=tej.te_value-te1.te_value-te2.te_value;
    return synergy>0?synergy:0.0;
}

void ifn_pid_print(const IFN_PIDResult* pid) {
    if (!pid) return;
    printf("PID: total=%.4f uniq1=%.4f uniq2=%.4f R=%.4f S=%.4f\n",
        pid->total_mi,pid->unique_y,pid->unique_z,pid->redundant,pid->synergistic);
}
