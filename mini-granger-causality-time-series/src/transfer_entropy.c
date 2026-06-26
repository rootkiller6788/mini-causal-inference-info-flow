#include "transfer_entropy.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
TransferEntropy* te_create(void){return calloc(1,sizeof(TransferEntropy));}
void te_free(TransferEntropy*te){free(te);}
double entropy_histogram(const double*data,int n,int n_bins){
    double dmin=data[0],dmax=data[0];
    for(int i=1;i<n;i++){if(data[i]<dmin)dmin=data[i];if(data[i]>dmax)dmax=data[i];}
    double bw=(dmax-dmin)/n_bins; if(bw<1e-10)bw=1e-10;
    int*counts=calloc((size_t)n_bins,sizeof(int));
    for(int i=0;i<n;i++){int b=(int)((data[i]-dmin)/bw);if(b<0)b=0;if(b>=n_bins)b=n_bins-1;counts[b]++;}
    double h=0; for(int i=0;i<n_bins;i++){if(counts[i]>0){double p=(double)counts[i]/n;h-=p*log(p);}}
    free(counts); return h;
}
double mutual_information_histogram(const double*x,const double*y,int n,int n_bins){
    double hx=entropy_histogram(x,n,n_bins);
    double hy=entropy_histogram(y,n,n_bins);
    double*xy=malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++)xy[i]=(x[i]+y[i])/2.0;
    double hxy=entropy_histogram(xy,n,n_bins);
    free(xy); return hx+hy-hxy;
}
int te_compute(TransferEntropy*te,const double*x,const double*y,int length,int n_bins,int k_history){
    if(!te||!x||!y||length<k_history+2)return-1;
    int n=length-k_history;
    double*hx=malloc((size_t)n*sizeof(double));
    double*hy=malloc((size_t)n*sizeof(double));
    double*hy_future=malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++){
        hy[i]=y[i+k_history-1];
        hy_future[i]=y[i+k_history];
        double sx=0; for(int k=0;k<k_history;k++)sx+=x[i+k];
        hx[i]=sx/k_history;
    }
    te->hx=entropy_histogram(hx,n,n_bins);
    te->hy=entropy_histogram(hy,n,n_bins);
    te->mi_xy=mutual_information_histogram(hx,hy,n,n_bins);
    double*h_both=malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++)h_both[i]=(hy_future[i]+hy[i])/2.0;
    double h_joint=entropy_histogram(h_both,n,n_bins);
    te->te_xy=entropy_histogram(hy_future,n,n_bins)-(h_joint-entropy_histogram(hy,n,n_bins));
    if(te->te_xy<0)te->te_xy=0;
    te->n_bins=n_bins; te->n_samples=n;
    free(hx);free(hy);free(hy_future);free(h_both);
    return 0;
}
int te_significance(TransferEntropy*te,const double*x,const double*y,int length,int n_bins,int k_history,int n_surrogates,double alpha){
    if(!te)return-1;
    te_compute(te,x,y,length,n_bins,k_history);
    double orig=te->te_xy; int count=0;
    double*surr=malloc((size_t)length*sizeof(double));
    for(int s=0;s<n_surrogates;s++){
        for(int i=0;i<length;i++)surr[i]=y[(i+length/3)%length];
        TransferEntropy tes; tes.te_xy=0;
        te_compute(&tes,x,surr,length,n_bins,k_history);
        if(tes.te_xy>=orig)count++;
    }
    te->p_value=(double)(count+1)/(n_surrogates+1);
    te->is_significant=(te->p_value<alpha);
    free(surr); return 0;
}
void te_print(const TransferEntropy*te){if(!te)return;printf("TE: xy=%.4f yx=%.4f sig=%s p=%.4f\n",te->te_xy,te->te_yx,te->is_significant?"YES":"NO",te->p_value);}
