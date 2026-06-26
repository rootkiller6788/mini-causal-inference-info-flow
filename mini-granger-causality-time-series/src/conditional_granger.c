#include "conditional_granger.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
ConditionalGranger* cg_create(void){return calloc(1,sizeof(ConditionalGranger));}
void cg_free(ConditionalGranger*cg){free(cg);}
int cg_test(ConditionalGranger*cg,const double*x,const double*y,const double*z,int length,int max_lag,double alpha){
    if(!cg||!x||!y||!z||length<max_lag+2)return-1;
    int p=max_lag,n=length-p;
    int k_with=1+3*p,k_without=1+2*p;
    Matrix*Xw=mat_create(n,k_with);
    Matrix*Xwo=mat_create(n,k_without);
    double*Y=malloc((size_t)n*sizeof(double));
    double*beta=malloc((size_t)k_with*sizeof(double));
    double*resid_w=malloc((size_t)n*sizeof(double));
    double*resid_wo=malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++){
        Xw->data[i*k_with+0]=1.0; Xwo->data[i*k_without+0]=1.0;
        for(int lag=0;lag<p;lag++){
            Xw->data[i*k_with+1+lag]=y[i+p-1-lag];
            Xw->data[i*k_with+1+p+lag]=x[i+p-1-lag];
            Xw->data[i*k_with+1+2*p+lag]=z[i+p-1-lag];
            Xwo->data[i*k_without+1+lag]=y[i+p-1-lag];
            Xwo->data[i*k_without+1+p+lag]=z[i+p-1-lag];
        }
        Y[i]=y[i+p];
    }
    double rss_w=ols_estimate(Xw,Y,beta,resid_w,n,k_with);
    double rss_wo=ols_estimate(Xwo,Y,beta,resid_wo,n,k_without);
    cg->rss_with_x=rss_w; cg->rss_without_x=rss_wo;
    cg->lag_order=p; cg->n_obs=n; cg->n_conditioning=1;
    cg->significance_level=alpha;
    double f=((rss_wo-rss_w)/p)/(rss_w/(n-k_with));
    if(f<0)f=0;
    cg->f_statistic=f;
    double xv=(double)(n-k_with)/(n-k_with+p*f);
    double a=(n-k_with)/2.0,b=p/2.0;
    double ib=pow(xv,a)*pow(1-xv,b-1);
    cg->p_value=1.0-ib;
    cg->is_significant=(cg->p_value<alpha);
    mat_free(Xw);mat_free(Xwo);free(Y);free(beta);free(resid_w);free(resid_wo);
    return 0;
}
int cg_test_multivariate(ConditionalGranger*cg,const double*x,const double*y,const TimeSeries*Z,int length,int max_lag,double alpha){
    if(!Z||Z->dim<1)return cg_test(cg,x,y,NULL,length,max_lag,alpha);
    return cg_test(cg,x,y,Z->values,length,max_lag,alpha);
}
void cg_print(const ConditionalGranger*cg){if(!cg)return;printf("CondGranger: F=%.4f p=%.4f sig=%s\n",cg->f_statistic,cg->p_value,cg->is_significant?"YES":"NO");}
