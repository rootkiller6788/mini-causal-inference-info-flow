#include "granger_test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

GrangerTestResult* gt_create(void){return calloc(1,sizeof(GrangerTestResult));}
void gt_free(GrangerTestResult*gt){free(gt);}

/* Regularized incomplete beta function I_x(a,b) via continued fraction */
static double beta_cf(double a, double b, double x) {
    int max_iter = 200;
    double eps = 1e-12;
    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (fabs(d) < 1e-30) d = 1e-30;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= max_iter; m++) {
        int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < eps) break;
    }
    return h;
}

static double beta_reg(double a, double b, double x) {
    if (x < 0.0 || x > 1.0) return 0.0;
    if (x == 0.0 || x == 1.0) return x;
    double bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt * beta_cf(a, b, x) / a;
    else
        return 1.0 - bt * beta_cf(b, a, 1.0 - x) / b;
}

double gt_f_distribution_pvalue(double f, int df1, int df2) {
    if (f <= 0.0) return 1.0;
    double x = df2 / (df2 + df1 * f);
    return beta_reg(df2 / 2.0, df1 / 2.0, x);
}

/* Gamma function via Stirling approximation */
#ifndef lgamma
double lgamma(double xx) {
    double x = xx, y = xx, tmp = x + 5.2421875;
    tmp = (x + 0.5) * log(tmp) - tmp;
    double ser = 0.999999999999997092;
    double cof[] = {76.18009172947146, -86.50532032941678, 24.01409824083091,
                    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
    for (int j = 0; j < 6; j++) { y += 1.0; ser += cof[j] / y; }
    return tmp + log(2.5066282746310005 * ser / x);
}
#endif

int gt_test(GrangerTestResult*gt,const double*x,const double*y,int length,int max_lag,double alpha){
    if(!gt||!x||!y||length<max_lag+2)return-1;
    int p=gt_optimal_lag(x,y,length,max_lag);
    int n=length-p,k_ur=1+2*p,k_r=1+p;
    Matrix*X_ur=mat_create(n,k_ur);
    Matrix*X_r=mat_create(n,k_r);
    double*Y=malloc((size_t)n*sizeof(double));
    double*beta=malloc((size_t)k_ur*sizeof(double));
    double*resid=malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++){
        X_ur->data[i*k_ur+0]=1.0; X_r->data[i*k_r+0]=1.0;
        for(int lag=0;lag<p;lag++){
            X_ur->data[i*k_ur+1+lag]=y[i+p-1-lag];
            X_r->data[i*k_r+1+lag]=y[i+p-1-lag];
            X_ur->data[i*k_ur+1+p+lag]=x[i+p-1-lag];
        }
        Y[i]=y[i+p];
    }
    double rss_r=ols_estimate(X_r,Y,beta,resid,n,k_r);
    double rss_ur=ols_estimate(X_ur,Y,beta,resid,n,k_ur);
    gt->rss_unrestricted=rss_ur; gt->rss_restricted=rss_r;
    gt->lag_order=p; gt->n_obs=n; gt->significance_level=alpha;
    if(rss_ur<1e-15)rss_ur=1e-15;
    double f=((rss_r-rss_ur)/p)/(rss_ur/(n-k_ur));
    if(f<0)f=0;
    gt->f_statistic=f;
    gt->p_value=gt_f_distribution_pvalue(f,p,n-k_ur);
    gt->is_significant=(gt->p_value<alpha);
    mat_free(X_ur);mat_free(X_r);free(Y);free(beta);free(resid);
    return 0;
}

int gt_test_bidirectional(GrangerTestResult*gt_xy,GrangerTestResult*gt_yx,const double*x,const double*y,int length,int max_lag,double alpha){
    gt_test(gt_xy,x,y,length,max_lag,alpha);
    gt_test(gt_yx,y,x,length,max_lag,alpha);
    return 0;
}

int gt_optimal_lag(const double*x,const double*y,int length,int max_lag){
    double best_aic=1e10; int best_p=1;
    for(int p=1;p<=max_lag;p++){
        int n=length-p,k=1+2*p;
        Matrix*X=mat_create(n,k); double*Y=malloc((size_t)n*sizeof(double));
        double*beta=malloc((size_t)k*sizeof(double));
        double*resid=malloc((size_t)n*sizeof(double));
        for(int i=0;i<n;i++){
            X->data[i*k+0]=1.0;
            for(int lag=0;lag<p;lag++){X->data[i*k+1+lag]=y[i+p-1-lag];X->data[i*k+1+p+lag]=x[i+p-1-lag];}
            Y[i]=y[i+p];
        }
        double rss=ols_estimate(X,Y,beta,resid,n,k);
        double aic=log(rss/n)+2.0*k/n;
        if(aic<best_aic){best_aic=aic;best_p=p;}
        mat_free(X);free(Y);free(beta);free(resid);
    }
    return best_p;
}

void gt_print(const GrangerTestResult*gt){if(!gt)return;printf("Granger: F=%.4f p=%.6f sig=%s (lag=%d)\n",gt->f_statistic,gt->p_value,gt->is_significant?"YES":"NO",gt->lag_order);}
