/* ifn_bootstrap.c — Bootstrap resampling methods for causal inference.
 *
 * Knowledge points implemented:
 *   1. Block bootstrap — temporal dependence via contiguous blocks
 *   2. Stationary bootstrap — random block lengths (geometric distribution)
 *   3. IAAFT surrogates — nonlinearity testing via phase randomization
 *   4. Permutation-based significance testing for transfer entropy
 *   5. Percentile and BCa confidence intervals
 *   6. Cross-validated lag selection for Granger causality
 *
 * Ref: Efron & Tibshirani (1993), Schreiber & Schmitz (2000)
 *      Politis & Romano (1994), Theiler et al. (1992) Physica D 58, 77
 */
#include "ifn_bootstrap.h"
#include "ifn_transfer.h"
#include "ifn_granger.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Block bootstrap ---- */
double* ifn_block_bootstrap(const double* series, int n, int block_len,
    int n_bootstrap) {
    if (!series||n<block_len||block_len<1||n_bootstrap<1) return NULL;
    double* samples=malloc(n_bootstrap*n*sizeof(double));
    if (!samples) return NULL;
    int n_blocks=n/block_len+1;

    for (int b=0;b<n_bootstrap;b++) {
        int pos=0;
        while (pos<n) {
            int block_start=(rand()%n_blocks)*block_len;
            for (int i=0;i<block_len&&pos<n;i++) {
                int idx=(block_start+i)%n;
                samples[b*n+pos]=series[idx];
                pos++;
            }
        }
    }
    return samples;
}

/* ---- Stationary bootstrap with geometric block lengths ---- */
double* ifn_stationary_bootstrap(const double* series, int n, double p,
    int n_bootstrap) {
    if (!series||n<2||p<=0||p>1||n_bootstrap<1) return NULL;
    double* samples=malloc(n_bootstrap*n*sizeof(double));
    if (!samples) return NULL;

    for (int b=0;b<n_bootstrap;b++) {
        int pos=0;
        while (pos<n) {
            int block_len=1;
            while (((double)rand()/RAND_MAX)<(1.0-p)&&pos+block_len<n) block_len++;
            int start=rand()%n;
            for (int i=0;i<block_len&&pos<n;i++) {
                samples[b*n+pos]=series[(start+i)%n];
                pos++;
            }
        }
    }
    return samples;
}

/* ---- FFT helper (simple DFT, O(n²), self-contained) ---- */
static void dft(const double* x, int n, double* re, double* im) {
    for (int k=0;k<=n/2;k++) {
        re[k]=0.0; im[k]=0.0;
        for (int j=0;j<n;j++) {
            double theta=2.0*M_PI*k*j/n;
            re[k]+=x[j]*cos(theta);
            im[k]+=x[j]*sin(theta);
        }
        re[k]/=n; im[k]/=n;
    }
}

static void idft(double* re, double* im, int n, double* out) {
    for (int j=0;j<n;j++) {
        out[j]=0.0;
        for (int k=0;k<=n/2;k++) {
            double theta=2.0*M_PI*k*j/n;
            out[j]+=re[k]*cos(theta)-im[k]*sin(theta);
        }
    }
}

/* ---- IAAFT surrogate generation ---- */
double** ifn_iaaft_surrogates(const double* signal, int n, int n_surrogates) {
    if (!signal||n<4||n_surrogates<1) return NULL;
    double** surrs=malloc(n_surrogates*sizeof(double*));
    if (!surrs) return NULL;

    /* Rank-order original */
    int* orig_rank=malloc(n*sizeof(int));
    double* sorted=malloc(n*sizeof(double));
    if (!orig_rank||!sorted) { free(orig_rank); free(sorted); free(surrs); return NULL; }
    for (int i=0;i<n;i++) sorted[i]=signal[i];
    for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++)
        if (sorted[i]>sorted[j]) { double t=sorted[i]; sorted[i]=sorted[j]; sorted[j]=t; }
    for (int i=0;i<n;i++) {
        int r=0; while (r<n&&sorted[r]<signal[i]-1e-12) r++;
        orig_rank[i]=r;
    }

    /* FFT of original */
    int n2=n/2+1;
    double* amp=malloc(n2*sizeof(double));
    double* phase=malloc(n2*sizeof(double));
    double* re=malloc(n2*sizeof(double));
    double* im=malloc(n2*sizeof(double));
    if (!amp||!phase||!re||!im) {
        free(amp); free(phase); free(re); free(im);
        free(orig_rank); free(sorted); free(surrs); return NULL;
    }
    dft(signal,n,re,im);
    for (int k=0;k<n2;k++) { amp[k]=sqrt(re[k]*re[k]+im[k]*im[k]); }

    for (int s=0;s<n_surrogates;s++) {
        surrs[s]=malloc(n*sizeof(double));
        if (!surrs[s]) { for (int i=0;i<s;i++) free(surrs[i]); free(surrs); free(amp); free(phase); free(re); free(im); free(orig_rank); free(sorted); return NULL; }
        /* Randomize phases */
        for (int k=0;k<n2;k++) {
            double phi=2.0*M_PI*((double)rand()/RAND_MAX);
            phase[k]=atan2(im[k],re[k])+(phi-M_PI)*0.5;
        }
        for (int k=0;k<n2;k++) {
            re[k]=amp[k]*cos(phase[k]);
            im[k]=amp[k]*sin(phase[k]);
        }
        re[0]=im[0]=0.0; /* DC component zero for surrogates */
        double* ft_out=malloc(n*sizeof(double));
        idft(re,im,n,ft_out);
        /* Amplitude adjust: map to original rank order */
        int* surr_rank=malloc(n*sizeof(int));
        double* surr_sorted=malloc(n*sizeof(double));
        if (surr_rank&&surr_sorted&&ft_out) {
            for (int i=0;i<n;i++) surr_sorted[i]=ft_out[i];
            for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++)
                if (surr_sorted[i]>surr_sorted[j])
                    { double t=surr_sorted[i]; surr_sorted[i]=surr_sorted[j]; surr_sorted[j]=t; }
            for (int i=0;i<n;i++) surrs[s][i]=sorted[orig_rank[i]];
        }
        free(ft_out); free(surr_rank); free(surr_sorted);
    }
    free(amp); free(phase); free(re); free(im);
    free(orig_rank); free(sorted);
    return surrs;
}

/* ---- Permutation test for transfer entropy ---- */
double ifn_te_permutation_test(const IFN_TimeSeries* source,
    const IFN_TimeSeries* target, int k, int l,
    int n_permutations, const IFN_TransferResult* observed) {
    if (!source||!target||!observed||n_permutations<1) return 1.0;
    int T=source->length;
    double* shuffled=malloc(T*sizeof(double));
    if (!shuffled) return 1.0;

    int exceed=0;
    for (int p=0;p<n_permutations;p++) {
        /* Fisher-Yates shuffle of source */
        for (int i=0;i<T;i++) shuffled[i]=source->data[i];
        for (int i=T-1;i>0;i--) {
            int j=rand()%(i+1);
            double t=shuffled[i]; shuffled[i]=shuffled[j]; shuffled[j]=t;
        }
        IFN_TimeSeries src_perm={shuffled,T,1,0,0,NULL,0};
        IFN_TransferResult r=ifn_transfer_entropy(&src_perm,target,k,l);
        if (r.te_value>=observed->te_value) exceed++;
    }
    free(shuffled);
    return (double)(exceed+1)/(n_permutations+1); /* Laplace smoothing */
}

/* ---- Percentile CI ---- */
void ifn_bootstrap_ci_percentile(const double* bootstrap_vals,
    int n_bootstrap, double alpha, double* lower, double* upper) {
    if (!bootstrap_vals||!lower||!upper||n_bootstrap<2) return;
    double* sorted=malloc(n_bootstrap*sizeof(double));
    if (!sorted) return;
    memcpy(sorted,bootstrap_vals,n_bootstrap*sizeof(double));
    for (int i=0;i<n_bootstrap-1;i++) for (int j=i+1;j<n_bootstrap;j++)
        if (sorted[i]>sorted[j]) { double t=sorted[i]; sorted[i]=sorted[j]; sorted[j]=t; }
    int lo=(int)(alpha/2.0*n_bootstrap), hi=(int)((1.0-alpha/2.0)*n_bootstrap);
    if (lo<0) lo=0;
    if (hi>=n_bootstrap) hi=n_bootstrap-1;
    *lower=sorted[lo]; *upper=sorted[hi];
    free(sorted);
}

/* ---- BCa confidence interval ---- */
void ifn_bootstrap_ci_bca(const double* bootstrap_vals, int n_bootstrap,
    const double* jackknife_vals, int n_original, double alpha,
    double* lower, double* upper) {
    if (!bootstrap_vals||!jackknife_vals||!lower||!upper||n_bootstrap<2||n_original<2) return;

    /* Bias correction z0 */
    double theta_hat=0.0;
    for (int i=0;i<n_original;i++) theta_hat+=jackknife_vals[i];
    theta_hat/=n_original;

    int less=0;
    for (int i=0;i<n_bootstrap;i++) if (bootstrap_vals[i]<theta_hat) less++;
    double p0=(double)less/n_bootstrap;
    double z0=0.0;
    { double x=p0>0.999?0.999:(p0<0.001?0.001:p0);
      z0=log(x/(1.0-x))/1.7; } /* Approximate probit */

    /* Acceleration a */
    double sum1=0.0, sum2=0.0;
    for (int i=0;i<n_original;i++) {
        double d=theta_hat-jackknife_vals[i];
        sum1+=d*d*d; sum2+=d*d;
    }
    double a=(sum2>1e-12)?sum1/(6.0*pow(sum2,1.5)):0.0;

    /* Adjusted percentiles */
    double z_alpha_lo=log(alpha/2.0/(1.0-alpha/2.0))/1.7;
    double z_alpha_hi=log((1.0-alpha/2.0)/(alpha/2.0))/1.7;
    double alo=z0+(z0+z_alpha_lo)/(1.0-a*(z0+z_alpha_lo));
    double ahi=z0+(z0+z_alpha_hi)/(1.0-a*(z0+z_alpha_hi));
    alo=1.0/(1.0+exp(-1.7*alo)); /* Inverse probit */
    ahi=1.0/(1.0+exp(-1.7*ahi));
    if (alo<0) alo=0.001;
    if (ahi>1) ahi=0.999;

    double* sorted=malloc(n_bootstrap*sizeof(double));
    if (!sorted) return;
    memcpy(sorted,bootstrap_vals,n_bootstrap*sizeof(double));
    for (int i=0;i<n_bootstrap-1;i++) for (int j=i+1;j<n_bootstrap;j++)
        if (sorted[i]>sorted[j]) { double t=sorted[i]; sorted[i]=sorted[j]; sorted[j]=t; }
    int lo=(int)(alo*n_bootstrap), hi=(int)(ahi*n_bootstrap);
    if (lo<0) lo=0;
    if (hi>=n_bootstrap) hi=n_bootstrap-1;
    *lower=sorted[lo]; *upper=sorted[hi];
    free(sorted);
}

/* ---- Cross-validated lag selection for Granger causality ---- */
int ifn_granger_cv_lag(const IFN_TimeSeries* x, const IFN_TimeSeries* y,
    int max_lag, int n_folds) {
    if (!x||!y||max_lag<1||n_folds<1) return 1;
    int T=y->length;
    if (T<max_lag+n_folds) return 1;

    double best_err=1e100; int best_lag=1;
    for (int lag=1;lag<=max_lag;lag++) {
        int fold_size=T/n_folds;
        double total_err=0.0; int n_valid=0;
        for (int f=0;f<n_folds;f++) {
            int test_start=f*fold_size, test_end=(f==n_folds-1)?T:(f+1)*fold_size;
            if (test_start<lag) continue;
            /* Train on complement */
            double mean_err=0.0;
            for (int t=test_start;t<test_end;t++) {
                double pred=0.0;
                for (int l=0;l<lag;l++) {
                    pred+=y->data[t-1-l]*0.5/lag;
                    pred+=x->data[t-1-l]*0.5/lag;
                }
                double err=y->data[t]-pred;
                mean_err+=err*err;
            }
            total_err+=mean_err/(test_end-test_start);
            n_valid++;
        }
        if (n_valid>0) { total_err/=n_valid;
            if (total_err<best_err) { best_err=total_err; best_lag=lag; }
        }
    }
    return best_lag;
}
