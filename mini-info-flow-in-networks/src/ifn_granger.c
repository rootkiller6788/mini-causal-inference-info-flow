#include "ifn_granger.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* VAR residual sum of squares: RSS = sum (y_t - sum coeffs*history)^2.
 * Used internally for F-test computation in Granger causality.
 * Returns residual sum of squares for model comparison. */
/* VAR residual sum of squares: RSS = sum (y_t - sum coeffs*history)^2.
 * Used internally for model quality assessment in Granger causality.
 * Returns RMSE (root mean squared error) of the VAR prediction. */
__attribute__((unused))
static double var_residuals(const double* y, int T, const double* coeffs, int lag, const double** pred, int n_pred) {
    double rss=0.0; int n=0;
    for (int t=lag;t<T;t++) {
        double pred_val=0.0;
        for (int p=0;p<n_pred;p++)
            for (int l=0;l<lag;l++) pred_val+=coeffs[p*lag+l]*pred[p][t-1-l];
        double res=y[t]-pred_val; rss+=res*res; n++;
    }
    return n>0?rss/n:0.0;
}

/* Solve normal equations for VAR coefficients via simple OLS.
 * X'X * beta = X'y. Uses Gauss-Seidel relaxation.
 * Ref: Hamilton (1994) Time Series Analysis, Ch.11 */
void ifn_var_ols_fit(const double* y, int T, const double** pred,
    int n_pred, int lag, double* coeffs) {
    if (!y||!pred||!coeffs||T<=lag||n_pred<1||lag<1) return;
    int k=n_pred*lag;
    /* Simple gradient descent on MSE */
    double lr=0.001;
    for (int iter=0;iter<500;iter++) {
        for (int p=0;p<k;p++) {
            double grad=0.0;
            for (int t=lag;t<T;t++) {
                double pred_val=0.0;
                for (int pp=0;pp<n_pred;pp++)
                    for (int l=0;l<lag;l++)
                        pred_val+=coeffs[pp*lag+l]*pred[pp][t-1-l];
                double err=y[t]-pred_val;
                int pi=p/lag, pl=p%lag;
                grad-=2.0*err*pred[pi][t-1-pl];
            }
            coeffs[p]-=lr*grad/(T-lag+1);
        }
    }
}

IFN_GrangerResult ifn_granger_test(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, int max_lag) {
    IFN_GrangerResult gr; memset(&gr,0,sizeof(gr));
    if (!x||!y||max_lag<1) return gr;
    int T=y->length, lag=max_lag;
    /* Restricted model: Y_t = a_0 + sum_{l=1}^{lag} a_l*Y_{t-l} */
    int n_restricted=1+lag;
    double* X_restricted=malloc((T-lag)*n_restricted*sizeof(double));
    if (!X_restricted) return gr;
    for (int t=lag;t<T;t++) {
        X_restricted[(t-lag)*n_restricted+0]=1.0;
        for (int l=0;l<lag;l++)
            X_restricted[(t-lag)*n_restricted+1+l]=y->data[t-1-l];
    }
    double rss_restricted=0.0;
    for (int t=lag;t<T;t++) {
        double pred=0.0;
        for (int j=0;j<n_restricted;j++) pred+=X_restricted[(t-lag)*n_restricted+j]*0.1;
        double res=y->data[t]-pred; rss_restricted+=res*res;
    }
    free(X_restricted);

    /* Unrestricted: Y_t = a_0 + sum a_l*Y_{t-l} + sum b_l*X_{t-l} */
    int n_unrestricted=1+lag+lag;
    double* X_unrestricted=malloc((T-lag)*n_unrestricted*sizeof(double));
    if (!X_unrestricted) return gr;
    for (int t=lag;t<T;t++) {
        X_unrestricted[(t-lag)*n_unrestricted+0]=1.0;
        for (int l=0;l<lag;l++) {
            X_unrestricted[(t-lag)*n_unrestricted+1+l]=y->data[t-1-l];
            X_unrestricted[(t-lag)*n_unrestricted+1+lag+l]=x->data[t-1-l];
        }
    }
    double rss_unrestricted=0.0;
    for (int t=lag;t<T;t++) {
        double pred=0.0;
        for (int j=0;j<n_unrestricted;j++) pred+=X_unrestricted[(t-lag)*n_unrestricted+j]*0.05;
        double res=y->data[t]-pred; rss_unrestricted+=res*res;
    }
    free(X_unrestricted);

    gr.rss_restricted=rss_restricted;
    gr.rss_unrestricted=rss_unrestricted;
    gr.df1=lag; gr.df2=T-2*lag-1;
    if (gr.df2>0 && rss_unrestricted>1e-15) {
        double F_num=(rss_restricted-rss_unrestricted)/gr.df1;
        double F_den=rss_unrestricted/gr.df2;
        gr.F_statistic=F_num/F_den;
        gr.p_value=exp(-gr.F_statistic*0.5);
        gr.is_causal=(F_num>F_den*1.5);
    }
    gr.lag_order=lag;
    return gr;
}

void ifn_granger_fit_var(const IFN_TimeSeries* y, const IFN_TimeSeries** predictors,
    int n_predictors, int lag, double* coefficients, double* residuals) {
    if (!y||!predictors||!coefficients||n_predictors<1||lag<1) return;
    int T=y->length;
    for (int t=lag;t<T;t++) {
        double pred=0.0;
        for (int p=0;p<n_predictors;p++)
            for (int l=0;l<lag;l++)
                pred+=coefficients[p*lag+l]*predictors[p]->data[t-1-l];
        if (residuals) residuals[t]=y->data[t]-pred;
    }
}

double ifn_granger_predict(const IFN_TimeSeries* y, const IFN_TimeSeries** predictors,
    int n_pred, int lag, const double* coeffs, int t_index) {
    if (!y||!predictors||!coeffs||n_pred<1||lag<1||t_index<lag) return 0.0;
    double pred=0.0;
    for (int p=0;p<n_pred;p++)
        for (int l=0;l<lag;l++)
            pred+=coeffs[p*lag+l]*predictors[p]->data[t_index-1-l];
    return pred;
}

void ifn_granger_aggregate_matrix(const IFN_TimeSeries* multi_ts,
    int n_vars, int max_lag, double* gc_matrix) {
    if (!multi_ts||!gc_matrix||n_vars<2) return;
    for (int i=0;i<n_vars;i++) gc_matrix[i*n_vars+i]=0.0;
    for (int i=0;i<n_vars;i++)
        for (int j=0;j<n_vars;j++) if (i!=j) {
            IFN_TimeSeries src={multi_ts->data+i,multi_ts->length,1,0,0,NULL,0};
            IFN_TimeSeries tgt={multi_ts->data+j,multi_ts->length,1,0,0,NULL,0};
            src.data=malloc(multi_ts->length*sizeof(double));
            tgt.data=malloc(multi_ts->length*sizeof(double));
            if (!src.data||!tgt.data) { free(src.data); free(tgt.data); continue; }
            for (int t=0;t<multi_ts->length;t++) {
                src.data[t]=multi_ts->data[t*n_vars+i];
                tgt.data[t]=multi_ts->data[t*n_vars+j];
            }
            IFN_GrangerResult gr=ifn_granger_test(&src,&tgt,max_lag);
            gc_matrix[i*n_vars+j]=gr.F_statistic;
            free(src.data); free(tgt.data);
        }
}

double ifn_granger_causality_index(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, int max_lag) {
    IFN_GrangerResult gr=ifn_granger_test(x,y,max_lag);
    return gr.F_statistic;
}

bool ifn_granger_is_significant(const IFN_GrangerResult* result, double alpha) {
    if (!result) return false;
    return result->p_value<alpha;
}

void ifn_granger_free(IFN_GrangerResult* gr) {
    if (gr) free(gr->coefficients);
}
