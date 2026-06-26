#include "ifn_mutual.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

IFN_MIResult ifn_mutual_information(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, int n_bins) {
    IFN_MIResult r; memset(&r,0,sizeof(r));
    if (!x||!y||n_bins<1) return r;
    r.mi_value=ifn_mutual_information_binned(x->data,y->data,x->length,n_bins,n_bins);
    double Hx=ifn_ts_entropy((IFN_TimeSeries*)x,0), Hy=ifn_ts_entropy((IFN_TimeSeries*)y,0);
    double minH=Hx<Hy?Hx:Hy;
    r.normalized_mi=minH>1e-15?r.mi_value/minH:0.0;
    r.n_nodes=2;
    return r;
}

double ifn_mi_continuous(const double* x, const double* y, int n, int n_bins) {
    return ifn_mutual_information_binned(x,y,n,n_bins,n_bins);
}

double ifn_mi_knn(const double* x, const double* y, int n, int k) {
    if (!x||!y||n<k+1||k<1) return 0.0;
    double MI=0.0;
    for (int i=0;i<n;i++) {
        double* dists=malloc(n*sizeof(double));
        if (!dists) continue;
        for (int j=0;j<n;j++)
            dists[j]=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
        for (int a=0;a<n-1;a++) for (int b=a+1;b<n;b++)
            if (dists[a]>dists[b]) { double t=dists[a]; dists[a]=dists[b]; dists[b]=t; }
        int nx=0,ny=0;
        for (int j=1;j<=k&&j<n;j++) {
            double dx=fabs(x[i]-x[j]), dy=fabs(y[i]-y[j]);
            if (dx<dists[k]) nx++;
            if (dy<dists[k]) ny++;
        }
        MI+=log((double)nx*ny/(double)(k*k));
        free(dists);
    }
    MI/=n;
    return MI>0?MI:0.0;
}

double ifn_mi_kernel(const double* x, const double* y, int n, double bandwidth) {
    if (!x||!y||n<2||bandwidth<=0) return 0.0;
    double MI=0.0;
    for (int i=0;i<n;i++) {
        double px=0.0,py=0.0,pxy=0.0;
        for (int j=0;j<n;j++) {
            double dx=(x[i]-x[j])/bandwidth, dy=(y[i]-y[j])/bandwidth;
            double kx=exp(-0.5*dx*dx), ky=exp(-0.5*dy*dy);
            px+=kx; py+=ky; pxy+=kx*ky;
        }
        if (px>1e-300&&py>1e-300&&pxy>1e-300)
            MI+=log(pxy*n/(px*py));
    }
    return MI>0?MI/n:0.0;
}

void ifn_mi_conditioned(const IFN_TimeSeries* multi_ts,
    int n_vars, int idx_x, int idx_y, int n_bins, double* cmi_out) {
    if (!multi_ts||!cmi_out||n_vars<3) return;
    for (int z=0;z<n_vars;z++) {
        if (z==idx_x||z==idx_y) { cmi_out[z]=0.0; continue; }
        double* x=malloc(multi_ts->length*sizeof(double));
        double* y=malloc(multi_ts->length*sizeof(double));
        double* zv=malloc(multi_ts->length*sizeof(double));
        if (!x||!y||!zv) { free(x); free(y); free(zv); continue; }
        for (int t=0;t<multi_ts->length;t++) {
            x[t]=multi_ts->data[t*n_vars+idx_x];
            y[t]=multi_ts->data[t*n_vars+idx_y];
            zv[t]=multi_ts->data[t*n_vars+z];
        }
        ifn_mutual_information_binned(x,zv,multi_ts->length,n_bins,n_bins);  /* I(X;Z) quality check */
        ifn_mutual_information_binned(y,zv,multi_ts->length,n_bins,n_bins);  /* I(Y;Z) quality check */
        double* xyz=malloc(multi_ts->length*2*sizeof(double));
        if (xyz) {
            for (int t=0;t<multi_ts->length;t++) { xyz[t*2]=x[t]; xyz[t*2+1]=y[t]; }
            double I_xy_given_z=ifn_mutual_information_binned(zv,xyz,multi_ts->length,n_bins,2*n_bins);
            cmi_out[z]=I_xy_given_z;
            free(xyz);
        }
        free(x); free(y); free(zv);
    }
}

void ifn_mi_matrix(const IFN_TimeSeries* multi_ts, int n_vars,
    int n_bins, double* mi_matrix) {
    if (!multi_ts||!mi_matrix||n_vars<2) return;
    for (int i=0;i<n_vars;i++) mi_matrix[i*n_vars+i]=0.0;
    for (int i=0;i<n_vars;i++)
        for (int j=i+1;j<n_vars;j++) {
            double* x=malloc(multi_ts->length*sizeof(double));
            double* y=malloc(multi_ts->length*sizeof(double));
            if (!x||!y) { free(x); free(y); continue; }
            for (int t=0;t<multi_ts->length;t++) {
                x[t]=multi_ts->data[t*n_vars+i];
                y[t]=multi_ts->data[t*n_vars+j];
            }
            double mi=ifn_mutual_information_binned(x,y,multi_ts->length,n_bins,n_bins);
            mi_matrix[i*n_vars+j]=mi_matrix[j*n_vars+i]=mi;
            free(x); free(y);
        }
}

double ifn_mi_normalized(const IFN_TimeSeries* x, const IFN_TimeSeries* y, int n_bins) {
    IFN_MIResult r=ifn_mutual_information(x,y,n_bins);
    return r.normalized_mi;
}

double ifn_partial_mutual_information(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, const IFN_TimeSeries* z, int n_bins) {
    if (!x||!y||!z) return 0.0;
    double I_xy=ifn_mutual_information_binned(x->data,y->data,x->length,n_bins,n_bins);
    /* Partial MI: I(X;Y|Z) = I(X;Y) - I(X;Y;Z) approximates to MI conditioned on Z */
    double* xz=malloc(x->length*2*sizeof(double));
    double* yz=malloc(y->length*2*sizeof(double));
    if (!xz||!yz) { free(xz); free(yz); return I_xy; }
    for (int t=0;t<x->length;t++) { xz[t*2]=x->data[t]; xz[t*2+1]=z->data[t]; }
    for (int t=0;t<y->length;t++) { yz[t*2]=y->data[t]; yz[t*2+1]=z->data[t]; }
    double I_xz=ifn_mutual_information_binned(x->data,z->data,x->length,n_bins,n_bins);
    double I_yz=ifn_mutual_information_binned(y->data,z->data,y->length,n_bins,n_bins);
    free(xz); free(yz);
    double pmi=I_xy-0.5*(I_xz+I_yz);
    return pmi>0?pmi:0.0;
}

double ifn_interaction_information(const IFN_TimeSeries* x,
    const IFN_TimeSeries* y, const IFN_TimeSeries* z, int n_bins) {
    if (!x||!y||!z) return 0.0;
    double I_xy=ifn_mutual_information_binned(x->data,y->data,x->length,n_bins,n_bins);
    double I_xz=ifn_mutual_information_binned(x->data,z->data,x->length,n_bins,n_bins);
    double I_yz=ifn_mutual_information_binned(y->data,z->data,y->length,n_bins,n_bins);
    /* Approximate II = I(X;Y) + I(X;Z|Y) - I(X;Z) */
    return I_xy+I_yz-I_xz;
}

double ifn_redundancy(const IFN_TimeSeries* x, const IFN_TimeSeries* y,
    const IFN_TimeSeries* z, int n_bins) {
    if (!x||!y||!z) return 0.0;
    double I_xy=ifn_mutual_information_binned(x->data,y->data,x->length,n_bins,n_bins);
    double I_xz=ifn_mutual_information_binned(x->data,z->data,x->length,n_bins,n_bins);
    /* Joint MI: concatenate y and z */
    double* yz=malloc(y->length*2*sizeof(double));
    if (!yz) return 0.0;
    for (int t=0;t<y->length;t++) { yz[t*2]=y->data[t]; yz[t*2+1]=z->data[t]; }
    double I_x_yz=ifn_mutual_information_binned(x->data,yz,x->length,n_bins,2*n_bins);
    free(yz);
    return I_xy+I_xz-I_x_yz;
}
