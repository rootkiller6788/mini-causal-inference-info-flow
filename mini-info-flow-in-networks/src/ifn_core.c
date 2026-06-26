#include "ifn_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

IFN_TimeSeries* ifn_ts_create(int length, int dim) {
    if (length<1||dim<1) return NULL;
    IFN_TimeSeries* ts=calloc(1,sizeof(IFN_TimeSeries));
    if (!ts) return NULL;
    ts->length=length; ts->dim=dim;
    ts->data=calloc(length*dim,sizeof(double));
    if (!ts->data) { free(ts); return NULL; }
    return ts;
}

void ifn_ts_free(IFN_TimeSeries* ts) {
    if (ts) { free(ts->data); free(ts->histogram); free(ts); }
}

void ifn_ts_normalize(IFN_TimeSeries* ts) {
    if (!ts) return;
    for (int d=0;d<ts->dim;d++) {
        double m=ifn_ts_mean(ts,d), s=ifn_ts_variance(ts,d);
        if (s<1e-12) s=1.0;
        for (int t=0;t<ts->length;t++)
            ts->data[t*ts->dim+d]=(ts->data[t*ts->dim+d]-m)/sqrt(s);
    }
}

double ifn_ts_mean(const IFN_TimeSeries* ts, int dim_idx) {
    if (!ts||dim_idx<0||dim_idx>=ts->dim) return 0.0;
    double s=0.0;
    for (int t=0;t<ts->length;t++) s+=ts->data[t*ts->dim+dim_idx];
    return s/ts->length;
}

double ifn_ts_variance(const IFN_TimeSeries* ts, int dim_idx) {
    if (!ts||dim_idx<0||dim_idx>=ts->dim) return 0.0;
    double m=ifn_ts_mean(ts,dim_idx), s=0.0;
    for (int t=0;t<ts->length;t++) {
        double d=ts->data[t*ts->dim+dim_idx]-m; s+=d*d;
    }
    return s/(ts->length-1);
}

void ifn_ts_histogram(IFN_TimeSeries* ts, int dim_idx, int n_bins) {
    if (!ts||dim_idx<0||dim_idx>=ts->dim||n_bins<1) return;
    free(ts->histogram);
    ts->histogram=calloc(n_bins+2,sizeof(double)); /* +2 for min/max */
    ts->n_bins=n_bins;
    if (!ts->histogram) return;
    double mn=ts->data[dim_idx], mx=mn;
    for (int t=1;t<ts->length;t++) {
        double v=ts->data[t*ts->dim+dim_idx];
        if (v<mn) mn=v;
        if (v>mx) mx=v;
    }
    double bw=(mx-mn+1e-10)/n_bins;
    ts->histogram[0]=mn; ts->histogram[1]=bw;
    for (int t=0;t<ts->length;t++) {
        int b=(int)((ts->data[t*ts->dim+dim_idx]-mn)/bw);
        if (b>=n_bins) b=n_bins-1;
        if (b<0) b=0;
        ts->histogram[2+b]++;
    }
}

double ifn_ts_entropy(IFN_TimeSeries* ts, int dim_idx) {
    if (!ts||dim_idx<0||dim_idx>=ts->dim) return 0.0;
    if (!ts->histogram) ifn_ts_histogram(ts,dim_idx,20);
    double H=0.0; int n=ts->length, nb=ts->n_bins;
    for (int b=0;b<nb;b++) {
        double p=ts->histogram[2+b]/n;
        if (p>1e-300) H-=p*log(p);
    }
    return H;
}

void ifn_ts_delay_embed(const IFN_TimeSeries* ts, int dim_idx,
    int tau, int m, double** embedded) {
    if (!ts||!embedded||dim_idx<0||dim_idx>=ts->dim||tau<1||m<1) return;
    int N=ts->length-(m-1)*tau;
    *embedded=malloc(N*m*sizeof(double));
    if (!*embedded) return;
    for (int i=0;i<N;i++)
        for (int j=0;j<m;j++)
            (*embedded)[i*m+j]=ts->data[(i+j*tau)*ts->dim+dim_idx];
}

double ifn_estimate_pdf_bins(const double* x, int n, int n_bins,
    double* bins, double* pdf) {
    if (!x||n<1||n_bins<1) return 0.0;
    double mn=x[0], mx=x[0];
    for (int i=1;i<n;i++) { if (x[i]<mn) mn=x[i]; if (x[i]>mx) mx=x[i]; }
    double bw=(mx-mn+1e-10)/n_bins;
    for (int b=0;b<n_bins;b++) { bins[b]=mn+(b+0.5)*bw; pdf[b]=0.0; }
    for (int i=0;i<n;i++) {
        int b=(int)((x[i]-mn)/bw);
        if (b>=n_bins) b=n_bins-1;
        if (b<0) b=0;
        pdf[b]++;
    }
    for (int b=0;b<n_bins;b++) pdf[b]/=n;
    return bw;
}

double ifn_estimate_joint_pdf(const double* x, const double* y, int n,
    int nx, int ny, double* joint) {
    if (!x||!y||!joint||n<1||nx<1||ny<1) return 0.0;
    memset(joint,0,nx*ny*sizeof(double));
    double mx=x[0],Mx=x[0]; for (int i=1;i<n;i++) { if (x[i]<mx) mx=x[i]; if (x[i]>Mx) Mx=x[i]; }
    double my=y[0],My=y[0]; for (int i=1;i<n;i++) { if (y[i]<my) my=y[i]; if (y[i]>My) My=y[i]; }
    double bwx=(Mx-mx+1e-10)/nx, bwy=(My-my+1e-10)/ny;
    for (int i=0;i<n;i++) {
        int bx=(int)((x[i]-mx)/bwx); if (bx>=nx) bx=nx-1; if (bx<0) bx=0;
        int by=(int)((y[i]-my)/bwy); if (by>=ny) by=ny-1; if (by<0) by=0;
        joint[bx*ny+by]++;
    }
    for (int i=0;i<nx*ny;i++) joint[i]/=n;
    return 0.0;
}

double ifn_estimate_conditional_pdf(const double* x, const double* y,
    const double* z, int n, int nx, int ny, int nz, double* cond) {
    if (!x||!y||!z||!cond||n<1) return 0.0;
    double* joint_xyz=calloc(nx*ny*nz,sizeof(double));
    double* joint_yz=calloc(ny*nz,sizeof(double));
    if (!joint_xyz||!joint_yz) { free(joint_xyz); free(joint_yz); return 0.0; }
    double mx=x[0],Mx=x[0]; for (int i=1;i<n;i++) { if (x[i]<mx) mx=x[i]; if (x[i]>Mx) Mx=x[i]; }
    double my=y[0],My=y[0]; for (int i=1;i<n;i++) { if (y[i]<my) my=y[i]; if (y[i]>My) My=y[i]; }
    double mz=z[0],Mz=z[0]; for (int i=1;i<n;i++) { if (z[i]<mz) mz=z[i]; if (z[i]>Mz) Mz=z[i]; }
    double bwx=(Mx-mx+1e-10)/nx, bwy=(My-my+1e-10)/ny, bwz=(Mz-mz+1e-10)/nz;
    for (int i=0;i<n;i++) {
        int bx=(int)((x[i]-mx)/bwx); if (bx>=nx) bx=nx-1; if (bx<0) bx=0;
        int by=(int)((y[i]-my)/bwy); if (by>=ny) by=ny-1; if (by<0) by=0;
        int bz=(int)((z[i]-mz)/bwz); if (bz>=nz) bz=nz-1; if (bz<0) bz=0;
        joint_xyz[(bx*ny+by)*nz+bz]++;
        joint_yz[by*nz+bz]++;
    }
    for (int i=0;i<nx*ny*nz;i++) joint_xyz[i]/=n;
    for (int i=0;i<ny*nz;i++) joint_yz[i]/=n;
    for (int bx=0;bx<nx;bx++)
        for (int by=0;by<ny;by++)
            for (int bz=0;bz<nz;bz++) {
                double p_yz=joint_yz[by*nz+bz];
                cond[(bx*ny+by)*nz+bz]=p_yz>1e-300?joint_xyz[(bx*ny+by)*nz+bz]/p_yz:0.0;
            }
    free(joint_xyz); free(joint_yz);
    return 0.0;
}

double ifn_mutual_information_binned(const double* x, const double* y, int n,
    int nx, int ny) {
    double* joint=calloc(nx*ny,sizeof(double));
    if (!joint) return 0.0;
    ifn_estimate_joint_pdf(x,y,n,nx,ny,joint);
    double* px=calloc(nx,sizeof(double)),*py=calloc(ny,sizeof(double));
    for (int i=0;i<nx;i++) for (int j=0;j<ny;j++) {
        px[i]+=joint[i*ny+j]; py[j]+=joint[i*ny+j];
    }
    double MI=0.0;
    for (int i=0;i<nx;i++) for (int j=0;j<ny;j++) {
        if (joint[i*ny+j]>1e-300 && px[i]>1e-300 && py[j]>1e-300)
            MI+=joint[i*ny+j]*log(joint[i*ny+j]/(px[i]*py[j]));
    }
    free(joint); free(px); free(py);
    return MI>0?MI:0.0;
}

double ifn_conditional_entropy_binned(const double* x, const double* y, int n,
    int nx, int ny) {
    double* joint=calloc(nx*ny,sizeof(double));
    ifn_estimate_joint_pdf(x,y,n,nx,ny,joint);
    double* py=calloc(ny,sizeof(double));
    for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) py[j]+=joint[i*ny+j];
    double H=0.0;
    for (int i=0;i<nx;i++) for (int j=0;j<ny;j++) {
        if (joint[i*ny+j]>1e-300 && py[j]>1e-300)
            H-=joint[i*ny+j]*log(joint[i*ny+j]/py[j]);
    }
    free(joint); free(py);
    return H>0?H:0.0;
}
