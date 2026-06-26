#ifndef TS_CORE_H
#define TS_CORE_H
#include <stdbool.h>
#include <stddef.h>

typedef struct { double* data; int n; } Vector;
typedef struct { double* data; int rows; int cols; } Matrix;

typedef struct { double* values; int length; int dim; } TimeSeries;

typedef struct {
    Matrix** A;
    double* c;
    int dim;
    int p;
    double aic;
    double bic;
    double sigma2;
} VARModel;

/* Vector ops */
Vector* vec_create(int n);
void vec_free(Vector* v);
double vec_norm(const Vector* v);
double vec_dot(const Vector* a, const Vector* b);
double vec_l1_norm(const Vector* v);
double vec_inf_norm(const Vector* v);
Vector* vec_add(const Vector* a, const Vector* b);
void vec_scale(Vector* v, double alpha);
double vec_mean(const Vector* v);
double vec_variance(const Vector* v);

/* Matrix ops */
Matrix* mat_create(int r, int c);
void mat_free(Matrix* m);
Matrix* mat_identity(int n);
Matrix* mat_copy(const Matrix* m);
Matrix* mat_mul(const Matrix* a, const Matrix* b);
Matrix* mat_transpose(const Matrix* m);
Matrix* mat_inverse(const Matrix* m);
double mat_det(const Matrix* m);
double mat_trace(const Matrix* m);
double mat_frobenius_norm(const Matrix* m);
double mat_condition_number(const Matrix* m);
Vector* mat_vec_mul(const Matrix* A, const Vector* x);
void mat_set(Matrix* m, int i, int j, double v);
double mat_get(const Matrix* m, int i, int j);

/* Time series */
TimeSeries* ts_create(int length, int dim);
void ts_free(TimeSeries* ts);
TimeSeries* ts_slice(const TimeSeries* ts, int start, int len);
double ols_estimate(const Matrix* X, const double* y, double* beta, double* resid, int n, int k);
void ts_print(const TimeSeries* ts);

/* VAR model */
VARModel* var_create(int dim, int p);
void var_free(VARModel* var);
int var_fit(VARModel* var, const TimeSeries* ts);
int var_forecast(const VARModel* var, const double* hist, double* fc, int steps);
int var_order_select(const TimeSeries* ts, int max_p, int dim, int* best_p, double* best_aic);

/* Test generators */
TimeSeries* ts_simulate_ar1(double phi, double sigma, int length);
TimeSeries* ts_simulate_var1_coupled(double a11,double a12,double a21,double a22,int len);
double* ts_extract_dim(const TimeSeries* ts, int dim_idx);

/* Autocorrelation */
double ts_autocorrelation(const double* x, int n, int k);
double* ts_acf(const double* x, int n, int max_lag);

#endif
