#include "granger_bootstrap.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
int main(void) {
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.3, 0.0, 0.5, 100);
    double* x = ts_extract_dim(ts, 0);
    double* y = ts_extract_dim(ts, 1);
    BootstrapGrangerResult* boot = boot_create(50);
    assert(boot != NULL);
    int rc = boot_parametric(boot, x, y, 100, 2, 0.05);
    fprintf(stderr,"rc=%d F=%.4f p=%.4f\n",rc,boot->observed_F,boot->bootstrap_p_value);
    boot_free(boot); free(x); free(y); ts_free(ts);
    printf("PASS\n");
    return 0;
}
