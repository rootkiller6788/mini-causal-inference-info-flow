#include "ts_core.h"
#include "granger_test.h"
#include <stdio.h>
#include <stdlib.h>
int main(void){
    printf("=== Granger Causality Test Demo ===\n");
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.4, 0.0, 0.5, 500);
    double* x = ts_extract_dim(ts, 0);
    double* y = ts_extract_dim(ts, 1);
    GrangerTestResult* g_xy = gt_create();
    GrangerTestResult* g_yx = gt_create();
    gt_test_bidirectional(g_xy, g_yx, x, y, 500, 10, 0.05);
    printf("X -> Y: "); gt_print(g_xy);
    printf("Y -> X: "); gt_print(g_yx);
    gt_free(g_xy); gt_free(g_yx); free(x); free(y); ts_free(ts);
    return 0;
}
/* Reference note 1: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 2: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 3: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 4: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 5: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 6: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 7: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 8: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 9: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 10: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 11: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 12: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 13: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 14: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
/* Reference note 15: Pearl (2009) Causality; Morgan & Winship (2015) Counterfactuals and Causal Inference. */
