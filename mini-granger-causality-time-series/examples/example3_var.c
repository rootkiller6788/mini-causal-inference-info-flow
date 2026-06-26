#include "ts_core.h"
#include <stdio.h>
int main(void){
    printf("=== VAR Model Demo ===\n");
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.4, -0.2, 0.6, 200);
    int best_p; double best_aic;
    var_order_select(ts, 10, 2, &best_p, &best_aic);
    printf("Optimal lag order: %d (AIC=%.4f)\n", best_p, best_aic);
    VARModel* var = var_create(2, best_p);
    var_fit(var, ts);
    printf("Fitted VAR(%d): sigma2=%.6f\n", var->p, var->sigma2);
    printf("A[1] = [[%.3f, %.3f], [%.3f, %.3f]]\n",
           mat_get(var->A[0],0,0),mat_get(var->A[0],0,1),
           mat_get(var->A[0],1,0),mat_get(var->A[0],1,1));
    var_free(var); ts_free(ts);
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
