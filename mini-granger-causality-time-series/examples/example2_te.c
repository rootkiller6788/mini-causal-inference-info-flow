#include "ts_core.h"
#include "transfer_entropy.h"
#include <stdio.h>
#include <stdlib.h>
int main(void){
    printf("=== Transfer Entropy Demo ===\n");
    TimeSeries* ts = ts_simulate_var1_coupled(0.5, 0.4, 0.0, 0.5, 1000);
    double* x = ts_extract_dim(ts, 0);
    double* y = ts_extract_dim(ts, 1);
    TransferEntropy* te = te_create();
    te_significance(te, x, y, 1000, 10, 3, 100, 0.05);
    printf("X -> Y TE: %.6f (sig=%s, p=%.4f)\n",
           te->te_xy, te->is_significant?"YES":"NO", te->p_value);
    te_significance(te, y, x, 1000, 10, 3, 100, 0.05);
    printf("Y -> X TE: %.6f (sig=%s, p=%.4f)\n",
           te->te_xy, te->is_significant?"YES":"NO", te->p_value);
    te_free(te); free(x); free(y); ts_free(ts);
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
