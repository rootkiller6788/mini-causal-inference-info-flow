/* Example 1: Sensitivity Analysis - Assessing Robustness to Unmeasured Confounding
 *
 * Demonstrates Rosenbaum bounds and E-value computation for a
 * matched-pair observational study. Shows how to determine if
 * a causal conclusion is robust to hidden bias.
 *
 * Scenario: Does a new teaching method improve test scores?
 * Matched pairs: 30 students, paired by pre-test scores.
 * Outcome: post-test score difference (treated - control).
 *
 * References:
 *   Rosenbaum, "Observational Studies", 2nd ed, 2002
 *   VanderWeele & Ding, Annals of Internal Medicine, 2017
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sensitivity_analysis.h"

int main(void) {
    srand((unsigned)time(NULL));
    printf("=== Example 1: Sensitivity Analysis ===\n\n");

    /* Generate matched-pair data: 30 pairs with true effect +2.0 */
    int n_pairs = 30;
    double *y_treated = calloc((size_t)n_pairs, sizeof(double));
    double *y_control = calloc((size_t)n_pairs, sizeof(double));
    if (!y_treated || !y_control) {
        fprintf(stderr, "Allocation failed\n");
        free(y_treated); free(y_control);
        return 1;
    }

    double true_effect = 2.0;
    for (int i = 0; i < n_pairs; i++) {
        double base = 50.0 + 15.0 * ((double)rand() / RAND_MAX);
        /* Unmeasured confounder effect: some pairs more affected */
        double conf = 1.5 * sin((double)i * 0.5);
        y_control[i] = base + conf;
        y_treated[i] = base + conf + true_effect +
                       2.0 * ((double)rand() / RAND_MAX - 0.5);
    }

    /* Compute Rosenbaum bounds for various Gamma values */
    double gammas[] = {1.0, 1.5, 2.0, 3.0, 5.0};
    int n_g = 5;
    printf("Rosenbaum Bounds Analysis:\n");
    printf("Gamma | p-upper | p-lower | ATE-range\n");
    printf("------+----------+---------+----------\n");
    for (int g = 0; g < n_g; g++) {
        RosenbaumBounds *rb = rosenbaum_bounds(y_treated, y_control,
                                                n_pairs, gammas[g]);
        if (rb) {
            printf(" %.1f  | %8.4f | %7.4f | [%6.2f, %6.2f]\n",
                   rb->gamma, rb->upper_bound, rb->lower_bound,
                   rb->lower_ate, rb->upper_ate);
            rosenbaum_free(rb);
        }
    }

    /* Critical Gamma */
    double crit_gamma = rosenbaum_critical_gamma(y_treated, y_control, n_pairs);
    printf("\nCritical Gamma: %.2f\n", crit_gamma);
    printf("Interpretation: Result remains significant for Gamma < %.2f\n",
           crit_gamma);

    /* E-value */
    double ate = 0.0;
    for (int i = 0; i < n_pairs; i++)
        ate += (y_treated[i] - y_control[i]);
    ate /= (double)n_pairs;

    double var_ate = 0.0;
    for (int i = 0; i < n_pairs; i++) {
        double d = (y_treated[i] - y_control[i]) - ate;
        var_ate += d * d;
    }
    var_ate /= (double)(n_pairs - 1);
    double se_ate = sqrt(var_ate / (double)n_pairs);

    EValue *ev = compute_evalue(ate, se_ate);
    if (ev) {
        printf("\nE-value Analysis:\n");
        printf("  Observed ATE: %.3f (SE=%.3f)\n", ate, se_ate);
        printf("  Risk ratio (approx): %.3f\n", ev->risk_ratio);
        printf("  E-value (point): %.2f\n", ev->e_value);
        printf("  E-value (CI): %.2f\n", ev->e_value_ci);
        printf("  Robust (E>1.5): %s\n", ev->is_robust ? "YES" : "NO");
        printf("\nInterpretation: An unmeasured confounder would need\n");
        printf("  RR >= %.2f with both treatment and outcome to explain away\n",
               ev->e_value);
        printf("  the observed effect.\n");
        evalue_free(ev);
    }

    free(y_treated); free(y_control);
    printf("\n=== Example 1 PASSED ===\n");
    return 0;
}
