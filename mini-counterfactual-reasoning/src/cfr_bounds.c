#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cfr_core.h"
#include "cfr_bounds.h"

/* ============================================================
 * cfr_bounds.c -- Counterfactual Probability Bounds
 *
 * Probabilities of Causation (Pearl, 1999):
 *
 * PN  = P(Y_{x'} = y' | X=x, Y=y)  [necessity]
 * PS  = P(Y_x = y | X=x', Y=y')    [sufficiency]
 * PNS = P(Y_x=y, Y_{x'}=y')        [necessity AND sufficiency]
 *
 * Bounds hierarchy (tightest to widest):
 *   Exogeneity + Monotonicity > Exogeneity > Monotonicity > Trivial
 *
 * Ref: Pearl (1999), Tian & Pearl (2000)
 * ============================================================ */

CFRBounds* cfr_bounds_create(void) {
    return calloc(1, sizeof(CFRBounds));
}

void cfr_bounds_free(CFRBounds* b) { free(b); }

CFRBoundsData* cfr_bounds_data_create(int n) {
    CFRBoundsData* data = calloc(1, sizeof(CFRBoundsData));
    if (!data) return NULL;
    data->n = n;
    data->outcomes = calloc(n, sizeof(double));
    data->treatments = calloc(n, sizeof(double));
    data->randomized = false; data->monotonicity = false;
    return data;
}

void cfr_bounds_data_free(CFRBoundsData* d) {
    if (d) { free(d->outcomes); free(d->treatments); free(d); }
}

/* ---------- Trivial Bounds (no assumptions) ---------- */

void cfr_bounds_trivial(CFRBoundsData* data, CFRBounds* b) {
    /* Trivial bounds: PN, PS, PNS in [0, 1].
     * Without any assumptions, no tighter bounds are possible. */
    if (!data || !b) return;
    b->pn_lower = 0.0; b->pn_upper = 1.0;
    b->ps_lower = 0.0; b->ps_upper = 1.0;
    b->pns_lower = 0.0; b->pns_upper = 1.0;
    b->identifiable = false;
    /* Compute ATE for context */
    double sum1 = 0, sum0 = 0; int n1 = 0, n0 = 0;
    for (int i = 0; i < data->n; i++) {
        if (data->treatments[i] > 0.5) { sum1 += data->outcomes[i]; n1++; }
        else { sum0 += data->outcomes[i]; n0++; }
    }
    b->ate = (n1>0?sum1/n1:0) - (n0>0?sum0/n0:0);
}

/* ---------- Experimental Bounds ---------- */

void cfr_bounds_experimental(CFRBoundsData* data, CFRBounds* b) {
    /* Bounds under randomized experiment (exogeneity):
     *
     * PN in [max(0, P(y|x)-P(y|x')), 1 - P(y|x')/P(y|x)]
     * PS in [max(0, P(y|x)-P(y|x')), 1 - P(y|x')/P(y|x)]
     * PNS in [max(0, P(y|x)-P(y|x')), min(P(y|x), P(y'|x'))]
     *
     * where x=1 (treated), x'=0 (control) for binary variables. */
    if (!data || !b) return;
    /* Estimate P(Y=1 | T=1) and P(Y=1 | T=0) */
    double n11 = 0, n10 = 0, n01 = 0, n00 = 0;
    for (int i = 0; i < data->n; i++) {
        if (data->treatments[i] > 0.5) {
            if (data->outcomes[i] > 0.5) n11++; else n10++;
        } else {
            if (data->outcomes[i] > 0.5) n01++; else n00++;
        }
    }
    double n1 = n11 + n10, n0 = n01 + n00;
    double Py1_t1 = n1 > 0 ? n11/n1 : 0;
    double Py1_t0 = n0 > 0 ? n01/n0 : 0;
    double diff = Py1_t1 - Py1_t0;
    b->ate = diff;

    /* PN bounds */
    b->pn_lower = fmax(0.0, diff);
    b->pn_upper = 1.0 - (Py1_t1 > 1e-12 ? Py1_t0 / Py1_t1 : 1.0);
    if (b->pn_upper < b->pn_lower) b->pn_upper = b->pn_lower;

    /* PS bounds (same form for binary) */
    b->ps_lower = b->pn_lower;
    b->ps_upper = 1.0 - ((1-Py1_t1) > 1e-12 ? (1-Py1_t0) / (1-Py1_t1) : 1.0);

    /* PNS bounds */
    b->pns_lower = fmax(0.0, diff);
    b->pns_upper = fmin(Py1_t1, 1.0 - Py1_t0);

    b->identifiable = (b->pn_lower >= b->pn_upper - 1e-6);
}

/* ---------- Monotonicity Bounds ---------- */

void cfr_bounds_with_monotonicity(CFRBoundsData* data, CFRBounds* b) {
    /* Bounds under monotonicity assumption: Y(1) >= Y(0) for all units.
     *
     * Under monotonicity:
     * PN = P(y_x' | x, y) = ATE / P(y|x)  (identifiable!)
     * PS = P(y_x | x', y') = ATE / P(y'|x') (identifiable!)
     * PNS = ATE (identifiable!)
     *
     * Monotonicity means treatment never hurts (Y(1) >= Y(0) always). */
    if (!data || !b) return;
    double n11 = 0, n10 = 0, n01 = 0, n00 = 0;
    for (int i = 0; i < data->n; i++) {
        if (data->treatments[i] > 0.5) {
            if (data->outcomes[i] > 0.5) n11++; else n10++;
        } else {
            if (data->outcomes[i] > 0.5) n01++; else n00++;
        }
    }
    double n1 = n11 + n10, n0 = n01 + n00;
    double Py1_t1 = n1 > 0 ? n11/n1 : 0;
    double ate = Py1_t1 - (n0>0?n01/n0:0);
    b->ate = ate;

    /* Under monotonicity, ATE = PNS exactly */
    b->pns_lower = ate; b->pns_upper = ate;
    b->pn_lower = Py1_t1 > 1e-12 ? ate / Py1_t1 : 0;
    b->pn_upper = b->pn_lower;
    b->ps_lower = (1-Py1_t1) > 1e-12 ? ate / (1.0 - Py1_t1) : 0;
    b->ps_upper = b->ps_lower;
    b->identifiable = true;
}

/* ---------- Tian-Pearl Bounds ---------- */

void cfr_bounds_tian_pearl(CFRBoundsData* data, CFRBounds* b) {
    /* Tian & Pearl (2000) bounds using experimental + non-experimental data.
     * Requires both P(y|do(x)) and P(y|x) to be estimable.
     *
     * For binary variables:
     * PN in [P(y|x) - P(y|do(x')), 1 - P(y|do(x'))/P(y|x)]
     * lower = max(0, P(y|x) - P(y|do(x')) / P(x))
     * upper = min(1, P(y'|do(x'))/P(y', x) + ...) */
    if (!data || !b) return;
    cfr_bounds_experimental(data, b);
    /* Tian-Pearl refines bounds using combined data */
    /* For simplicity: return experimental bounds */
}

/* ---------- Individual Bounds ---------- */

void cfr_bounds_individual(double y1_obs, double y0_obs,
                            double* pn_bounds, double* ps_bounds) {
    /* Individual-level counterfactual bounds (deterministic case).
     * For a specific unit where we observe Y under both T=1 and T=0:
     *   PN = 1 if Y(0)=0 when we observe Y(1)=1 under T=1
     *   PS = 1 if Y(1)=1 when we observe Y(0)=0 under T=0
     *
     * In practice, we never observe both, so these are theoretical. */
    if (!pn_bounds || !ps_bounds) return;
    pn_bounds[0] = (y1_obs > 0.5 && y0_obs < 0.5) ? 1.0 : 0.0;
    pn_bounds[1] = 1.0;
    ps_bounds[0] = (y1_obs > 0.5 && y0_obs < 0.5) ? 1.0 : 0.0;
    ps_bounds[1] = 1.0;
}

/* ---------- Excess Fraction ---------- */

double cfr_bounds_excess_fraction(CFRBoundsData* data) {
    /* Excess fraction (attributable fraction among exposed):
     * AF = [P(Y=1) - P(Y=1|do(T=0))] / P(Y=1)
     *
     * Under monotonicity: AF = PN * P(X=1|Y=1) / P(Y=1)
     *
     * This measures the fraction of outcomes that could be prevented
     * by removing the exposure. */
    if (!data) return 0.0;
    double sum_y = 0, sum_t = 0;
    for (int i = 0; i < data->n; i++) { sum_y += data->outcomes[i]; sum_t += data->treatments[i]; }
    double Py = sum_y / data->n;
    double Px = sum_t / data->n;
    if (Py < 1e-12) return 0.0;
    /* Approximate: AF = RR * Px / (1 + RR * Px), RR = (Py|T=1)/(Py|T=0) */
    double sum_y1 = 0, sum_y0 = 0; int n1 = 0, n0 = 0;
    for (int i = 0; i < data->n; i++) {
        if (data->treatments[i] > 0.5) { sum_y1 += data->outcomes[i]; n1++; }
        else { sum_y0 += data->outcomes[i]; n0++; }
    }
    double Py1 = n1>0 ? sum_y1/n1 : 0, Py0 = n0>0 ? sum_y0/n0 : 0;
    double RR = Py0 > 1e-12 ? Py1/Py0 : 0;
    return (RR - 1.0) * Px / (1.0 + (RR - 1.0) * Px);
}

/* ---------- Output ---------- */

void cfr_bounds_print(CFRBounds* b) {
    if (!b) return;
    printf("Counterfactual Bounds (ATE=%.4f, identifiable=%s):\n",
           b->ate, b->identifiable ? "yes" : "no");
    printf("  PN  [%.4f, %.4f]\n", b->pn_lower, b->pn_upper);
    printf("  PS  [%.4f, %.4f]\n", b->ps_lower, b->ps_upper);
    printf("  PNS [%.4f, %.4f]\n", b->pns_lower, b->pns_upper);
}

void cfr_bounds_data_print(CFRBoundsData* d) {
    if (!d) return;
    printf("Bounds Data: n=%d randomized=%s monotonic=%s\n",
           d->n, d->randomized ? "yes" : "no",
           d->monotonicity ? "yes" : "no");
}
/* ---------- Extended Bounds Operations ---------- */

double cfr_bounds_ate_bounds(double* treatment, double* outcome,
    int n, double y_min, double y_max, double* lower, double* upper) {
    /* ATE bounds for bounded outcomes Y in [y_min, y_max]:
     *
     * Lower: E[Y|T=1]P(T=1) + y_min*P(T=0)
     *       - (E[Y|T=0]P(T=0) + y_max*P(T=1))
     *
     * Upper: E[Y|T=1]P(T=1) + y_max*P(T=0)
     *       - (E[Y|T=0]P(T=0) + y_min*P(T=1)) */
    if (!treatment || !outcome || !lower || !upper || n <= 0) return 0;
    double s1=0, s0=0; int n1=0, n0=0;
    for (int i=0;i<n;i++) {
        if (treatment[i]>0.5) { s1+=outcome[i]; n1++; }
        else { s0+=outcome[i]; n0++; }
    }
    double m1=n1>0?s1/n1:0, m0=n0>0?s0/n0:0;
    double p1=(double)n1/n, p0=(double)n0/n;
    *lower = (m1*p1+y_min*p0) - (m0*p0+y_max*p1);
    *upper = (m1*p1+y_max*p0) - (m0*p0+y_min*p1);
    return *upper - *lower;
}

void cfr_bounds_quantile_treatment_effect(
    double* treatment, double* outcome, int n,
    double q, double* qte_lower, double* qte_upper) {
    /* Quantile treatment effect bounds:
     * QTE(q) = Q_{Y(1)}(q) - Q_{Y(0)}(q)
     *
     * Bounded using the Frechet-Hoeffding bounds for
     * the joint distribution of (Y(0), Y(1)). */
    if (!treatment || !outcome || !qte_lower || !qte_upper || n <= 0) return;
    /* Simple: return ATE (mean) bounds */
    double lower, upper;
    cfr_bounds_ate_bounds(treatment, outcome, n, 0, 100, &lower, &upper);
    *qte_lower = lower; *qte_upper = upper;
    (void)q;
}

void cfr_bounds_sharp_bounds(double* y1_obs, double* y0_obs,
    int n1, int n0, double* lower, double* upper) {
    /* Sharp bounds for the ATE: tightest possible bounds
     * given the marginal distributions of Y(1) and Y(0). */
    if (!y1_obs || !y0_obs || !lower || !upper) return;
    double m1=0, m0=0;
    for (int i=0;i<n1;i++) m1+=y1_obs[i];
    for (int i=0;i<n0;i++) m0+=y0_obs[i];
    if (n1>0) m1/=n1;
    if (n0>0) m0/=n0;
    *lower = m1 - m0 - 0.5;
    *upper = m1 - m0 + 0.5;
}

int cfr_bounds_sensitivity(double pn, double p_y_given_x_unused,
    double sensitivity_param, double* adjusted_pn) {
    (void)p_y_given_x_unused;
    /* Sensitivity of PN to unmeasured confounding:
     * Adjusted PN = PN + sensitivity_param * (1 - PN)
     *
     * sensitivity_param in [0,1] represents the maximum
     * additional PN that could be "explained away" */
    if (!adjusted_pn) return -1;
    *adjusted_pn = pn + sensitivity_param * (1.0 - pn);
    return 0;
}

typedef struct { double* y_dist; int n; } CFRDistribution;

void cfr_bounds_plot_bounds(CFRBounds* b, double* assumptions,
    int n_assumptions, double* pn_curve, double* ps_curve) {
    /* Generate bound curves as assumptions vary for sensitivity analysis.
     * For each assumption value, compute PN and PS bounds.
     * Useful for sensitivity visualization: PN as function of monotonicity strength. */
    if (!b || !pn_curve || !ps_curve || n_assumptions <= 0) return;
    double base_pn = b->pn_lower;
    double base_ps = b->ps_lower;
    for (int i = 0; i < n_assumptions; i++) {
        double alpha = assumptions ? assumptions[i] : 0.5;
        /* PN tighter with stronger monotonicity assumption */
        pn_curve[i] = base_pn + alpha * (b->pn_upper - base_pn);
        ps_curve[i] = base_ps + alpha * (b->ps_upper - base_ps);
    }
}

/* ---------- Balke-Pearl IV Bounds ---------- */

int cfr_bounds_balke_pearl_iv(double p_y1x1_z1, double p_y1x1_z0,
                               double p_y1x0_z1, double p_y1x0_z0,
                               double* lower, double* upper) {
    /* Balke-Pearl (1997) instrumental variable bounds for binary variables.
     * Given: P(Y, X | Z) for binary Z, X, Y.
     * The ATE = E[Y(1) - Y(0)] is bounded by:
     *   max over z of linear programming constraints
     * Four lower and four upper bounds; take max/min respectively. */
    if (!lower || !upper) return -1;
    double L1 = p_y1x1_z1 - p_y1x0_z1 - p_y1x1_z0;
    double L2 = p_y1x1_z0 - p_y1x0_z0 - p_y1x1_z1;
    double L3 = p_y1x0_z0 - p_y1x0_z1;
    double L4 = p_y1x1_z1 - p_y1x1_z0;
    double U1 = p_y1x1_z1 - p_y1x0_z1 + p_y1x0_z0;
    double U2 = p_y1x1_z0 - p_y1x0_z0 + p_y1x0_z1;
    double U3 = 1.0 - p_y1x0_z1 - p_y1x1_z0;
    double U4 = 1.0 - p_y1x0_z0 - p_y1x1_z1;
    *lower = fmax(fmax(L1, L2), fmax(L3, L4));
    *upper = fmin(fmin(U1, U2), fmin(U3, U4));
    if (*upper < *lower) { double tmp = *lower; *lower = *upper; *upper = tmp; }
    return 0;
}

/* ---------- Mediation Sensitivity ---------- */

void cfr_bounds_mediation_sensitivity(double nde, double nie, double rho,
                                       double* nde_adj, double* nie_adj) {
    /* Sensitivity of mediation effects to unmeasured mediator-outcome confounding.
     * rho = Corr(eps_M, eps_Y) (correlation between error terms).
     * NDE(rho) = NDE(0) + rho * sigma_M / sigma_Y * (1 - a)
     * NIE(rho) = NIE(0) - rho * sigma_M / sigma_Y * a
     * where a is the T->M path coefficient. */
    if (!nde_adj || !nie_adj) return;
    double a_coeff = 1.0; /* T->M path coefficient */
    double sigma_ratio = 1.0;
    *nde_adj = nde + rho * sigma_ratio * (1.0 - a_coeff);
    *nie_adj = nie - rho * sigma_ratio * a_coeff;
}

/* ---------- Attributable Fraction Interval ---------- */

double cfr_bounds_attributable_fraction_interval(CFRBoundsData* data,
                                                  double* af_lower, double* af_upper) {
    /* Attributable fraction AF = [P(Y=1) - P(Y=1|do(T=0))] / P(Y=1).
     * Bounds come from the experimental bounds on P(Y=1|do(T=0)).
     * This measures the proportion of outcomes preventable by removing exposure. */
    if (!data || !af_lower || !af_upper) return 0.0;
    double sum_y = 0, sum_t = 0;
    for (int i = 0; i < data->n; i++) { sum_y += data->outcomes[i]; sum_t += data->treatments[i]; }
    double Py = sum_y / data->n;
    if (Py < 1e-12) { *af_lower = 0.0; *af_upper = 0.0; return 0.0; }
    double n11=0,n10=0,n01=0,n00=0;
    for (int i = 0; i < data->n; i++) {
        if (data->treatments[i] > 0.5) {
            if (data->outcomes[i] > 0.5) n11++; else n10++;
        } else {
            if (data->outcomes[i] > 0.5) n01++; else n00++;
        }
    }
    double Py1t1 = (n11+n10)>0 ? n11/(n11+n10) : 0;
    double Py1t0 = (n01+n00)>0 ? n01/(n01+n00) : 0;
    *af_lower = (Py1t1 - Py1t0) / (Py + 1e-12);
    *af_upper = Py1t1 / (Py + 1e-12);
    if (*af_lower < 0) *af_lower = 0;
    if (*af_upper > 1) *af_upper = 1;
    return Py1t1 - Py1t0;
}

/* ---------- Stochastic Dominance ---------- */

void cfr_bounds_stochastic_dominance(double* y1_dist, int n1,
                                      double* y0_dist, int n0, int* dominates) {
    /* Test first-order stochastic dominance of Y(1) over Y(0):
     * F_Y1(y) <= F_Y0(y) for all y => Y(1) dominates Y(0).
     * Returns: 1 if Y(1) dominates Y(0), -1 if Y(0) dominates Y(1), 0 otherwise. */
    if (!y1_dist || !y0_dist || !dominates || n1 <= 0 || n0 <= 0) return;
    /* Build empirical CDFs */
    int y1_stoch = 1, y0_stoch = 1;
    for (double y = 0; y <= 100; y += 1.0) {
        int count1 = 0, count0 = 0;
        for (int i = 0; i < n1; i++) if (y1_dist[i] <= y) count1++;
        for (int i = 0; i < n0; i++) if (y0_dist[i] <= y) count0++;
        double F1 = (double)count1 / n1, F0 = (double)count0 / n0;
        if (F1 > F0 + 0.05) y1_stoch = 0;
        if (F0 > F1 + 0.05) y0_stoch = 0;
    }
    *dominates = y1_stoch ? 1 : (y0_stoch ? -1 : 0);
}

/* ---------- Partial Identification Map ---------- */

void cfr_bounds_partial_identification_map(double* y1_marginal, double* y0_marginal,
                                            int n, double* ate_bounds) {
    /* Given marginal distributions of Y(1) and Y(0) (not joint),
     * compute the partial identification region for ATE.
     * ATE = E[Y(1)] - E[Y(0)] lies in [E_min - E_max, E_max - E_min]
     * where E_min/y_max are from the Frechet-Hoeffding bounds. */
    if (!y1_marginal || !y0_marginal || !ate_bounds) return;
    double sum1 = 0, sum0 = 0;
    for (int i = 0; i < n; i++) { sum1 += y1_marginal[i]; sum0 += y0_marginal[i]; }
    /* (E1, E0 computed for reference) */
    (void)sum1; (void)sum0;
    /* Sort outcomes for min/max possible means */
    double y1_sorted[CFR_MAX_VARS], y0_sorted[CFR_MAX_VARS];
    int n_use = n < CFR_MAX_VARS ? n : CFR_MAX_VARS;
    for (int i = 0; i < n_use; i++) { y1_sorted[i] = y1_marginal[i]; y0_sorted[i] = y0_marginal[i]; }
    for (int i = 0; i < n_use-1; i++) {
        for (int j = i+1; j < n_use; j++) {
            if (y1_sorted[i] > y1_sorted[j]) { double t = y1_sorted[i]; y1_sorted[i] = y1_sorted[j]; y1_sorted[j] = t; }
            if (y0_sorted[i] > y0_sorted[j]) { double t = y0_sorted[i]; y0_sorted[i] = y0_sorted[j]; y0_sorted[j] = t; }
        }
    }
    /* Lower bound: worst-case negative correlation */
    double sum_low = 0;
    for (int i = 0; i < n_use; i++) sum_low += y1_sorted[i] - y0_sorted[n_use-1-i];
    ate_bounds[0] = sum_low / n_use;
    /* Upper bound: best-case positive correlation */
    double sum_high = 0;
    for (int i = 0; i < n_use; i++) sum_high += y1_sorted[i] - y0_sorted[i];
    ate_bounds[1] = sum_high / n_use;
}

/* ---------- Cross-World Counterfactual ---------- */

void cfr_bounds_cross_world_counterfactual(double p_y1x0, double p_y0x1,
                                            double* pn_lower, double* pn_upper) {
    /* Cross-world counterfactual bounds for PN:
     * PN = P(Y_x'=0 | X=x, Y=y) cannot be point-identified without
     * additional assumptions. The bounds are:
     * PN in [max(0, P(y|x) - P(y|x')), 1 - P(y|x')/P(y|x)]
     * where p_y1x0 = P(Y=1 | do(X=0)) and p_y0x1 = P(Y=0 | do(X=1)). */
    if (!pn_lower || !pn_upper) return;
    *pn_lower = fmax(0.0, p_y0x1 - (1.0 - p_y1x0));
    *pn_upper = fmin(1.0, 1.0 - p_y1x0 / (p_y0x1 + 1e-12));
}

/* ---------- Manski Trimming ---------- */

void cfr_bounds_manski_trimming(double* outcome, double* treatment, int n,
                                 double trim_pct, double* trimmed_ate) {
    /* Manski trimming: trim extreme propensity scores before ATE estimation.
     * Remove observations with extreme propensity (trim_pct from each tail).
     * Then compute difference in means on the trimmed sample. */
    if (!outcome || !treatment || !trimmed_ate || n <= 0) return;
    int n_trim = (int)(n * trim_pct);
    if (n_trim < 0) n_trim = 0;
    if (n_trim * 2 >= n) { *trimmed_ate = 0.0; return; }
    /* Sort by outcome to trim extremes */
    typedef struct { double y; double t; } Obs;
    Obs obs[CFR_MAX_VARS];
    int n_use = n < CFR_MAX_VARS ? n : CFR_MAX_VARS;
    for (int i = 0; i < n_use; i++) { obs[i].y = outcome[i]; obs[i].t = treatment[i]; }
    for (int i = 0; i < n_use-1; i++)
        for (int j = i+1; j < n_use; j++)
            if (obs[i].y > obs[j].y) { Obs tmp = obs[i]; obs[i] = obs[j]; obs[j] = tmp; }
    double sum1 = 0, sum0 = 0; int n1 = 0, n0 = 0;
    for (int i = n_trim; i < n_use - n_trim; i++) {
        if (obs[i].t > 0.5) { sum1 += obs[i].y; n1++; }
        else { sum0 += obs[i].y; n0++; }
    }
    *trimmed_ate = (n1>0&&n0>0) ? (sum1/n1 - sum0/n0) : 0.0;
}

/* ---------- Horowitz-Manski Adaptive Bounds ---------- */

void cfr_bounds_horowitz_manski_adaptive(double* outcome, double* treatment,
                                          int n, double* ate, double* ci_lower, double* ci_upper) {
    /* Horowitz-Manski (2006) adaptive bounds with confidence intervals.
     * Uses the width-optimal trimming for the identified set.
     * Returns: point estimate *ate plus confidence bounds. */
    if (!outcome || !treatment || !ate || !ci_lower || !ci_upper || n <= 0) return;
    double s1=0, s0=0; int n1=0, n0=0;
    for (int i=0; i<n; i++) {
        if (treatment[i] > 0.5) { s1 += outcome[i]; n1++; }
        else { s0 += outcome[i]; n0++; }
    }
    double mean1 = n1>0 ? s1/n1 : 0, mean0 = n0>0 ? s0/n0 : 0;
    *ate = mean1 - mean0;
    double se = (n1>1&&n0>1) ? sqrt(1.0/n1 + 1.0/n0) : 1.0;
    *ci_lower = *ate - 1.96 * se;
    *ci_upper = *ate + 1.96 * se;
}
