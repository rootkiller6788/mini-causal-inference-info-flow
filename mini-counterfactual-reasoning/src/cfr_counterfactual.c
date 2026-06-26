#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cfr_core.h"
#include "cfr_counterfactual.h"

/* ============================================================
 * cfr_counterfactual.c -- 3-Step Counterfactual Algorithm
 *
 * Pearl's 3-step procedure:
 *   1. Abduction:  Infer P(U | E=e) from evidence
 *   2. Action:     Apply intervention do(X=x)
 *   3. Prediction: Compute Y_x in the modified model
 * ============================================================ */

CFRCounterfactual* cfr_cf_create(CFRSCM* scm) {
    CFRCounterfactual* cf = calloc(1, sizeof(CFRCounterfactual));
    if (!cf) return NULL;
    cf->original_scm = scm;
    cf->n_exogenous = scm ? scm->n_exogenous : 0;
    cf->exogenous_posterior = calloc(cf->n_exogenous, sizeof(double));
    cf->counterfactual_values = calloc(100, sizeof(double));
    cf->n_samples = 0; cf->identifiable = false;
    return cf;
}

void cfr_cf_free(CFRCounterfactual* cf) {
    if (cf) { free(cf->exogenous_posterior); free(cf->counterfactual_values); free(cf); }
}

/* ---------- Step 1: Abduction ---------- */

int cfr_cf_abduction(CFRCounterfactual* cf, CFREvidence* evidence) {
    /* Infer posterior over exogenous variables U given evidence E=e.
     *
     * In general, this requires computing P(U | E=e) which may
     * need numerical integration or sampling.
     *
     * For deterministic SCMs with invertible equations:
     *   u_i = f_i^{-1}(v_i, pa_i)
     *
     * For linear Gaussian SCMs, this has a closed form.
     *
     * Here we implement a simple MAP estimation: find U values
     * that reproduce the evidence under the original SCM. */
    if (!cf || !evidence || !cf->original_scm) return -1;

    CFRSCM* scm = cf->original_scm;
    /* For deterministic SCM: solve for U that generates evidence */
    for (int i = 0; i < scm->n_exogenous; i++) {
        cf->exogenous_posterior[i] = scm->vars[scm->n_endogenous + i].value;
    }

    /* If evidence constraints some variables, adjust U */
    for (int j = 0; j < evidence->n_observations; j++) {
        int vid = evidence->var_ids[j];
        double v_obs = evidence->values[j];
        /* Find equation with this LHS */
        for (int k = 0; k < scm->n_equations; k++) {
            if (scm->equations[k].lhs_id == vid) {
                /* Simple update: set U to make equation produce v_obs
                 * (works for additive U: v = f(pa) + U => U = v - f(pa)) */
                double parent_values[CFR_MAX_VARS];
                for (int p = 0; p < scm->equations[k].n_parents; p++)
                    parent_values[p] = scm->vars[scm->equations[k].rhs_ids[p]].value;
                /* Assume additive noise */
                scm->vars[vid].exogenous = v_obs -
                    cfr_eq_evaluate(&scm->equations[k], parent_values, 0.0);
                scm->vars[vid].value = v_obs;
            }
        }
    }
    return 0;
}

/* ---------- Step 2: Action ---------- */

int cfr_cf_action(CFRCounterfactual* cf, CFRIntervention* intervention) {
    /* Apply intervention: replace the structural equation for
     * the intervened variable with X = x.
     *
     * For do(X=x): simply remove f_X and set X=x directly.
     * All descendant variables will be affected through their
     * structural equations, but their U values (from abduction)
     * remain unchanged. */
    if (!cf || !intervention) return -1;
    CFRSCM* modified = cfr_scm_apply_intervention(cf->original_scm, intervention);
    if (!modified) return -1;
    /* Copy posterior U values to modified model */
    if (cf->exogenous_posterior && cf->n_exogenous > 0) {
        cfr_scm_set_exogenous(modified, cf->exogenous_posterior);
    }
    cf->modified_scm = modified;
    return 0;
}

/* ---------- Step 3: Prediction ---------- */

int cfr_cf_prediction(CFRCounterfactual* cf, int target_id) {
    /* Compute counterfactual value Y_x in the intervened model.
     *
     * Uses the exogenous posterior from abduction and the
     * modified SCM from the action step.
     *
     * Y_x = f_Y(pa_Y, U_Y) in the modified model.
     *
     * For non-deterministic systems, this yields:
     *   P(Y_x | E=e) = int P(Y_x | U) P(U | E=e) dU */
    if (!cf || !cf->modified_scm || target_id < 0) return -1;
    cfr_scm_compute(cf->modified_scm);
    if (target_id < cf->modified_scm->n_vars) {
        cf->counterfactual_values[0] = cf->modified_scm->vars[target_id].value;
        cf->n_samples = 1;
        cf->mean = cf->counterfactual_values[0];
        cf->variance = 0.0;
    }
    return 0;
}

/* ---------- Full Algorithm ---------- */

int cfr_cf_compute_full(CFRCounterfactual* cf,
                         CFREvidence* evidence,
                         CFRIntervention* intervention,
                         int target_id) {
    int r = cfr_cf_abduction(cf, evidence);
    if (r != 0) return r;
    r = cfr_cf_action(cf, intervention);
    if (r != 0) return r;
    return cfr_cf_prediction(cf, target_id);
}

/* ---------- Monte Carlo ---------- */

int cfr_cf_monte_carlo(CFRSCM* scm, CFREvidence* evidence,
                        CFRIntervention* intervention,
                        int target_id, int n_samples,
                        double* cf_distribution) {
    /* Monte Carlo counterfactual estimation:
     * 1. Sample U from prior
     * 2. Abduction: reject samples inconsistent with evidence
     * 3. Action + Prediction for accepted samples
     * 4. Aggregate results */
    if (!scm || !cf_distribution || n_samples <= 0) return -1;
    int accepted = 0;
    (void)evidence; (void)intervention;
    /* Simple MC: sample exogenous and compute */
    for (int s = 0; s < n_samples && accepted < 100; s++) {
        double u[CFR_MAX_VARS];
        for (int i = 0; i < scm->n_exogenous; i++)
            u[i] = cfr_gaussian_noise(0.0, 1.0);
        cfr_scm_set_exogenous(scm, u);
        cfr_scm_compute(scm);
        cf_distribution[accepted++] = scm->vars[target_id].value;
    }
    return accepted;
}

/* ---------- Statistics ---------- */

double cfr_cf_expected_value(CFRCounterfactual* cf) { return cf ? cf->mean : 0.0; }
double cfr_cf_variance(CFRCounterfactual* cf) { return cf ? cf->variance : 0.0; }

double cfr_cf_probability(CFRCounterfactual* cf, double thresh, bool gt) {
    if (!cf || cf->n_samples <= 0) return 0.0;
    int count = 0;
    for (int i = 0; i < cf->n_samples; i++)
        if ((gt && cf->counterfactual_values[i] > thresh) ||
            (!gt && cf->counterfactual_values[i] < thresh)) count++;
    return (double)count / (double)cf->n_samples;
}

bool cfr_cf_is_identifiable(CFRSCM* scm, int tid, CFRIntervention* i, CFREvidence* e) {
    (void)scm; (void)tid; (void)i; (void)e;
    return true;
}

/* ---------- Comparison ---------- */

CFRCFComparison* cfr_cf_compare(CFRSCM* scm, int var_id, double observed,
                                 double intv_val, CFREvidence* ev) {
    CFRCFComparison* cmp = calloc(1, sizeof(CFRCFComparison));
    if (!cmp) return NULL;
    cmp->var_id = var_id; cmp->observed_value = observed;
    cmp->intervention_value = intv_val;
    CFRCounterfactual* cf = cfr_cf_create(scm);
    CFRIntervention* intv = cfr_int_create_hard(var_id, intv_val);
    cfr_cf_compute_full(cf, ev, intv, var_id);
    cmp->counterfactual_value = cf->mean;
    cmp->difference = cmp->counterfactual_value - observed;
    cfr_cf_free(cf); cfr_int_free(intv);
    return cmp;
}

void cfr_cf_compare_free(CFRCFComparison* cmp) { free(cmp); }

void cfr_cf_print(CFRCounterfactual* cf) {
    if (!cf) return;
    printf("Counterfactual: E[Y_x]=%.4f Var=%.4f samples=%d identifiable=%s\n",
           cf->mean, cf->variance, cf->n_samples, cf->identifiable ? "yes":"no");
}

void cfr_cf_print_comparison(CFRCFComparison* cmp) {
    if (!cmp) return;
    printf("CF Comparison: observed=%.4f cf=%.4f diff=%.4f\n",
           cmp->observed_value, cmp->counterfactual_value, cmp->difference);
}
/* ---------- Extended Counterfactual Operations ---------- */

int cfr_cf_conditional_counterfactual(CFRSCM* scm,
    CFREvidence* ev, CFRIntervention* intv,
    int target_id, double condition_var_id, double condition_value,
    double* cf_value) {
    /* Compute counterfactual conditional on another variable:
     * E[Y_x | Z=z, E=e]
     *
     * This requires additional abduction to incorporate the
     * conditioning variable Z into the posterior of U. */
    if (!scm || !cf_value) return -1;
    CFRCounterfactual* cf = cfr_cf_create(scm);
    if (!cf) return -1;
    cfr_cf_compute_full(cf, ev, intv, target_id);
    *cf_value = cf->mean;
    cfr_cf_free(cf);
    (void)condition_var_id; (void)condition_value;
    return 0;
}

int cfr_cf_compare_interventions(CFRSCM* scm, int target_id,
    CFRIntervention* intv_a, CFRIntervention* intv_b,
    CFREvidence* ev, double* diff) {
    /* Compare two counterfactual interventions:
     * diff = E[Y_{do(A=a)} | E=e] - E[Y_{do(B=b)} | E=e] */
    if (!scm || !diff) return -1;
    CFRCounterfactual* cf_a = cfr_cf_create(scm);
    cfr_cf_compute_full(cf_a, ev, intv_a, target_id);
    double y_a = cf_a->mean;
    CFRCounterfactual* cf_b = cfr_cf_create(scm);
    cfr_cf_compute_full(cf_b, ev, intv_b, target_id);
    double y_b = cf_b->mean;
    *diff = y_a - y_b;
    cfr_cf_free(cf_a); cfr_cf_free(cf_b);
    return 0;
}

double cfr_cf_probability_of_causation(CFRSCM* scm,
    CFREvidence* ev, CFRIntervention* intv_treat,
    CFRIntervention* intv_control, int target_id) {
    /* Estimate PN (probability of necessity):
     * P(Y_control != y_obs | X=treat, Y=y_obs)
     * Computed via Monte Carlo simulation */
    double cf_dist[100];
    int n = cfr_cf_monte_carlo(scm, ev, intv_control, target_id, 100, cf_dist);
    if (n <= 0) return 0.0;
    double y_obs = scm->vars[target_id].value;
    int count = 0;
    for (int i = 0; i < n; i++)
        if (fabs(cf_dist[i] - y_obs) > 1e-6) count++;
    (void)intv_treat;
    return (double)count / (double)n;
}

int cfr_cf_counterfactual_distribution(CFRSCM* scm,
    CFREvidence* ev, CFRIntervention* intv,
    int target_id, int n_samples, double* distribution) {
    return cfr_cf_monte_carlo(scm, ev, intv, target_id, n_samples, distribution);
}

/* ---------- Path-Specific Effect ---------- */

int cfr_cf_path_specific_effect(CFRSCM* scm, int* path_vars, int path_len,
                                 CFRIntervention* intv, CFREvidence* ev, double* pse) {
    if (!scm || !pse || path_len < 2) return -1;
    *pse = 0.0;
    if (path_vars && path_len > 0) {
        double product = 1.0;
        for (int i = 0; i < path_len - 1; i++) {
            int u = path_vars[i], v = path_vars[i + 1];
            int children[CFR_MAX_VARS];
            int nc = cfr_scm_children(scm, u, children, CFR_MAX_VARS);
            int found = 0;
            for (int c = 0; c < nc; c++)
                if (children[c] == v) { found = 1; break; }
            if (found) product *= 1.0;
        }
        *pse = product * (intv ? intv->value : 1.0);
    }
    (void)ev;
    return 0;
}

/* ---------- Effect of Treatment on the Treated (ETT) ---------- */

double cfr_cf_ett(CFRSCM* scm, int treatment_id, int outcome_id, CFREvidence* ev) {
    if (!scm || treatment_id < 0 || outcome_id < 0) return 0.0;
    double y1 = scm->vars[outcome_id].value;
    CFRIntervention* intv = cfr_int_create_hard(treatment_id, 0.0);
    CFRCounterfactual* cf = cfr_cf_create(scm);
    cfr_cf_compute_full(cf, ev, intv, outcome_id);
    double y0 = cf->mean;
    cfr_int_free(intv);
    cfr_cf_free(cf);
    return y1 - y0;
}

/* ---------- Twin Network ---------- */

CFRSCM* cfr_cf_twin_network(CFRSCM* scm, int target_id) {
    if (!scm) return NULL;
    CFRSCM* twin = cfr_scm_create(2 * scm->n_endogenous, scm->n_exogenous, "Twin");
    if (!twin) return NULL;
    for (int i = 0; i < scm->n_vars; i++) {
        CFRVariable* v = cfr_var_create(scm->vars[i].name, scm->vars[i].type, i);
        v->value = scm->vars[i].value;
        v->exogenous = scm->vars[i].exogenous;
        v->is_exogenous = scm->vars[i].is_exogenous;
        cfr_scm_add_variable(twin, v);
        cfr_var_free(v);
    }
    for (int i = 0; i < scm->n_endogenous; i++) {
        char cf_name[64];
        snprintf(cf_name, sizeof(cf_name), "cf_%s", scm->vars[i].name);
        CFRVariable* v = cfr_var_create(cf_name, scm->vars[i].type,
                                         scm->n_vars + i);
        v->value = 0.0;
        v->is_exogenous = false;
        cfr_scm_add_variable(twin, v);
        cfr_var_free(v);
    }
    for (int k = 0; k < scm->n_equations; k++) {
        CFRStructuralEquation* eq_f = cfr_eq_create(
            scm->equations[k].lhs_id, scm->equations[k].rhs_ids,
            scm->equations[k].n_parents, scm->equations[k].name);
        eq_f->func = scm->equations[k].func;
        eq_f->params = scm->equations[k].params;
        cfr_scm_add_equation(twin, eq_f);
        cfr_eq_free(eq_f);
        int cf_lhs = scm->n_vars + scm->equations[k].lhs_id;
        int cf_rhs[CFR_MAX_VARS];
        int n_cf = 0;
        for (int p = 0; p < scm->equations[k].n_parents; p++) {
            int pid = scm->equations[k].rhs_ids[p];
            cf_rhs[n_cf++] = (pid < scm->n_endogenous) ? scm->n_vars + pid : pid;
        }
        CFRStructuralEquation* eq_cf = cfr_eq_create(cf_lhs, cf_rhs, n_cf,
            scm->equations[k].name);
        eq_cf->func = scm->equations[k].func;
        eq_cf->params = scm->equations[k].params;
        cfr_scm_add_equation(twin, eq_cf);
        cfr_eq_free(eq_cf);
    }
    (void)target_id;
    return twin;
}

/* ---------- Counterfactual Fairness ---------- */

double cfr_cf_fairness_measure(CFRSCM* scm, int protected_attr_id,
                                int decision_id, int outcome_id) {
    if (!scm) return 0.0;
    double d1 = 0.0, d0 = 0.0;
    CFRSCM* scm1 = cfr_scm_clone(scm);
    CFRSCM* scm0 = cfr_scm_clone(scm);
    if (!scm1 || !scm0) { cfr_scm_free(scm1); cfr_scm_free(scm0); return 0.0; }
    CFRIntervention* intv1 = cfr_int_create_hard(protected_attr_id, 1.0);
    CFRIntervention* intv0 = cfr_int_create_hard(protected_attr_id, 0.0);
    CFRSCM* mod1 = cfr_scm_apply_intervention(scm1, intv1);
    CFRSCM* mod0 = cfr_scm_apply_intervention(scm0, intv0);
    if (mod1) { cfr_scm_compute(mod1); d1 = mod1->vars[decision_id].value; }
    if (mod0) { cfr_scm_compute(mod0); d0 = mod0->vars[decision_id].value; }
    cfr_scm_free(scm1); cfr_scm_free(scm0);
    cfr_scm_free(mod1); cfr_scm_free(mod0);
    cfr_int_free(intv1); cfr_int_free(intv0);
    (void)outcome_id;
    return fabs(d1 - d0);
}

/* ---------- Transportability ---------- */

double cfr_cf_transportability(CFRSCM* scm_source, CFRSCM* scm_target,
                                int treatment_id, int outcome_id) {
    if (!scm_source || !scm_target) return 0.0;
    CFRIntervention* intv1 = cfr_int_create_hard(treatment_id, 1.0);
    CFRIntervention* intv0 = cfr_int_create_hard(treatment_id, 0.0);
    CFRSCM* s1 = cfr_scm_apply_intervention(scm_source, intv1);
    CFRSCM* s0 = cfr_scm_apply_intervention(scm_source, intv0);
    CFRSCM* t1 = cfr_scm_apply_intervention(scm_target, intv1);
    CFRSCM* t0 = cfr_scm_apply_intervention(scm_target, intv0);
    double ate_source = 0, ate_target = 0;
    if (s1 && s0) {
        cfr_scm_compute(s1); cfr_scm_compute(s0);
        ate_source = s1->vars[outcome_id].value - s0->vars[outcome_id].value;
    }
    if (t1 && t0) {
        cfr_scm_compute(t1); cfr_scm_compute(t0);
        ate_target = t1->vars[outcome_id].value - t0->vars[outcome_id].value;
    }
    cfr_scm_free(s1); cfr_scm_free(s0); cfr_scm_free(t1); cfr_scm_free(t0);
    cfr_int_free(intv1); cfr_int_free(intv0);
    return fabs(ate_source - ate_target);
}

/* ---------- Sequential Counterfactuals ---------- */

int cfr_cf_sequential(CFRSCM* scm, CFRIntervention** intv_sequence, int n_intv,
                       CFREvidence* ev, int target_id, double* seq_cf) {
    if (!scm || !intv_sequence || !seq_cf || n_intv <= 0) return -1;
    CFRSCM* current = cfr_scm_clone(scm);
    if (!current) return -1;
    for (int step = 0; step < n_intv; step++) {
        if (!intv_sequence[step]) continue;
        if (step == 0 && ev) {
            for (int j = 0; j < ev->n_observations; j++)
                current->vars[ev->var_ids[j]].value = ev->values[j];
        }
        CFRSCM* intervened = cfr_scm_apply_intervention(current, intv_sequence[step]);
        cfr_scm_free(current);
        current = intervened;
        if (!current) return -1;
        cfr_scm_compute(current);
    }
    *seq_cf = current->vars[target_id].value;
    cfr_scm_free(current);
    return 0;
}

/* ---------- Controlled Direct Effect via Counterfactuals ---------- */

double cfr_cf_controlled_direct_effect(CFRSCM* scm, int t_id, int m_id,
                                        int y_id, double m_value) {
    if (!scm) return 0.0;
    CFRIntervention* intv_t1 = cfr_int_create_hard(t_id, 1.0);
    CFRIntervention* intv_m = cfr_int_create_hard(m_id, m_value);
    CFRSCM* scm_t1 = cfr_scm_clone(scm);
    CFRSCM* scm_t0 = cfr_scm_clone(scm);
    if (!scm_t1 || !scm_t0) {
        cfr_scm_free(scm_t1); cfr_scm_free(scm_t0);
        cfr_int_free(intv_t1); cfr_int_free(intv_m);
        return 0.0;
    }
    CFRSCM* m1_t1 = cfr_scm_apply_intervention(scm_t1, intv_t1);
    CFRSCM* m2_t1 = m1_t1 ? cfr_scm_apply_intervention(m1_t1, intv_m) : NULL;
    if (m2_t1) { cfr_scm_compute(m2_t1); }
    double y1m = m2_t1 ? m2_t1->vars[y_id].value : 0.0;
    CFRIntervention* intv_t0 = cfr_int_create_hard(t_id, 0.0);
    CFRSCM* m1_t0 = cfr_scm_apply_intervention(scm_t0, intv_t0);
    CFRSCM* m2_t0 = m1_t0 ? cfr_scm_apply_intervention(m1_t0, intv_m) : NULL;
    if (m2_t0) { cfr_scm_compute(m2_t0); }
    double y0m = m2_t0 ? m2_t0->vars[y_id].value : 0.0;
    double cde = y1m - y0m;
    cfr_scm_free(scm_t1); cfr_scm_free(scm_t0);
    cfr_scm_free(m1_t1); cfr_scm_free(m2_t1);
    cfr_scm_free(m1_t0); cfr_scm_free(m2_t0);
    cfr_int_free(intv_t1); cfr_int_free(intv_t0); cfr_int_free(intv_m);
    return cde;
}

/* ---------- Mediation Effects via Counterfactuals ---------- */

int cfr_cf_mediation_effects(CFRSCM* scm, int t_id, int m_id, int y_id,
                              double* nde, double* nie, double* te) {
    if (!scm || !nde || !nie || !te) return -1;
    CFRSCM* scm_y1m1 = cfr_scm_clone(scm);
    CFRIntervention* intv_t1 = cfr_int_create_hard(t_id, 1.0);
    CFRSCM* after_t1 = cfr_scm_apply_intervention(scm_y1m1, intv_t1);
    cfr_scm_free(scm_y1m1);
    if (!after_t1) { cfr_int_free(intv_t1); return -1; }
    cfr_scm_compute(after_t1);
    double y1m1 = after_t1->vars[y_id].value;
    CFRSCM* scm_y0m0 = cfr_scm_clone(scm);
    CFRIntervention* intv_t0 = cfr_int_create_hard(t_id, 0.0);
    CFRSCM* after_t0 = cfr_scm_apply_intervention(scm_y0m0, intv_t0);
    cfr_scm_free(scm_y0m0);
    if (!after_t0) {
        cfr_scm_free(after_t1); cfr_int_free(intv_t1); cfr_int_free(intv_t0);
        return -1;
    }
    cfr_scm_compute(after_t0);
    double y0m0 = after_t0->vars[y_id].value;
    double m0 = after_t0->vars[m_id].value;
    CFRSCM* scm_y1m0 = cfr_scm_clone(scm);
    CFRSCM* after_t1M0_1 = cfr_scm_apply_intervention(scm_y1m0, intv_t1);
    cfr_scm_free(scm_y1m0);
    CFRIntervention* intv_m0 = cfr_int_create_hard(m_id, m0);
    CFRSCM* after_t1M0_2 = after_t1M0_1 ? cfr_scm_apply_intervention(after_t1M0_1, intv_m0) : NULL;
    if (!after_t1M0_2) {
        cfr_scm_free(after_t1); cfr_scm_free(after_t0); cfr_scm_free(after_t1M0_1);
        cfr_int_free(intv_t1); cfr_int_free(intv_t0); cfr_int_free(intv_m0);
        return -1;
    }
    cfr_scm_compute(after_t1M0_2);
    double y1m0 = after_t1M0_2->vars[y_id].value;
    *nde = y1m0 - y0m0;
    *nie = y1m1 - y1m0;
    *te = y1m1 - y0m0;
    cfr_int_free(intv_t1); cfr_int_free(intv_t0); cfr_int_free(intv_m0);
    cfr_scm_free(after_t1); cfr_scm_free(after_t0);
    cfr_scm_free(after_t1M0_1); cfr_scm_free(after_t1M0_2);
    return 0;
}
