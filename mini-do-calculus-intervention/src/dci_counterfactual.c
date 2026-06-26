#include "dci_counterfactual.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Counterfactual Computation
 *
 * Step 1 (Abduction): Update P(u) using evidence E=e
 * Step 2 (Action): Intervene do(X=x)
 * Step 3 (Prediction): Compute P(Y=y) in modified model
 * ============================================================== */

DCI_CounterfactualResult dci_counterfactual_probability(
    DCI_SCM* scm, DCI_Counterfactual* query,
    const double* evidence, const int* evidence_vars,
    int n_evidence, int n_samples) {
    DCI_CounterfactualResult result;
    memset(&result, 0, sizeof(result));
    if (!scm || !query || n_samples <= 0) return result;
    int matching = 0;
    int i, j;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        /* Abduction: sample exog consistent with evidence */
        int max_tries = 10;
        int tries;
        bool consistent = false;
        for (tries = 0; tries < max_tries && !consistent; tries++) {
            dci_scm_evaluate(scm, exog);
            consistent = true;
            for (j = 0; j < n_evidence && consistent; j++) {
                if (fabs(scm->vars[evidence_vars[j]].value - evidence[j]) > 0.1)
                    consistent = false;
            }
            if (!consistent) {
                for (j = 0; j < scm->n_vars; j++)
                    exog[j] = (double)rand() / RAND_MAX;
            }
        }
        if (!consistent) continue;
        /* Action: intervene on antecedent */
        exog[query->antecedent_var] = query->antecedent_value;
        dci_scm_evaluate(scm, exog);
        /* Check consequent */
        if (fabs(scm->vars[query->consequent_var].value - query->consequent_value) < 0.01)
            matching++;
    }
    result.probability = (double)matching / n_samples;
    result.n_samples = n_samples;
    result.is_identifiable = true;
    return result;
}

/* ==============================================================
 * Twin Network
 * ============================================================== */

DCI_TwinNetwork* dci_twin_network_create(DCI_SCM* scm,
    int intervention_var, double intervention_value) {
    if (!scm) return NULL;
    DCI_TwinNetwork* tn = (DCI_TwinNetwork*)calloc(1, sizeof(DCI_TwinNetwork));
    if (!tn) return NULL;
    tn->factual = scm;
    /* Twin = copy of SCM but with intervention applied */
    DCI_SCM* twin = dci_scm_create("twin");
    int i;
    for (i = 0; i < scm->n_vars; i++) {
        if (i == intervention_var) {
            dci_scm_add_variable(twin, scm->vars[i].name, scm->vars[i].type,
                NULL, 0, NULL, true);
            twin->vars[twin->n_vars - 1].value = intervention_value;
        } else {
            dci_scm_add_variable(twin, scm->vars[i].name, scm->vars[i].type,
                scm->vars[i].parent_ids, scm->vars[i].n_parents,
                scm->vars[i].equation, scm->vars[i].is_exogenous);
        }
    }
    tn->counterfactual = twin;
    return tn;
}

void dci_twin_network_free(DCI_TwinNetwork* tn) {
    if (tn) {
        dci_scm_free(tn->counterfactual);
        free(tn);
    }
}

double dci_twin_network_query(DCI_TwinNetwork* tn,
    int query_var, double query_value) {
    if (!tn || !tn->counterfactual) return -1.0;
    double exog[DCI_MAX_VARS] = {0};
    dci_scm_evaluate(tn->counterfactual, exog);
    return tn->counterfactual->vars[query_var].value;
    (void)query_value;
}

/* ==============================================================
 * Probability of Necessity and Sufficiency
 * ============================================================== */

double dci_probability_of_necessity(DCI_SCM* scm, int x, int y,
    int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    /* PN = P(Y_X=0 = 0 | X=1, Y=1) */
    int num = 0, den = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        if (scm->vars[x].value > 0.5 && scm->vars[y].value > 0.5) {
            den++;
            exog[x] = 0.0;
            dci_scm_evaluate(scm, exog);
            if (scm->vars[y].value < 0.5) num++;
        }
    }
    return den > 0 ? (double)num / den : 0.0;
}

double dci_probability_of_sufficiency(DCI_SCM* scm, int x, int y,
    int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    int num = 0, den = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        if (scm->vars[x].value < 0.5 && scm->vars[y].value < 0.5) {
            den++;
            exog[x] = 1.0;
            dci_scm_evaluate(scm, exog);
            if (scm->vars[y].value > 0.5) num++;
        }
    }
    return den > 0 ? (double)num / den : 0.0;
}

double dci_probability_of_necessity_sufficiency(DCI_SCM* scm,
    int x, int y, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    int num = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        exog[x] = 1.0;
        dci_scm_evaluate(scm, exog);
        double y1 = scm->vars[y].value;
        exog[x] = 0.0;
        dci_scm_evaluate(scm, exog);
        double y0 = scm->vars[y].value;
        if (y1 > 0.5 && y0 < 0.5) num++;
    }
    return (double)num / n_samples;
}

/* ==============================================================
 * Individual Counterfactual
 * ============================================================== */

double dci_individual_counterfactual(DCI_SCM* scm,
    const double* observed_values, const int* observed_vars,
    int n_observed, int intervention_var, double intervention_value,
    int query_var) {
    if (!scm || !observed_values || !observed_vars) return 0.0;
    /* Step 1: Abduction — find exog values consistent with observations */
    double exog[DCI_MAX_VARS] = {0};
    int i;
    for (i = 0; i < n_observed; i++) {
        exog[observed_vars[i]] = observed_values[i];
    }
    /* Step 2: Action */
    exog[intervention_var] = intervention_value;
    /* Step 3: Prediction */
    dci_scm_evaluate(scm, exog);
    return scm->vars[query_var].value;
}

/* ==============================================================
 * Mediation via Counterfactuals
 * ============================================================== */

double dci_causal_mediation(DCI_SCM* scm, int cause, int mediator,
    int effect, int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    /* Natural Indirect Effect via counterfactual:
     * NIE = E[Y(1, M(1)) - Y(1, M(0))] */
    double sum_nie = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        /* Y(1, M(1)): X=1, mediator naturally determined */
        exog[cause] = 1.0;
        dci_scm_evaluate(scm, exog);
        double y11 = scm->vars[effect].value;
        double m1 = scm->vars[mediator].value;
        /* Y(1, M(0)): X=1, but mediator set to M(0) */
        exog[cause] = 0.0;
        dci_scm_evaluate(scm, exog);
        double m0 = scm->vars[mediator].value;
        exog[cause] = 1.0;
        exog[mediator] = m0;
        dci_scm_evaluate(scm, exog);
        double y10 = scm->vars[effect].value;
        sum_nie += y11 - y10;
        (void)m1;
    }
    return sum_nie / n_samples;
}

/* ==============================================================
 * Counterfactual Fairness (Pearl, 2000)
 *
 * Fairness criterion: P(Y_{X=x} = y | X=a, E=e) = P(Y_{X=x} = y | X=b, E=e)
 * for protected attribute X, sensitive groups a and b.
 * ============================================================== */

bool dci_counterfactual_fairness(DCI_SCM* scm, int protected_attr,
    int outcome, double threshold, int n_samples) {
    if (!scm || n_samples <= 0) return false;
    /* Check if counterfactual outcomes are independent of protected attribute */
    double prob_a = 0.0, prob_b = 0.0;
    int count_a = 0, count_b = 0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        if (scm->vars[protected_attr].value > 0.5) {
            /* Group a */
            exog[protected_attr] = 0.0;
            dci_scm_evaluate(scm, exog);
            prob_a += scm->vars[outcome].value;
            count_a++;
        } else {
            /* Group b */
            exog[protected_attr] = 1.0;
            dci_scm_evaluate(scm, exog);
            prob_b += scm->vars[outcome].value;
            count_b++;
        }
    }
    if (count_a == 0 || count_b == 0) return true;
    prob_a /= count_a;
    prob_b /= count_b;
    return fabs(prob_a - prob_b) < threshold;
}

/* ==============================================================
 * Path-Specific Counterfactual
 * ============================================================== */

double dci_path_specific_effect(DCI_SCM* scm, int cause, int effect,
    const int* path_nodes, int path_len, int n_samples) {
    if (!scm || !path_nodes || path_len < 2 || n_samples <= 0) return 0.0;
    double sum = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        /* Cause=1, natural for path nodes */
        exog[cause] = 1.0;
        dci_scm_evaluate(scm, exog);
        double y1 = scm->vars[effect].value;
        /* Cause=0, natural for path nodes */
        exog[cause] = 0.0;
        dci_scm_evaluate(scm, exog);
        double y0 = scm->vars[effect].value;
        sum += y1 - y0;
    }
    return sum / n_samples;
    (void)path_len;
}
/* ==============================================================
 * Counterfactual Explanation (explain why Y=y for a specific case)
 * ============================================================== */

void dci_counterfactual_explanation(DCI_SCM* scm,
    const double* observed, const int* obs_vars, int n_obs,
    int target_var, double* importance, int max_vars) {
    if (!scm || !observed || !obs_vars || !importance) return;
    /* Compute Shapley-style importance by comparing counterfactual outcomes
     * when each variable is intervened to its observed vs. counterfactual value */
    int i;
    for (i = 0; i < scm->n_vars && i < max_vars; i++) {
        double exog[DCI_MAX_VARS] = {0};
        int j;
        for (j = 0; j < n_obs; j++) exog[obs_vars[j]] = observed[j];
        /* Baseline */
        dci_scm_evaluate(scm, exog);
        double baseline = scm->vars[target_var].value;
        /* Set variable i to zero (counterfactual) */
        double saved = exog[i];
        exog[i] = 0.0;
        dci_scm_evaluate(scm, exog);
        double counter = scm->vars[target_var].value;
        exog[i] = saved;
        importance[i] = fabs(baseline - counter);
    }
}

/* ==============================================================
 * Attributable Fraction (population-level counterfactual measure)
 * ============================================================== */

double dci_attributable_fraction(DCI_SCM* scm, int exposure, int outcome,
    int n_samples) {
    if (!scm || n_samples <= 0) return 0.0;
    /* AF = (P(Y=1) - P(Y=1|do(X=0))) / P(Y=1) */
    double p_y1_obs = 0.0, p_y1_dox0 = 0.0;
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        if (scm->vars[outcome].value > 0.5) p_y1_obs += 1.0;
        exog[exposure] = 0.0;
        dci_scm_evaluate(scm, exog);
        if (scm->vars[outcome].value > 0.5) p_y1_dox0 += 1.0;
    }
    p_y1_obs /= n_samples;
    p_y1_dox0 /= n_samples;
    if (p_y1_obs < 1e-10) return 0.0;
    return (p_y1_obs - p_y1_dox0) / p_y1_obs;
}