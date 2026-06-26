#include "dci_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* SCM create */
DCI_SCM* dci_scm_create(const char* name) {
    DCI_SCM* scm = (DCI_SCM*)calloc(1, sizeof(DCI_SCM));
    if (!scm) return NULL;
    snprintf(scm->name, sizeof(scm->name), "%s", name ? name : "unnamed");
    return scm;
}

void dci_scm_free(DCI_SCM* scm) { free(scm); }

/*
 * Add variable to SCM with structural equation.
 * Variables are assigned IDs 0, 1, 2, ... in order of addition.
 */
int dci_scm_add_variable(DCI_SCM* scm, const char* name, DCI_VarType type,
    int* parent_ids, int n_parents, DCI_StructuralEq eq, bool exogenous) {
    if (!scm || scm->n_vars >= DCI_MAX_VARS) return -1;
    int id = scm->n_vars;
    DCI_Variable* v = &scm->vars[id];
    v->id = id;
    strncpy(v->name, name, 63);
    v->type = type;
    v->equation = eq;
    v->is_exogenous = exogenous;
    v->n_parents = n_parents;
    int i;
    for (i = 0; i < n_parents && i < DCI_MAX_PARENTS; i++) {
        v->parent_ids[i] = parent_ids[i];
        /* Update child list in parent */
        int pid = parent_ids[i];
        if (pid >= 0 && pid < id) {
            DCI_Variable* p = &scm->vars[pid];
            if (p->n_children < DCI_MAX_VARS) {
                p->child_ids[p->n_children++] = id;
            }
        }
    }
    scm->n_vars++;
    return id;
}

void dci_scm_add_edge(DCI_SCM* scm, int from_id, int to_id) {
    if (!scm) return;
    if (from_id < 0 || from_id >= scm->n_vars) return;
    if (to_id < 0 || to_id >= scm->n_vars) return;
    DCI_Variable* from = &scm->vars[from_id];
    if (from->n_children < DCI_MAX_VARS) {
        from->child_ids[from->n_children++] = to_id;
    }
    DCI_Variable* to = &scm->vars[to_id];
    if (to->n_parents < DCI_MAX_PARENTS) {
        to->parent_ids[to->n_parents++] = from_id;
    }
}

/*
 * Check if SCM graph is a DAG via DFS cycle detection.
 */
bool dci_scm_is_dag(DCI_SCM* scm) {
    if (!scm) return false;
    int n = scm->n_vars;
    int visited[DCI_MAX_VARS] = {0};
    int rec_stack[DCI_MAX_VARS] = {0};
    int stack[DCI_MAX_VARS];
    int top;
    int i;
    for (i = 0; i < n; i++) {
        if (visited[i]) continue;
        top = 0;
        stack[top++] = i;
        while (top > 0) {
            int node = stack[top - 1];
            if (!visited[node]) {
                visited[node] = 1;
                rec_stack[node] = 1;
                DCI_Variable* v = &scm->vars[node];
                int j;
                for (j = 0; j < v->n_children; j++) {
                    int child = v->child_ids[j];
                    if (!visited[child]) stack[top++] = child;
                    else if (rec_stack[child]) return false;
                }
            } else {
                top--;
                rec_stack[node] = 0;
            }
        }
    }
    return true;
}

/*
 * Evaluate SCM given exogenous values.
 * Topological order assumed; evaluates each variable in sequence.
 */
double dci_scm_evaluate(DCI_SCM* scm, const double* exog_vals) {
    if (!scm) return 0.0;
    int i;
    for (i = 0; i < scm->n_vars; i++) {
        DCI_Variable* v = &scm->vars[i];
        if (v->is_exogenous && exog_vals) {
            v->value = exog_vals[i];
        } else if (v->equation) {
            double parent_vals[DCI_MAX_PARENTS];
            int j;
            for (j = 0; j < v->n_parents; j++) {
                int pid = v->parent_ids[j];
                parent_vals[j] = scm->vars[pid].value;
            }
            double u = exog_vals ? exog_vals[i] : 0.0;
            v->value = v->equation(parent_vals, &u);
        }
    }
    return scm->vars[scm->n_vars - 1].value;
}

void dci_scm_print(DCI_SCM* scm) {
    if (!scm) return;
    printf("SCM: %s (%d variables)\n", scm->name, scm->n_vars);
    int i;
    for (i = 0; i < scm->n_vars; i++) {
        DCI_Variable* v = &scm->vars[i];
        printf("  V%d [%s] type=%d parents={", v->id, v->name, v->type);
        int j;
        for (j = 0; j < v->n_parents; j++) {
            printf("%d%s", v->parent_ids[j], j < v->n_parents-1 ? "," : "");
        }
        printf("} children={");
        for (j = 0; j < v->n_children; j++) {
            printf("%d%s", v->child_ids[j], j < v->n_children-1 ? "," : "");
        }
        printf("} exog=%s\n", v->is_exogenous ? "yes" : "no");
    }
}

DCI_Variable* dci_var_by_name(DCI_SCM* scm, const char* name) {
    if (!scm || !name) return NULL;
    int i;
    for (i = 0; i < scm->n_vars; i++) {
        if (strcmp(scm->vars[i].name, name) == 0) return &scm->vars[i];
    }
    return NULL;
}

DCI_Variable* dci_var_by_id(DCI_SCM* scm, int id) {
    if (!scm || id < 0 || id >= scm->n_vars) return NULL;
    return &scm->vars[id];
}

int dci_var_count_parents(DCI_SCM* scm, int var_id) {
    DCI_Variable* v = dci_var_by_id(scm, var_id);
    return v ? v->n_parents : -1;
}

int dci_var_count_children(DCI_SCM* scm, int var_id) {
    DCI_Variable* v = dci_var_by_id(scm, var_id);
    return v ? v->n_children : -1;
}

/* ---- Graph Operations ---- */

DCI_Graph* dci_graph_create(int n_nodes) {
    DCI_Graph* g = (DCI_Graph*)calloc(1, sizeof(DCI_Graph));
    if (!g) return NULL;
    g->n_nodes = n_nodes;
    return g;
}

void dci_graph_free(DCI_Graph* g) { free(g); }

void dci_graph_add_edge(DCI_Graph* g, int from, int to, double weight) {
    if (!g || from < 0 || from >= g->n_nodes || to < 0 || to >= g->n_nodes)
        return;
    g->adjacency[from][to] = weight;
}

void dci_graph_from_scm(DCI_SCM* scm, DCI_Graph* g) {
    if (!scm || !g) return;
    g->n_nodes = scm->n_vars;
    int i;
    for (i = 0; i < scm->n_vars; i++) {
        strncpy(g->names[i], scm->vars[i].name, 63);
        int j;
        for (j = 0; j < scm->vars[i].n_children; j++) {
            int child = scm->vars[i].child_ids[j];
            g->adjacency[i][child] = 1.0;
        }
    }
}

void dci_graph_print(DCI_Graph* g) {
    if (!g) return;
    printf("Graph: %d nodes\n", g->n_nodes);
    int i, j;
    for (i = 0; i < g->n_nodes; i++) {
        printf("  %d -> {", i);
        int first = 1;
        for (j = 0; j < g->n_nodes; j++) {
            if (g->adjacency[i][j] > 0.5) {
                printf("%s%d", first ? "" : ",", j);
                first = 0;
            }
        }
        printf("}\n");
    }
}

/* ==============================================================
 * SCM Simulation
 * ============================================================== */

void dci_scm_simulate(DCI_SCM* scm, int n_samples, double** results) {
    if (!scm || !results || n_samples <= 0) return;
    int i, j;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        dci_scm_evaluate(scm, exog);
        for (j = 0; j < scm->n_vars; j++) {
            if (results[j]) results[j][i] = scm->vars[j].value;
        }
    }
}

/* ==============================================================
 * SCM from adjacency matrix (linear structural equations)
 * ============================================================== */

static double linear_eq(const double* parents, const double* u) {
    /* Linear: v = Σ β_i * pa_i + u */
    double sum = u ? *u : 0.0;
    return parents ? sum + parents[0] : sum;
}

DCI_SCM* dci_scm_from_adjacency(const double* A, int n_vars,
    const char** var_names) {
    if (!A || n_vars <= 0 || n_vars > DCI_MAX_VARS) return NULL;
    DCI_SCM* scm = dci_scm_create("from_adjacency");
    int i, j;
    /* Create all variables first */
    for (i = 0; i < n_vars; i++) {
        char name[64];
        snprintf(name, 63, "%s", var_names ? var_names[i] : "V");
        dci_scm_add_variable(scm, name, DCI_CONTINUOUS,
            NULL, 0, linear_eq, false);
    }
    /* Add edges */
    for (i = 0; i < n_vars; i++) {
        for (j = 0; j < n_vars; j++) {
            if (A[i * n_vars + j] != 0.0) {
                dci_scm_add_edge(scm, i, j);
            }
        }
    }
    return scm;
}

/* ==============================================================
 * Structural equation helpers
 * ============================================================== */

double dci_linear_equation(const double* parents, int n_parents,
    const double* coeffs, double intercept) {
    double sum = intercept;
    int i;
    for (i = 0; i < n_parents; i++) sum += coeffs[i] * parents[i];
    return sum;
}

double dci_logistic_equation(const double* parents, int n_parents,
    const double* coeffs, double intercept) {
    double z = dci_linear_equation(parents, n_parents, coeffs, intercept);
    return 1.0 / (1.0 + exp(-z));
}

/* ==============================================================
 * Graph comparison (SHD - Structural Hamming Distance)
 * ============================================================== */

int dci_graph_shd(DCI_Graph* g1, DCI_Graph* g2) {
    if (!g1 || !g2) return -1;
    int n = g1->n_nodes < g2->n_nodes ? g1->n_nodes : g2->n_nodes;
    int shd = 0;
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            bool e1 = g1->adjacency[i][j] > 0.5;
            bool e2 = g2->adjacency[i][j] > 0.5;
            if (e1 != e2) shd++;
        }
    }
    return shd;
}

/* ==============================================================
 * SCM intervention simulation (batch)
 * ============================================================== */

typedef struct {
    double* y1_samples;
    double* y0_samples;
    int n_samples;
    double mean_y1;
    double mean_y0;
} DCI_SimResult;

DCI_SimResult dci_simulate_intervention(DCI_SCM* scm,
    int cause, int effect, int n_samples) {
    DCI_SimResult r;
    memset(&r, 0, sizeof(r));
    if (!scm || n_samples <= 0) return r;
    r.n_samples = n_samples;
    r.y1_samples = (double*)malloc(n_samples * sizeof(double));
    r.y0_samples = (double*)malloc(n_samples * sizeof(double));
    int i;
    for (i = 0; i < n_samples; i++) {
        double exog[DCI_MAX_VARS];
        int j;
        for (j = 0; j < scm->n_vars; j++) exog[j] = (double)rand() / RAND_MAX;
        exog[cause] = 1.0;
        dci_scm_evaluate(scm, exog);
        r.y1_samples[i] = scm->vars[effect].value;
        r.mean_y1 += r.y1_samples[i];
        exog[cause] = 0.0;
        dci_scm_evaluate(scm, exog);
        r.y0_samples[i] = scm->vars[effect].value;
        r.mean_y0 += r.y0_samples[i];
    }
    r.mean_y1 /= n_samples;
    r.mean_y0 /= n_samples;
    return r;
}
/* ==============================================================
 * SCM Query: Evaluate expression in interventional distribution
 * ============================================================== */

double dci_scm_query(DCI_SCM* scm, int query_var, const int* interventions,
    const double* intervention_values, int n_interventions) {
    if (!scm) return 0.0;
    double exog[DCI_MAX_VARS] = {0};
    int i;
    for (i = 0; i < n_interventions; i++) {
        exog[interventions[i]] = intervention_values[i];
    }
    dci_scm_evaluate(scm, exog);
    return scm->vars[query_var].value;
}

/* ==============================================================
 * SCM Copy (deep copy)
 * ============================================================== */

DCI_SCM* dci_scm_copy(DCI_SCM* scm) {
    if (!scm) return NULL;
    DCI_SCM* copy = dci_scm_create(scm->name);
    int i;
    for (i = 0; i < scm->n_vars; i++) {
        dci_scm_add_variable(copy, scm->vars[i].name, scm->vars[i].type,
            scm->vars[i].parent_ids, scm->vars[i].n_parents,
            scm->vars[i].equation, scm->vars[i].is_exogenous);
        copy->vars[i].value = scm->vars[i].value;
    }
    return copy;
}

/* ==============================================================
 * Standard SCM Examples
 * ============================================================== */

/*
 * Simpsons Paradox SCM:
 *   Z (gender) → X (treatment), Z → Y (outcome), X → Y
 */
DCI_SCM* dci_scm_simpsons_paradox(void) {
    DCI_SCM* scm = dci_scm_create("simpsons_paradox");
    int empty[] = {0};
    /* Z: gender (binary, exogenous) */
    dci_scm_add_variable(scm, "Z", DCI_BINARY, NULL, 0, NULL, true);
    /* X: treatment (binary, depends on Z) */
    dci_scm_add_variable(scm, "X", DCI_BINARY, empty, 1, linear_eq, false);
    /* Y: outcome (depends on X and Z) */
    int parents_y[] = {0, 1}; /* Z and X */
    dci_scm_add_variable(scm, "Y", DCI_CONTINUOUS, parents_y, 2, linear_eq, false);
    return scm;
}

/*
 * Instrumental Variable SCM:
 *   Z → X, X → Y, U → X, U → Y (U is unobserved confounder)
 */
DCI_SCM* dci_scm_iv(int n_samples) {
    (void)n_samples;
    DCI_SCM* scm = dci_scm_create("iv_model");
    /* Z: instrument (binary, exogenous) */
    dci_scm_add_variable(scm, "Z", DCI_BINARY, NULL, 0, NULL, true);
    /* U: unobserved confounder */
    dci_scm_add_variable(scm, "U", DCI_CONTINUOUS, NULL, 0, NULL, true);
    /* X: treatment (depends on Z and U) */
    int px[] = {0, 1}; /* Z, U */
    dci_scm_add_variable(scm, "X", DCI_CONTINUOUS, px, 2, linear_eq, false);
    /* Y: outcome (depends on X and U) */
    int py[] = {2, 1}; /* X, U */
    dci_scm_add_variable(scm, "Y", DCI_CONTINUOUS, py, 2, linear_eq, false);
    return scm;
}
/* ==============================================================
 * Mediation SCM (Baron-Kenny, 1986)
 *   X → M → Y, X → Y (direct + indirect paths)
 * ============================================================== */

DCI_SCM* dci_scm_mediation(void) {
    DCI_SCM* scm = dci_scm_create("mediation");
    /* X: cause (binary, exogenous) */
    dci_scm_add_variable(scm, "X", DCI_BINARY, NULL, 0, NULL, true);
    /* M: mediator (depends on X) */
    int pm[] = {0}; /* X */
    dci_scm_add_variable(scm, "M", DCI_CONTINUOUS, pm, 1, linear_eq, false);
    /* Y: outcome (depends on X and M) */
    int py[] = {0, 1}; /* X, M */
    dci_scm_add_variable(scm, "Y", DCI_CONTINUOUS, py, 2, linear_eq, false);
    return scm;
}

/* ==============================================================
 * M-bias SCM (conditioning on collider creates spurious association)
 *   U1 → X, U1 → Z, U2 → Y, U2 → Z
 * ============================================================== */

DCI_SCM* dci_scm_m_bias(void) {
    DCI_SCM* scm = dci_scm_create("m_bias");
    dci_scm_add_variable(scm, "U1", DCI_CONTINUOUS, NULL, 0, NULL, true);
    dci_scm_add_variable(scm, "U2", DCI_CONTINUOUS, NULL, 0, NULL, true);
    int px[] = {0}; /* U1 → X */
    dci_scm_add_variable(scm, "X", DCI_CONTINUOUS, px, 1, linear_eq, false);
    int py[] = {1}; /* U2 → Y */
    dci_scm_add_variable(scm, "Y", DCI_CONTINUOUS, py, 1, linear_eq, false);
    int pz[] = {0, 1}; /* U1 → Z ← U2 (collider) */
    dci_scm_add_variable(scm, "Z", DCI_CONTINUOUS, pz, 2, linear_eq, false);
    return scm;
}

/* ==============================================================
 * Front-door SCM
 *   U → X, X → M → Y, U → Y (U is unobserved confounder)
 * ============================================================== */

DCI_SCM* dci_scm_frontdoor(void) {
    DCI_SCM* scm = dci_scm_create("frontdoor");
    dci_scm_add_variable(scm, "U", DCI_CONTINUOUS, NULL, 0, NULL, true);
    int px[] = {0}; /* U → X */
    dci_scm_add_variable(scm, "X", DCI_CONTINUOUS, px, 1, linear_eq, false);
    int pm[] = {1}; /* X → M */
    dci_scm_add_variable(scm, "M", DCI_CONTINUOUS, pm, 1, linear_eq, false);
    int py[] = {2, 0}; /* M → Y, U → Y */
    dci_scm_add_variable(scm, "Y", DCI_CONTINUOUS, py, 2, linear_eq, false);
    return scm;
}
/* Quick utility: get SCM variable count */
int dci_scm_var_count(const DCI_SCM* scm) { return scm ? scm->n_vars : 0; }
/* Get variable value */
double dci_scm_get_value(const DCI_SCM* scm, int var_id) {
    if (!scm || var_id < 0 || var_id >= scm->n_vars) return 0.0;
    return scm->vars[var_id].value;
}