/* ==============================================================
 * Example 2: Do-Calculus Rules & Graph Surgery
 *
 * Demonstrates Pearl's three do-calculus rules (1995) and
 * graph surgery for causal effect identification:
 *
 * Rule 1: P(y|do(x),z,w) = P(y|do(x),w) if Y ⟂ Z | X,W in G_X̅
 * Rule 2: P(y|do(x),do(z),w) = P(y|do(x),z,w) if Y ⟂ Z | X,W in G_X̅Z̲
 * Rule 3: P(y|do(x),do(z),w) = P(y|do(x),w) if Y ⟂ Z | X,W in G_X̅Z(W)̅
 *
 * Uses IV model: Instrument Z → Treatment X → Outcome Y
 *                Confounder U → X, U → Y
 *
 * References:
 *   Pearl, J. (1995). Biometrika 82(4):669-710.
 *   Pearl, J. (2009). Causality, Ch. 3, §3.4.
 *   Shpitser & Pearl (2006). JMLR 9:1941-1967.
 * ============================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dci_core.h"
#include "dci_graph.h"
#include "dci_do.h"
#include "dci_backdoor.h"
#include "dci_effect.h"

int main(void) {
    printf("=== Example 2: Do-Calculus & Graph Surgery ===\n");
    printf("IV Model: Z(instrument) → X(treatment) → Y(outcome)\n");
    printf("          U(confounder) → X,  U → Y\n\n");

    /* 1. Create IV SCM */
    DCI_SCM* scm = dci_scm_iv(100);
    printf("--- SCM Structure ---\n");
    dci_scm_print(scm);

    /* 2. Original causal graph */
    DCI_Graph g;
    dci_graph_from_scm(scm, &g);
    printf("\n--- Original Causal Graph ---\n");
    dci_graph_print(&g);

    /* 3. Topological order */
    DCI_TopoOrder topo = dci_topological_sort(&g);
    printf("\nTopological order: ");
    for (int i = 0; i < topo.n_nodes; i++)
        printf("%d ", topo.order[i]);
    printf("(DAG: %s)\n", dci_graph_is_dag(&g) ? "YES" : "NO");

    /* 4. Graph surgery: mutilate X (remove all incoming edges) */
    printf("\n--- Graph Surgery: do(X=x) ---\n");
    printf("Removing all incoming edges to X (node 2)...\n");
    DCI_Graph g_mut;
    dci_graph_mutilate(&g, (int[]){2}, 1, &g_mut);
    printf("After mutilation (G_{X̅}):\n");
    dci_graph_print(&g_mut);

    /* 5. Verify d-separation after surgery */
    printf("\n--- d-Separation After Surgery ---\n");
    bool z_indep_y = dci_is_d_separated(&g_mut, 0, 3, NULL, 0);
    printf("Z d-separated from Y in G_X̅: %s\n",
        z_indep_y ? "YES" : "NO");

    /* 6. Apply do-calculus Rule 1 */
    printf("\n--- Do-Calculus Rules ---\n");
    DCI_Derivation r1 = dci_apply_rule1(scm, &g, 3, (int[]){2}, 1,
        (int[]){0}, 1, NULL, 0);
    printf("Rule 1 (obs insertion/deletion): %s\n",
        r1.success ? "APPLICABLE" : "not applicable");

    /* 7. Apply do-calculus Rule 2 */
    DCI_Derivation r2 = dci_apply_rule2(scm, &g, 3, (int[]){2}, 1,
        (int[]){2}, 1, NULL, 0);
    printf("Rule 2 (action↔obs exchange):  %s\n",
        r2.success ? "APPLICABLE" : "not applicable");

    /* 8. Apply do-calculus Rule 3 */
    DCI_Derivation r3 = dci_apply_rule3(scm, &g, 3, (int[]){2}, 1,
        (int[]){0}, 1, NULL, 0);
    printf("Rule 3 (action insertion/deletion): %s\n",
        r3.success ? "APPLICABLE" : "not applicable");

    /* 9. Full do-calculus derivation */
    printf("\n--- Full Derivation ---\n");
    bool identifiable = false;
    DCI_Derivation full_der = dci_do_calculus(scm,
        "P(y|do(x))", &identifiable);
    printf("Target: P(y|do(x))\n");
    printf("Identifiable: %s\n", identifiable ? "YES" : "NO");
    printf("Steps: %d\n", full_der.n_steps);
    dci_derivation_print(&full_der);

    /* 10. Identifiability verification */
    printf("\n--- Identifiability Check ---\n");
    char derivation[256];
    bool id_check = dci_identify_effect(scm, 2, 3, derivation, 255);
    printf("P(y|do(x)) identifiable: %s\n",
        id_check ? "YES" : "NO");
    printf("Derivation path: %s\n", derivation);

    /* 11. Truncated factorization */
    printf("\n--- Truncated Factorization ---\n");
    printf("P(v|do(x)) = ∏_{i∉X} P(v_i|pa_i)\n");
    double obs[] = {0.5, 0.3, 0.5, 0.7};
    double tf = dci_truncated_factorization(scm, (int[]){2}, (double[]){1.0},
        1, obs);
    printf("Truncated factorization value: %.4f\n", tf);

    /* 12. Compute causal effect via do-calculus */
    printf("\n--- Causal Effect via do(X=1) ---\n");
    DCI_CausalEffect ce = dci_compute_ate(scm, 2, 3, 2000);
    dci_effect_print(&ce);

    /* 13. Check if Z is a valid instrument */
    printf("\n--- Instrumental Variable Check ---\n");
    bool is_iv = dci_is_instrument(&g, 0, 2, 3);
    printf("Is Z a valid instrument for X→Y? %s\n",
        is_iv ? "YES" : "NO");

    dci_scm_free(scm);
    printf("\n=== Example 2 Complete ===\n");
    return 0;
}
