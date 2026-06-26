#ifndef CDF_ORIENTATION_H
#define CDF_ORIENTATION_H

/*
 * cdf_orientation.h — Edge Orientation Rules (PC & FCI)
 *
 * PC algorithm orientation rules (Meek, 1995):
 *   R0: Orient v-structures (unshielded colliders)
 *   R1: Avoid new v-structures (acyclicity)
 *   R2: Avoid cycles (transitive closure)
 *   R3: Avoid new v-structures by transitivity
 *   R4: Discriminating paths
 *
 * FCI orientation rules (Zhang, 2008):
 *   R0-R4: Same as PC but applied to PAG
 *   R5-R7: Additional rules for latent variable structures
 *   R8-R10: Rules for selection bias
 *
 * The rules are applied exhaustively (fixed-point iteration)
 * until no more edges can be oriented.
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "cdf_core.h"

/* ── Rule Application Result ─────────────────────────────────────── */

typedef struct {
    int  rules_applied[10];  /* count of each rule application */
    int  total_oriented;     /* total number of edges oriented */
    bool changed;            /* true if any orientation occurred */
} CdfOrientResult;

/* ── PC Orientation Rules ────────────────────────────────────────── */

/**
 * Apply all PC orientation rules exhaustively (R0-R4) to a graph
 * containing a skeleton with v-structures already oriented.
 *
 * R1: If a→b—c with a and c not adjacent, orient b→c.
 *     (Otherwise b would be a collider on (a,b,c))
 *
 * R2: If a→b→c and a—c, orient a→c.
 *     (Preserve acyclicity: transitive closure)
 *
 * R3: If a—b→c, a—d→c, b and d not adjacent, and a—c, orient a→c.
 *     (Complex v-structure avoidance)
 *
 * R4: Discriminating path rule for definite orientation.
 *
 * @param g   graph with skeleton + v-structures oriented
 * @return    orientation result with rule counts
 */
CdfOrientResult cdf_orient_pc_rules(CdfGraph *g);

/** Apply only R1 (no new v-structures). */
int cdf_orient_rule_r1(CdfGraph *g);

/** Apply only R2 (acyclicity/transitive closure). */
int cdf_orient_rule_r2(CdfGraph *g);

/** Apply only R3. */
int cdf_orient_rule_r3(CdfGraph *g);

/** Apply only R4 (discriminating path). */
int cdf_orient_rule_r4(CdfGraph *g);

/* ── FCI Orientation Rules ───────────────────────────────────────── */

/**
 * Apply all FCI orientation rules exhaustively (R0-R10) to a PAG.
 *
 * FCI-specific rules:
 *   R5: If a ∘→ b — c and a, c not adjacent, orient b → c.
 *   R6: If a ∘→ b ∘→ c with a—c unoriented, orient b → c.
 *   R7: If a ∘→ b → c with a—c, orient a → c.
 *   R8: If a → b → c and a ∘→ c, orient a → c.
 *   R9: If a ∘→ c and the path a ∘→ b → c, orient a → c.
 *   R10: If a ∘→ c and a → b ← c, orient a → c.
 *
 * @param g   PAG (will be modified)
 * @return    orientation result
 */
CdfOrientResult cdf_orient_fci_rules(CdfGraph *g);

/** Apply FCI rules R5-R7. */
int cdf_orient_fci_rules_5_7(CdfGraph *g);

/** Apply FCI rules R8-R10. */
int cdf_orient_fci_rules_8_10(CdfGraph *g);

/* ── Discriminating Path ─────────────────────────────────────────── */

/**
 * Check if a path ⟨v_0, ..., v_k=w⟩ is a discriminating path for (u,v,w).
 *
 * A discriminating path π = ⟨x, q_1, ..., q_m, b, y⟩ is an uncovered
 * path where x is not adjacent to y and every vertex q_i is a collider
 * on π and a parent of y.
 *
 * @param g      graph
 * @param u,v,w  triple (u—v—w in skeleton, u,w not adjacent)
 * @param path   [len] node sequence
 * @param len    path length
 * @return       true if the path is discriminating
 */
bool cdf_orient_is_discriminating_path(const CdfGraph *g,
                                        int u, int v, int w,
                                        const int *path, int len);

/* ── Orientation Utilities ───────────────────────────────────────── */

/** Check if orienting u→v would create a directed cycle.
 *  Uses DFS from v back to u along directed edges. */
bool cdf_orient_would_create_cycle(const CdfGraph *g, int u, int v);

/** Orient edge u→v (also sets v←u appropriately). */
void cdf_orient_edge(CdfGraph *g, int u, int v, CdfEdgeType type);

/** Count unoriented edges in the graph. */
int cdf_orient_count_unoriented(const CdfGraph *g);

/** Check if the graph is maximally oriented (all orientable edges done). */
bool cdf_orient_is_maximal(const CdfGraph *g);

/** Print orientation result summary. */
void cdf_orient_print_result(const CdfOrientResult *res);

#endif /* CDF_ORIENTATION_H */