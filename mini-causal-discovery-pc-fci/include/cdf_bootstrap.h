#ifndef CDF_BOOTSTRAP_H
#define CDF_BOOTSTRAP_H

/*
 * cdf_bootstrap.h — Bootstrap-Based Edge Stability Assessment
 *
 * In finite-sample causal discovery, the learned graph can vary with
 * sampling noise. Bootstrap resampling (Efron, 1979) provides
 * non-parametric estimates of edge appearance frequency and
 * orientation stability.
 *
 * For constraint-based methods (PC, FCI), bootstrap-based stability
 * selection (Meinshausen & Bühlmann, 2010) controls false positives:
 *   - Resample data with replacement B times
 *   - Run PC/FCI on each bootstrap sample
 *   - For each edge (i,j): compute frequency of appearance
 *   - Select edges exceeding a stability threshold π_thr
 *
 * Reference: Friedman, Goldszmidt & Wyner (1999), "Data Analysis with
 *   Bayesian Networks: A Bootstrap Approach". UAI.
 * Reference: Scutari & Nagarajan (2013), "Identifying significant edges
 *   in graphical models of molecular networks". AI in Medicine.
 */

#include <stddef.h>
#include <stdbool.h>
#include "cdf_core.h"
#include "cdf_pc.h"
#include "cdf_fci.h"

/* ── Bootstrap Configuration ───────────────────────────────────────── */

typedef struct {
    int     n_bootstrap;      /* number of bootstrap resamples (>=50) */
    double  stability_thresh; /* edge stability threshold (e.g. 0.7) */
    int     random_seed;      /* random seed for reproducibility */
    bool    use_pc;           /* true=PC, false=FCI */
    bool    verbose;          /* print progress */
} CdfBootstrapConfig;

/* ── Bootstrap Edge Record ─────────────────────────────────────────── */

/** Per-edge stability statistics across bootstrap replicates. */
typedef struct {
    int     u;               /* source node */
    int     v;               /* target node */
    int     n_appearances;   /* how many times edge appeared */
    int     n_directed_uv;   /* times oriented u->v */
    int     n_directed_vu;   /* times oriented v->u */
    int     n_undirected;    /* times undirected */
    int     n_bidirected;    /* times bidirected (FCI only) */
    double  freq_stable;     /* appearance frequency = n_appearances / B */
    bool    is_stable;       /* above stability threshold */
} CdfBootstrapEdge;

/* ── Bootstrap Result ──────────────────────────────────────────────── */

typedef struct {
    CdfBootstrapEdge *edges;  /* [max_pairs] edge records */
    int              n_edges; /* number of edge records */
    int              n_bootstrap; /* B */
    int              p;       /* number of variables */
    double           elapsed_sec; /* runtime */
} CdfBootstrapResult;

/* ── Bootstrap API ─────────────────────────────────────────────────── */

/**
 * Run bootstrap edge stability analysis on a dataset.
 *
 * Algorithm:
 *   1. Run PC/FCI on the full dataset to get reference graph G_ref
 *   2. For b = 1..B:
 *      a. Sample N rows with replacement -> D_b
 *      b. Run PC/FCI on D_b -> graph G_b
 *      c. Record each edge (u,v,type) present in G_b
 *   3. Compute appearance frequencies f_{uv} = count_{uv} / B
 *   4. Filter edges: include only if f_{uv} >= stability_thresh
 *
 * Complexity: O(B * T) where T is PC/FCI runtime on N samples.
 *
 * @param ds      dataset (N samples x p variables)
 * @param config  bootstrap configuration
 * @return        bootstrap result (caller must free with
 *                cdf_bootstrap_result_free)
 */
CdfBootstrapResult* cdf_bootstrap_edges(const CdfDataset *ds,
                                         const CdfBootstrapConfig *config);

/** Free bootstrap result. */
void cdf_bootstrap_result_free(CdfBootstrapResult *res);

/** Create default bootstrap configuration. */
CdfBootstrapConfig cdf_bootstrap_config_default(void);

/** Print bootstrap result summary. */
void cdf_bootstrap_print(const CdfBootstrapResult *res);

/* ── Stability Selection ───────────────────────────────────────────── */

/**
 * Stability selection: determine the optimal threshold for a given
 * target FDR (false discovery rate) using the method of
 * Meinshausen & Buhlmann (2010).
 *
 * For a threshold pi in [0.5, 1.0], the expected number of false
 * positives is bounded by: E[V] <= 1/(2*pi-1) * q^2/p
 * where q is the expected number of selected variables.
 *
 * @param edges      edge records from bootstrap
 * @param n_edges    number of records
 * @param fdr_target desired false discovery rate
 * @return           recommended stability threshold
 */
double cdf_bootstrap_fdr_threshold(const CdfBootstrapEdge *edges,
                                    int n_edges, double fdr_target);

/**
 * Extract the stable graph: build a CdfGraph containing only edges
 * that passed the stability threshold. Edge type is determined by
 * majority vote across bootstrap replicates.
 *
 * @param edges    edge records
 * @param n_edges  number of records
 * @param p        number of variables
 * @param thresh   stability threshold
 * @return         stable CdfGraph (caller must cdf_graph_free)
 */
CdfGraph* cdf_bootstrap_stable_graph(const CdfBootstrapEdge *edges,
                                      int n_edges, int p, double thresh);

/**
 * Generate one bootstrap sample by resampling N rows with replacement
 * from the dataset.
 *
 * @param ds  original dataset
 * @return    new dataset (caller must cdf_dataset_free)
 */
CdfDataset* cdf_bootstrap_resample(const CdfDataset *ds);

#endif /* CDF_BOOTSTRAP_H */
