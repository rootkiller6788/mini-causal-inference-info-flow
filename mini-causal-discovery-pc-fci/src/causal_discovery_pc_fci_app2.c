/* causal_discovery_pc_fci_app2.c -- L7: Supply chain causal discovery
 *
 * Application: Causal discovery in supply chain management.
 *
 * Uses the PC and FCI algorithms to discover causal relationships
 * among supply chain variables: supplier lead times, demand signals,
 * inventory levels, production rates, and delivery performance.
 *
 * In supply chains, causal discovery reveals hidden dependencies
 * that traditional correlation analysis misses - crucial for
 * identifying bottlenecks, propagating disruptions, and building
 * resilient supplier networks.
 *
 * Keywords: supplier, supply chain, causal discovery, PC, FCI
 */

#include "cdf_core.h"
#include "cdf_citest.h"
#include "cdf_graph.h"
#include "cdf_orientation.h"
#include "cdf_pc.h"
#include "cdf_fci.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N_SAMPLES 600
#define N_VARS    5

enum { DEMAND = 0, LEAD_TIME = 1, INVENTORY = 2, PRODUCTION = 3, DELIVERY = 4 };

/*
 * Generate synthetic supply chain data.
 *
 * Causal DAG:
 *   DEMAND     -> PRODUCTION    (demand drives production planning)
 *   LEAD_TIME  -> DELIVERY      (supplier lead time affects delivery speed)
 *   DEMAND     -> INVENTORY     (demand depletes inventory)
 *   PRODUCTION -> INVENTORY     (production replenishes inventory)
 *   INVENTORY  -> DELIVERY      (inventory availability affects delivery)
 *   LEAD_TIME  -> INVENTORY     (long lead times require higher buffer)
 *
 * This DAG captures key dynamics in a multi-echelon supply chain.
 */
static double* generate_supply_chain_data(int N)
{
    double *data = (double*)calloc((size_t)N_VARS * N, sizeof(double));
    if (!data) return NULL;

    for (int t = 0; t < N; t++) {
        double e_demand     = ((double)((t * 11  + 23)  % 991)  / 991.0  - 0.5) * 0.4;
        double e_lead       = ((double)((t * 17  + 41)  % 997)  / 997.0  - 0.5) * 0.3;
        double e_inv        = ((double)((t * 29  + 59)  % 1009) / 1009.0 - 0.5) * 0.25;
        double e_prod       = ((double)((t * 37  + 73)  % 1013) / 1013.0 - 0.5) * 0.2;
        double e_deliv      = ((double)((t * 43  + 97)  % 1019) / 1019.0 - 0.5) * 0.15;

        double demand = 1.0 + e_demand;
        data[DEMAND * N + t] = demand;

        double lead = 1.0 + e_lead;
        data[LEAD_TIME * N + t] = lead;

        double production = 0.7 * demand + e_prod;
        data[PRODUCTION * N + t] = production;

        double inventory = 0.5 * production - 0.3 * demand + 0.2 * lead + e_inv;
        data[INVENTORY * N + t] = inventory;

        double delivery = -0.3 * lead + 0.4 * inventory + e_deliv;
        data[DELIVERY * N + t] = delivery;
    }

    return data;
}

/* Print supply chain variable correlation matrix */
static void print_sc_correlations(const CdfDataset *ds)
{
    int p = ds->p;
    double *corr = (double*)malloc((size_t)p * p * sizeof(double));
    if (!corr) return;
    cdf_citest_corr_matrix(ds, corr);
    printf("\n  Supply chain correlation matrix:\n");
    printf("           Demand  LeadTm  Invtry  Prodct  Delivr\n");
    const char *names[] = {"Demand","LeadTm","Invtry","Prodct","Delivr"};
    for (int i = 0; i < p; i++) {
        printf("  %-6s", names[i]);
        for (int j = 0; j < p; j++)
            printf(" %+6.3f", corr[i * p + j]);
        printf("\n");
    }
    printf("  Note: Correlations do NOT imply causation.\n");
    printf("        Causal discovery separates true causes from associations.\n");
    free(corr);
}

/* Analyze discovered edges in supply chain context */
static void analyze_sc_edges(const CdfGraph *g, const char *alg_name)
{
    int p = g->p;
    const char *vars[] = {"Demand","LeadTm","Invtry","Prodct","Delivr"};
    printf("  %s discovered edges:\n", alg_name);
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            CdfEdgeType e = g->edges[i * g->p + j];
            const char *sym = cdf_graph_edge_char(e);
            if (e != CDF_EDGE_NONE) {
                printf("    %-6s %s %-6s", vars[i], sym, vars[j]);
                if (i == DEMAND && j == PRODUCTION && e == CDF_EDGE_DIRECTED)
                    printf("  [Demand -> production planning]");
                else if (i == DEMAND && j == INVENTORY && e == CDF_EDGE_DIRECTED)
                    printf("  [Demand depletes inventory]");
                else if (i == PRODUCTION && j == INVENTORY && e == CDF_EDGE_DIRECTED)
                    printf("  [Production replenishes inventory]");
                else if (i == LEAD_TIME && j == DELIVERY && e == CDF_EDGE_DIRECTED)
                    printf("  [Supplier lead time -> delivery delay]");
                else if (i == LEAD_TIME && j == INVENTORY && e == CDF_EDGE_DIRECTED)
                    printf("  [Long lead times -> higher buffer inventory]");
                else if (i == INVENTORY && j == DELIVERY && e == CDF_EDGE_DIRECTED)
                    printf("  [Inventory availability -> delivery performance]");
                else if (e == CDF_EDGE_BIDIRECTED)
                    printf("  [Latent confounder detected - check supplier dependencies]");
                printf("\n");
            }
        }
    }
}

/* Evaluate graph quality metrics */
static void evaluate_graph_quality(const CdfGraph *g, const char *name)
{
    int p = g->p;
    int n_directed = 0, n_undirected = 0, n_bidirected = 0;
    int n_partial = 0, n_nondir = 0;
    int n_edges = 0;

    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            CdfEdgeType e = g->edges[i * g->p + j];
            if (e == CDF_EDGE_NONE) continue;
            n_edges++;
            switch (e) {
            case CDF_EDGE_DIRECTED:   n_directed++; break;
            case CDF_EDGE_UNDIRECTED: n_undirected++; break;
            case CDF_EDGE_BIDIRECTED: n_bidirected++; break;
            case CDF_EDGE_PARTIAL_I:
            case CDF_EDGE_PARTIAL_J:  n_partial++; break;
            case CDF_EDGE_NONDIR:     n_nondir++; break;
            default: break;
            }
        }
    }

    printf("  %s quality metrics:\n", name);
    printf("    Total edges:    %d\n", n_edges);
    printf("    Directed:       %d (%.0f%%  - confirmed causal arrows)\n",
           n_directed, n_edges > 0 ? 100.0 * n_directed / n_edges : 0.0);
    printf("    Undirected:     %d (%.0f%%  - Markov equivalence ambiguity)\n",
           n_undirected, n_edges > 0 ? 100.0 * n_undirected / n_edges : 0.0);
    printf("    Bidirected:     %d (%.0f%%  - latent confounders)\n",
           n_bidirected, n_edges > 0 ? 100.0 * n_bidirected / n_edges : 0.0);
    printf("    Partially dir:  %d (FCI-specific, tail uncertain)\n", n_partial);
    printf("    Nondirected:    %d (FCI-specific, both ends uncertain)\n", n_nondir);

    /* Sparsity metric */
    int max_possible = p * (p - 1) / 2;
    printf("    Sparsity:       %.0f%% (%d/%d max)\n",
           100.0 - 100.0 * n_edges / max_possible, n_edges, max_possible);
}

/* Conditional independence analysis for key supply chain queries */
static void supply_chain_ci_queries(const CdfDataset *ds)
{
    double alpha = 0.05;

    printf("\n  === Supply chain conditional independence queries ===\n");

    /* Q1: Demand _|_ Delivery | {Inventory, Production}? */
    /* In the true DAG: Demand -> Production -> Inventory -> Delivery,
     * so controlling for Inventory and Production should block the path. */
    int Z1[] = {INVENTORY, PRODUCTION};
    CdfCITestResult r1 = cdf_citest_fisher_z(ds, DEMAND, DELIVERY, Z1, 2, alpha);
    printf("  Q1: Demand _|_ Delivery | {Inv,Prod}?\n");
    printf("      rho_partial=%.4f, p=%.4f, %s\n",
           cdf_citest_partial_corr(ds, DEMAND, DELIVERY, Z1, 2),
           r1.p_value, r1.is_independent ? "INDEP" : "DEP");

    /* Q2: LeadTime _|_ Inventory | {Demand}? */
    /* Lead -> Inv directly, so NOT independent even after controlling Demand. */
    int Z2[] = {DEMAND};
    CdfCITestResult r2 = cdf_citest_fisher_z(ds, LEAD_TIME, INVENTORY, Z2, 1, alpha);
    printf("  Q2: LeadTm _|_ Inventory | {Demand}?\n");
    printf("      rho_partial=%.4f, p=%.4f, %s (expect DEP - direct edge)\n",
           cdf_citest_partial_corr(ds, LEAD_TIME, INVENTORY, Z2, 1),
           r2.p_value, r2.is_independent ? "INDEP" : "DEP");

    /* Q3: LeadTime _|_ Production | {}? */
    /* These are marginally independent (no edge, no confounding path). */
    CdfCITestResult r3 = cdf_citest_fisher_z(ds, LEAD_TIME, PRODUCTION, NULL, 0, alpha);
    printf("  Q3: LeadTm _|_ Production | {}?\n");
    printf("      rho=%.4f, p=%.4f, %s (expect INDEP - no causal link)\n",
           cdf_citest_correlation(ds, LEAD_TIME, PRODUCTION),
           r3.p_value, r3.is_independent ? "INDEP" : "DEP");
}

int main(void)
{
    printf("==========================================\n");
    printf("  L7: Causal Discovery - Supply Chain     \n");
    printf("==========================================\n\n");

    printf("[1] Generating supply chain data...\n");
    printf("    Variables: Demand, LeadTime, Inventory, Production, Delivery\n");
    printf("    Samples: %d\n\n", N_SAMPLES);

    double *raw = generate_supply_chain_data(N_SAMPLES);
    if (!raw) { printf("ERROR: Data gen failed\n"); return 1; }

    CdfDataset *ds = cdf_dataset_create(raw, N_SAMPLES, N_VARS);
    free(raw);
    if (!ds) { printf("ERROR: Dataset failed\n"); return 1; }

    print_sc_correlations(ds);

    /* PC algorithm */
    printf("\n[2] Running PC algorithm...\n");
    CdfPCConfig pc_cfg = cdf_pc_config_default();
    pc_cfg.alpha = 0.05;

    CdfPCResult *pc_res = cdf_pc_run(ds, &pc_cfg);
    if (!pc_res) { printf("ERROR: PC failed\n"); cdf_dataset_free(ds); return 1; }
    cdf_pc_print_result(pc_res);
    analyze_sc_edges(pc_res->graph, "PC");
    evaluate_graph_quality(pc_res->graph, "PC CPDAG");

    /* FCI algorithm */
    printf("\n[3] Running FCI algorithm...\n");
    CdfFCIConfig fci_cfg = cdf_fci_config_default();
    fci_cfg.alpha = 0.05;
    fci_cfg.max_path_length = 3;

    CdfFCIResult *fci_res = cdf_fci_run(ds, &fci_cfg);
    if (!fci_res) {
        printf("ERROR: FCI failed\n");
        cdf_pc_result_free(pc_res); cdf_dataset_free(ds);
        return 1;
    }
    cdf_fci_print_result(fci_res);
    analyze_sc_edges(fci_res->graph, "FCI");
    evaluate_graph_quality(fci_res->graph, "FCI PAG");

    /* CI queries for supply chain decisions */
    supply_chain_ci_queries(ds);

    /* Compare PC and FCI results */
    printf("\n[4] Comparison:\n");
    int diff = cdf_fci_compare_pags(pc_res->graph, fci_res->graph);
    printf("    Edge differences: %d\n", diff);
    printf("    (Supply chains often have latent variables from\n");
    printf("     unobserved market factors; FCI captures these.)\n");

    /* PC-FCI agreement matrix */
    printf("\n    Edge-by-edge agreement (PC row vs FCI col):\n");
    printf("    Pair      PC        FCI       Match?\n");
    printf("    --------  --------  --------  ------\n");
    const char *pairs[] = {"Dem-Lead","Dem-Inv","Dem-Prod","Dem-Delv",
                           "Lead-Inv","Lead-Prod","Lead-Delv",
                           "Inv-Prod","Inv-Delv","Prod-Delv"};
    int u[] = {0,0,0,0, 1,1,1, 2,2, 3};
    int v[] = {1,2,3,4, 2,3,4, 3,4, 4};
    for (int k = 0; k < 10; k++) {
        CdfEdgeType epc = pc_res->graph->edges[u[k] * N_VARS + v[k]];
        CdfEdgeType efci = fci_res->graph->edges[u[k] * N_VARS + v[k]];
        printf("    %-9s %-9s %-9s %s\n",
               pairs[k],
               cdf_graph_edge_char(epc),
               cdf_graph_edge_char(efci),
               (epc == efci) ? "YES" : "DIFF");
    }

    /* Cleanup */
    cdf_pc_result_free(pc_res);
    cdf_fci_result_free(fci_res);
    cdf_dataset_free(ds);

    printf("\n==========================================\n");
    printf("  Application complete.\n");
    printf("==========================================\n");
    return 0;
}