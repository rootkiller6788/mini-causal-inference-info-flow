/* causal_discovery_pc_fci_app1.c -- L7: Climate causal discovery
 *
 * Application: Causal discovery in climate science.
 *
 * Uses PC and FCI algorithms to discover causal relationships
 * among climate variables from observational data.
 * Demonstrates how constraint-based causal discovery can reveal
 * the causal structure driving climate dynamics.
 *
 * Keywords: climate, causal discovery, PC, FCI, time series
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

#define N_SAMPLES 1000
#define N_VARS    5

enum { TEMP = 0, CO2 = 1, SOLAR = 2, OCEAN = 3, ICE = 4 };

/* Generate synthetic climate data with known causal structure:
 *
 * Causal DAG:
 *   SOLAR -> TEMP          (solar forcing warms surface)
 *   CO2   -> TEMP          (greenhouse effect)
 *   TEMP  -> ICE           (warming reduces ice extent)
 *   TEMP  -> OCEAN         (surface warming transfers to ocean)
 *   OCEAN -> ICE           (ocean heat content affects ice)
 */
static double* generate_climate_data(int N)
{
    double *data = (double*)calloc((size_t)N_VARS * N, sizeof(double));
    if (!data) return NULL;

    for (int t = 0; t < N; t++) {
        double e_solar  = ((double)((t * 7   + 13)  % 997)  / 997.0  - 0.5) * 0.3;
        double e_co2    = ((double)((t * 19  + 37)  % 1009) / 1009.0 - 0.5) * 0.2;
        double e_temp   = ((double)((t * 31  + 61)  % 1013) / 1013.0 - 0.5) * 0.15;
        double e_ocean  = ((double)((t * 53  + 89)  % 1019) / 1019.0 - 0.5) * 0.1;
        double e_ice    = ((double)((t * 73  + 127) % 1021) / 1021.0 - 0.5) * 0.12;

        double solar = e_solar;
        data[SOLAR * N + t] = solar;

        double co2_trend = 0.0002 * (double)t;
        double co2 = 1.0 + co2_trend + e_co2;
        data[CO2 * N + t] = co2;

        double temp = 0.4 * solar + 0.5 * (co2 - 1.0) + e_temp;
        data[TEMP * N + t] = temp;

        double ocean = 0.6 * temp + e_ocean;
        data[OCEAN * N + t] = ocean;

        double ice = -0.35 * temp - 0.25 * ocean + e_ice;
        data[ICE * N + t] = ice;
    }

    return data;
}

/* Correlation matrix summary for climate variables */
static void print_correlation_summary(const CdfDataset *ds)
{
    int p = ds->p;
    double *corr = (double*)malloc((size_t)p * p * sizeof(double));
    if (!corr) return;
    cdf_citest_corr_matrix(ds, corr);
    printf("\n  Climate variable correlation matrix:\n");
    printf("          Temp    CO2    Solar  Ocean   Ice\n");
    const char *names[] = {"Temp", "CO2", "Solar", "Ocean", "Ice"};
    for (int i = 0; i < p; i++) {
        printf("  %-5s", names[i]);
        for (int j = 0; j < p; j++)
            printf(" %+6.3f", corr[i * p + j]);
        printf("\n");
    }
    free(corr);
}

/* Interpret edges in climate context */
static void analyze_climate_edges(const CdfGraph *g, const char *alg_name)
{
    int p = g->p;
    const char *vars[] = {"Temp", "CO2", "Solar", "Ocean", "Ice"};
    printf("  %s learned edges:\n", alg_name);
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            CdfEdgeType e = g->edges[i * g->p + j];
            const char *sym = cdf_graph_edge_char(e);
            if (e != CDF_EDGE_NONE) {
                printf("    %-5s %s %-5s", vars[i], sym, vars[j]);
                if (i == SOLAR && j == TEMP && e == CDF_EDGE_DIRECTED)
                    printf("  [Solar forcing -> surface warming]");
                else if (i == CO2 && j == TEMP && e == CDF_EDGE_DIRECTED)
                    printf("  [Greenhouse: CO2 -> warmth]");
                else if (i == TEMP && j == ICE && e == CDF_EDGE_DIRECTED)
                    printf("  [Warming -> ice melt]");
                else if (i == TEMP && j == OCEAN && e == CDF_EDGE_DIRECTED)
                    printf("  [Surface warming -> ocean heat]");
                else if (i == OCEAN && j == ICE && e == CDF_EDGE_DIRECTED)
                    printf("  [Ocean heat -> ice loss]");
                printf("\n");
            }
        }
    }
}

/* Compare PC and FCI outputs */
static void compare_pc_fci(const CdfGraph *g_pc, const CdfGraph *g_fci)
{
    if (!g_pc || !g_fci || g_pc->p != g_fci->p) return;
    int diff = cdf_fci_compare_pags(g_pc, g_fci);
    printf("\n  === PC vs FCI comparison ===\n");
    printf("  Edge-type differences: %d\n", diff);
    int p = g_pc->p;
    int pc_dir = 0, fci_dir = 0, pc_undir = 0, fci_undir = 0;
    int fci_bidir = 0, fci_nondir = 0;
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            switch (g_pc->edges[i * p + j]) {
            case CDF_EDGE_DIRECTED: pc_dir++; break;
            case CDF_EDGE_UNDIRECTED: pc_undir++; break;
            default: break;
            }
            switch (g_fci->edges[i * p + j]) {
            case CDF_EDGE_DIRECTED: fci_dir++; break;
            case CDF_EDGE_UNDIRECTED: fci_undir++; break;
            case CDF_EDGE_BIDIRECTED: fci_bidir++; break;
            case CDF_EDGE_NONDIR: fci_nondir++; break;
            default: break;
            }
        }
    }
    printf("  PC:  %d directed, %d undirected\n", pc_dir, pc_undir);
    printf("  FCI: %d directed, %d undirected, %d bidir, %d nondir\n",
           fci_dir, fci_undir, fci_bidir, fci_nondir);
    if (diff == 0)
        printf("  -> PC and FCI agree perfectly (no latent confounders)\n");
    else if (diff <= 2)
        printf("  -> Minor differences, consistent with finite-sample noise\n");
    else
        printf("  -> Notable differences: possible latent structure\n");
}

/* Alpha sensitivity sweep */
static void sensitivity_analysis(const CdfDataset *ds,
                                  CdfPCConfig *cfg,
                                  double alphas[], int n_alphas)
{
    printf("\n  === Sensitivity analysis (alpha sweep) ===\n");
    printf("  %-8s %-12s %-12s %-12s\n", "Alpha", "Edges", "V-struct", "Oriented");
    CdfPCConfig cfg_copy = *cfg;
    for (int k = 0; k < n_alphas; k++) {
        cfg_copy.alpha = alphas[k];
        CdfPCResult *res = cdf_pc_run(ds, &cfg_copy);
        if (!res) continue;
        int n_edges = 0;
        for (int i = 0; i < ds->p; i++)
            for (int j = i + 1; j < ds->p; j++)
                if (cdf_graph_has_edge(res->graph, i, j))
                    n_edges++;
        printf("  %-8.4f %-12d %-12d %-12d\n",
               alphas[k], n_edges, res->n_vstructs, res->n_oriented);
        cdf_pc_result_free(res);
    }
}

int main(void)
{
    printf("==========================================\n");
    printf("  L7: Causal Discovery - Climate Analysis \n");
    printf("==========================================\n\n");

    printf("[1] Generating climate time series...\n");
    printf("    Variables: Temp, CO2, Solar, Ocean, Ice\n");
    printf("    Samples: %d\n\n", N_SAMPLES);

    double *raw = generate_climate_data(N_SAMPLES);
    if (!raw) { printf("ERROR: Data gen failed\n"); return 1; }

    CdfDataset *ds = cdf_dataset_create(raw, N_SAMPLES, N_VARS);
    free(raw);
    if (!ds) { printf("ERROR: Dataset failed\n"); return 1; }

    print_correlation_summary(ds);

    printf("\n[2] Running PC algorithm...\n");
    CdfPCConfig pc_cfg = cdf_pc_config_default();
    pc_cfg.alpha = 0.05;
    pc_cfg.verbose = false;
    CdfPCResult *pc_res = cdf_pc_run(ds, &pc_cfg);
    if (!pc_res) { printf("ERROR: PC failed\n"); cdf_dataset_free(ds); return 1; }
    cdf_pc_print_result(pc_res);
    analyze_climate_edges(pc_res->graph, "PC");

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
    analyze_climate_edges(fci_res->graph, "FCI");

    compare_pc_fci(pc_res->graph, fci_res->graph);

    double alphas[] = {0.001, 0.01, 0.05, 0.10, 0.20};
    sensitivity_analysis(ds, &pc_cfg, alphas, 5);

    printf("\n[4] d-separation verification:\n");
    int sep_set[] = {TEMP};
    bool dsep = cdf_graph_d_separated(pc_res->graph, SOLAR, ICE, sep_set, 1);
    printf("    Solar _|_ Ice | {Temp}? %s\n", dsep ? "YES" : "NO");
    dsep = cdf_graph_d_separated(pc_res->graph, CO2, OCEAN, sep_set, 1);
    printf("    CO2   _|_ Ocean | {Temp}? %s\n", dsep ? "YES" : "NO");

    /* Structural Hamming Distance from true graph */
    CdfGraph *true_g = cdf_graph_create(N_VARS);
    cdf_graph_add_edge(true_g, SOLAR, TEMP,  CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(true_g, CO2,   TEMP,  CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(true_g, TEMP,  ICE,   CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(true_g, TEMP,  OCEAN, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(true_g, OCEAN, ICE,   CDF_EDGE_DIRECTED);
    int shd = cdf_pc_shd(pc_res->graph, true_g);
    printf("\n[5] SHD (PC vs True DAG): %d\n", shd);
    printf("    (Lower = better. SHD=0 = perfect recovery.)\n");
    cdf_graph_free(true_g);

    cdf_pc_result_free(pc_res);
    cdf_fci_result_free(fci_res);
    cdf_dataset_free(ds);

    printf("\n==========================================\n");
    printf("  Application complete.\n");
    printf("==========================================\n");
    return 0;
}