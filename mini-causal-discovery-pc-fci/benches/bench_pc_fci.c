/* bench_pc_fci.c - Performance Benchmarks for PC/FCI Causal Discovery
 *
 * Benchmarks key operations to quantify computational scalability:
 *   - CI test throughput (tests/sec)
 *   - Skeleton search scaling (vs number of variables p)
 *   - d-separation check throughput
 *   - PC vs FCI runtime comparison
 *   - Bootstrap overhead
 *
 * Reports operations per second and asymptotic scaling behavior.
 *
 * Reference: Kalisch & Buhlmann (2007), "Estimating High-Dimensional
 *   Directed Acyclic Graphs with the PC-Algorithm". JMLR.
 */

#include "cdf_core.h"
#include "cdf_citest.h"
#include "cdf_graph.h"
#include "cdf_orientation.h"
#include "cdf_pc.h"
#include "cdf_fci.h"
#include "cdf_bootstrap.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Generate random correlated Gaussian data with p variables */
static double* gen_random_data(int N, int p)
{
    double *d = (double*)calloc((size_t)p * N, sizeof(double));
    if (!d) return NULL;
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < N; j++) {
            d[i * N + j] = ((double)((j * (i + 1) * 7 + 11) % 1000) / 1000.0 - 0.5);
        }
    }
    return d;
}

/* Generate chain data X1->X2->...->Xp */
static double* gen_chain_data(int N, int p)
{
    double *d = (double*)calloc((size_t)p * N, sizeof(double));
    if (!d) return NULL;
    for (int j = 0; j < N; j++) {
        double prev = 0.0;
        for (int i = 0; i < p; i++) {
            double e = ((double)((j * (i + 2) * 7 + i * 13) % 1000) / 1000.0 - 0.5);
            if (i == 0) prev = e;
            else prev = 0.7 * prev + e;
            d[i * N + j] = prev;
        }
    }
    return d;
}

static double elapsed_sec(clock_t start) {
    return (double)(clock() - start) / (double)CLOCKS_PER_SEC;
}

int main(void)
{
    printf("==========================================\n");
    printf("  PC/FCI Causal Discovery Benchmarks      \n");
    printf("==========================================\n\n");

    /* ── Benchmark 1: CI Test Throughput ── */
    printf("[Bench 1] CI Test Throughput\n");
    {
        int N = 500, p = 20;
        double *data = gen_random_data(N, p);
        CdfDataset *ds = cdf_dataset_create(data, N, p);
        free(data);
        if (ds) {
            int n_tests = 1000;
            clock_t t0 = clock();
            for (int k = 0; k < n_tests; k++) {
                cdf_citest_fisher_z(ds, k % p, (k + 1) % p, NULL, 0, 0.05);
            }
            double dt = elapsed_sec(t0);
            printf("  Fisher-z tests: %d in %.3f sec (%.0f tests/sec)\n",
                   n_tests, dt, dt > 0 ? n_tests / dt : 0.0);
            cdf_dataset_free(ds);
        }
    }

    /* ── Benchmark 2: Correlation Matrix ── */
    printf("\n[Bench 2] Correlation Matrix Computation\n");
    {
        int N = 500, p = 15;
        double *data = gen_random_data(N, p);
        CdfDataset *ds = cdf_dataset_create(data, N, p);
        free(data);
        if (ds) {
            double *corr = (double*)malloc((size_t)p * p * sizeof(double));
            clock_t t0 = clock();
            cdf_citest_corr_matrix(ds, corr);
            double dt = elapsed_sec(t0);
            printf("  %dx%d correlation matrix: %.4f sec\n", p, p, dt);
            free(corr);
            cdf_dataset_free(ds);
        }
    }

    /* ── Benchmark 3: d-Separation Throughput ── */
    printf("\n[Bench 3] d-Separation Throughput\n");
    {
        int p = 10;
        CdfGraph *g = cdf_graph_create(p);
        /* Build a random-ish DAG */
        for (int i = 0; i < p - 1; i++)
            cdf_graph_add_edge(g, i, i + 1, CDF_EDGE_DIRECTED);
        /* Add some extra edges */
        cdf_graph_add_edge(g, 0, 3, CDF_EDGE_DIRECTED);
        cdf_graph_add_edge(g, 2, 5, CDF_EDGE_DIRECTED);
        cdf_graph_add_edge(g, 4, 7, CDF_EDGE_DIRECTED);

        int Z[] = {3, 5};
        int n_dsep = 1000;
        clock_t t0 = clock();
        for (int k = 0; k < n_dsep; k++) {
            cdf_graph_d_separated(g, k % p, (k + 5) % p, Z, 2);
        }
        double dt = elapsed_sec(t0);
        printf("  d-sep checks: %d in %.3f sec (%.0f checks/sec)\n",
               n_dsep, dt, dt > 0 ? n_dsep / dt : 0.0);
        cdf_graph_free(g);
    }

    /* ── Benchmark 4: Skeleton Search Scaling ── */
    printf("\n[Bench 4] Skeleton Search Scaling vs p\n");
    {
        int N = 200;
        int p_sizes[] = {3, 5, 7, 10, 0};
        printf("  %-4s %-10s %-10s %-12s\n", "p", "Edges init", "CI tests", "Time(sec)");
        for (int k = 0; p_sizes[k] > 0; k++) {
            int p = p_sizes[k];
            double *data = gen_chain_data(N, p);
            CdfDataset *ds = cdf_dataset_create(data, N, p);
            free(data);
            if (!ds) continue;

            CdfGraph *skel = cdf_graph_create(p);
            CdfPCConfig cfg = cdf_pc_config_default();
            cfg.alpha = 0.05;

            clock_t t0 = clock();
            int n_ci = cdf_graph_skeleton_pc(skel, ds, &cfg);
            double dt = elapsed_sec(t0);

            printf("  %-4d %-10d %-10d %-12.4f\n",
                   p, p * (p - 1) / 2, n_ci, dt);

            cdf_graph_free(skel);
            cdf_dataset_free(ds);
        }
    }

    /* ── Benchmark 5: PC vs FCI Runtime ── */
    printf("\n[Bench 5] PC vs FCI Runtime Comparison\n");
    {
        int N = 300, p = 6;
        double *data = gen_chain_data(N, p);
        CdfDataset *ds = cdf_dataset_create(data, N, p);
        free(data);
        if (ds) {
            CdfPCConfig pc_cfg = cdf_pc_config_default();
            pc_cfg.alpha = 0.05;

            clock_t t0 = clock();
            CdfPCResult *pc = cdf_pc_run(ds, &pc_cfg);
            double dt_pc = elapsed_sec(t0);
            printf("  PC:  %.4f sec, CI tests=%d\n",
                   dt_pc, pc ? pc->n_ci_tests : -1);
            if (pc) cdf_pc_result_free(pc);

            CdfFCIConfig fci_cfg = cdf_fci_config_default();
            fci_cfg.alpha = 0.05;

            clock_t t1 = clock();
            CdfFCIResult *fci = cdf_fci_run(ds, &fci_cfg);
            double dt_fci = elapsed_sec(t1);
            printf("  FCI: %.4f sec, CI tests=%d\n",
                   dt_fci, fci ? fci->n_ci_tests : -1);
            if (fci) cdf_fci_result_free(fci);

            printf("  Ratio (FCI/PC): %.2fx\n",
                   dt_pc > 0 ? dt_fci / dt_pc : 0.0);
            cdf_dataset_free(ds);
        }
    }

    /* ── Benchmark 6: Orientation Rules ── */
    printf("\n[Bench 6] Orientation Rules Throughput\n");
    {
        int p = 8;
        CdfGraph *g = cdf_graph_create(p);
        /* Dense skeleton */
        for (int i = 0; i < p; i++)
            for (int j = i + 1; j < p; j++)
                cdf_graph_add_edge(g, i, j, CDF_EDGE_UNDIRECTED);
        /* Add some v-structures */
        cdf_graph_add_edge(g, 0, 1, CDF_EDGE_DIRECTED);
        cdf_graph_add_edge(g, 2, 1, CDF_EDGE_DIRECTED);

        clock_t t0 = clock();
        CdfOrientResult r = cdf_orient_pc_rules(g);
        double dt = elapsed_sec(t0);
        printf("  PC rules on %d nodes: %.4f sec, oriented %d edges\n",
               p, dt, r.total_oriented);
        cdf_graph_free(g);
    }

    /* ── Benchmark 7: Bootstrap Overhead ── */
    printf("\n[Bench 7] Bootstrap Edge Stability Overhead\n");
    {
        int N = 200, p = 4;
        double *data = gen_chain_data(N, p);
        CdfDataset *ds = cdf_dataset_create(data, N, p);
        free(data);
        if (ds) {
            CdfBootstrapConfig bs_cfg = cdf_bootstrap_config_default();
            bs_cfg.n_bootstrap = 20;
            bs_cfg.use_pc = true;
            bs_cfg.verbose = false;

            clock_t t0 = clock();
            CdfBootstrapResult *bs = cdf_bootstrap_edges(ds, &bs_cfg);
            double dt = elapsed_sec(t0);
            printf("  Bootstrap B=%d on p=%d, N=%d: %.3f sec\n",
                   bs_cfg.n_bootstrap, p, N, dt);
            printf("  Per-replicate: %.4f sec\n",
                   bs_cfg.n_bootstrap > 0 ? dt / bs_cfg.n_bootstrap : 0.0);
            if (bs) {
                printf("  Stable edges: ");
                int n_stable = 0;
                for (int k = 0; k < bs->n_edges; k++)
                    if (bs->edges[k].is_stable) n_stable++;
                printf("%d out of %d pairs\n", n_stable, bs->n_edges);
                cdf_bootstrap_result_free(bs);
            }
            cdf_dataset_free(ds);
        }
    }

    printf("\n==========================================\n");
    printf("  Benchmarks complete.\n");
    printf("==========================================\n");
    return 0;
}
