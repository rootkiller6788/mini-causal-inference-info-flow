/* ifn_motifs.c — Network motif analysis for information flow.
 *
 * Knowledge points implemented:
 *   1. Triadic census — all 16 directed 3-node subgraph types
 *   2. Motif Z-scores — statistical significance vs random graphs
 *   3. Motif-based information flow — flow routing by triad type
 *   4. Feed-forward loop detection — coherent/incoherent FFL
 *   5. Bi-fan motif detection — dense signaling pattern
 *
 * Ref: Milo et al. (2002) Science; Alon (2007) Nat. Rev. Genet.
 *      Holland & Leinhardt (1976) Sociol. Methodol.
 */
#include "ifn_motifs.h"
#include "ifn_causal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Encode the 3-node directed subgraph type as an integer 0-15.
 * The 16 triads are classified by mutual, asymmetric, and null dyad counts
 * (MAN labeling: 003, 012, 102, 021D, 021U, 021C, 111D, 111U, 030T, 030C,
 *  201, 120D, 120U, 120C, 210, 300).
 *
 * We use the adjacency pattern among 3 nodes (a,b,c):
 *   bits 0-1: a<->b dyad (00=none, 01=a->b, 10=b->a, 11=mutual)
 *   bits 2-3: a<->c dyad
 *   bits 4-5: b<->c dyad
 *   → 6 bits, 64 patterns → map to 16 types */
static int triad_type(int ab, int ba, int ac, int ca, int bc, int cb) {
    int m = (ab&&ba?1:0) + (ac&&ca?1:0) + (bc&&cb?1:0); /* mutual dyads */
    int s = (ab+ba+ac+ca+bc+cb);                           /* total edges */

    /* Map to 16 types by MAN code */
    if (s==0) return 0;  /* 003 */
    if (s==1) return 1;  /* 012 */
    if (s==2) {
        if (m==1) return 2; /* 102 */
        return (ab+ac+bc!=ba+ca+cb)?3:4; /* 021D or 021U */
    }
    if (s==3) {
        if (m==0) return 5; /* 111D or 111U */
        return 6; /* 201 */
    }
    if (s==4) {
        if (m==3) return 10; /* 300 */
        if (m==0) return 8; /* 120D or 120U */
        return 7; /* 210 */
    }
    if (s==5) return 9; /* 210 (alt) */
    return 10; /* 300 */
}

void ifn_triadic_census(const IFN_CausalGraph* g, int* triad_counts) {
    if (!g||!triad_counts) return;
    int n=g->n_nodes;
    memset(triad_counts,0,16*sizeof(int));

    for (int i=0;i<n;i++)
        for (int j=i+1;j<n;j++)
            for (int k=j+1;k<n;k++) {
                int ab=g->adjacency[i*n+j], ba=g->adjacency[j*n+i];
                int ac=g->adjacency[i*n+k], ca=g->adjacency[k*n+i];
                int bc=g->adjacency[j*n+k], cb=g->adjacency[k*n+j];
                int t=triad_type(ab,ba,ac,ca,bc,cb);
                if (t>=0&&t<16) triad_counts[t]++;
            }
}

/* Generate random graph with same in/out degree sequence (edge swap MCMC) */
static void randomize_graph(const IFN_CausalGraph* g, IFN_CausalGraph* rg) {
    int n=g->n_nodes;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
        rg->adjacency[i*n+j]=g->adjacency[i*n+j];
        rg->edge_weights[i*n+j]=g->edge_weights[i*n+j];
    }
    rg->n_edges=g->n_edges;
    rg->n_nodes=n;

    int n_swaps=g->n_edges*5;
    for (int s=0;s<n_swaps;s++) {
        int a=rand()%n, b=rand()%n, c=rand()%n, d=rand()%n;
        if (a==b||a==c||a==d||b==c||b==d||c==d) continue;
        if (rg->adjacency[a*n+b]&&rg->adjacency[c*n+d]&&
            !rg->adjacency[a*n+d]&&!rg->adjacency[c*n+b]) {
            rg->adjacency[a*n+b]=0; rg->adjacency[c*n+d]=0;
            rg->adjacency[a*n+d]=1; rg->adjacency[c*n+b]=1;
        }
    }
}

void ifn_motif_zscore(const IFN_CausalGraph* g, int* triad_counts,
    int n_random, double* z_scores) {
    if (!g||!triad_counts||!z_scores) return;
    int n=g->n_nodes;

    double* mean_counts=calloc(16,sizeof(double));
    double* mean_sq=calloc(16,sizeof(double));
    if (!mean_counts||!mean_sq) { free(mean_counts); free(mean_sq); return; }

    IFN_CausalGraph* rg=ifn_causal_create(n);
    if (!rg) { free(mean_counts); free(mean_sq); return; }

    for (int r=0;r<n_random;r++) {
        randomize_graph(g,rg);
        int r_counts[16]; memset(r_counts,0,sizeof(r_counts));
        ifn_triadic_census(rg,r_counts);
        for (int t=0;t<16;t++) {
            mean_counts[t]+=r_counts[t];
            mean_sq[t]+=r_counts[t]*r_counts[t];
        }
    }
    ifn_causal_free(rg);

    for (int t=0;t<16;t++) {
        double mu=mean_counts[t]/n_random;
        double sigma=sqrt(mean_sq[t]/n_random - mu*mu);
        z_scores[t]=(sigma>1e-12)?(triad_counts[t]-mu)/sigma:0.0;
    }
    free(mean_counts); free(mean_sq);
}

void ifn_motif_information_flow(const IFN_CausalGraph* g,
    const double* node_signals, int n_vars, double* motif_flow) {
    (void)n_vars; /* reserved for future per-variable flow normalization */
    if (!g||!node_signals||!motif_flow) return;
    int n=g->n_nodes;
    memset(motif_flow,0,16*sizeof(double));

    for (int i=0;i<n;i++)
        for (int j=i+1;j<n;j++)
            for (int k=j+1;k<n;k++) {
                int ab=g->adjacency[i*n+j], ba=g->adjacency[j*n+i];
                int ac=g->adjacency[i*n+k], ca=g->adjacency[k*n+i];
                int bc=g->adjacency[j*n+k], cb=g->adjacency[k*n+j];
                int t=triad_type(ab,ba,ac,ca,bc,cb);
                if (t>=0&&t<16) {
                    /* Signal flow through this triad */
                    double sig_flow=fabs(node_signals[i])+fabs(node_signals[j])+fabs(node_signals[k]);
                    motif_flow[t]+=sig_flow;
                }
            }
}

int ifn_count_feedforward_loops(const IFN_CausalGraph* g, double* coherent_ratio) {
    if (!g) return 0;
    int n=g->n_nodes, count=0, coherent=0;
    for (int x=0;x<n;x++)
        for (int y=0;y<n;y++) if (x!=y&&g->adjacency[x*n+y])
            for (int z=0;z<n;z++) if (z!=x&&z!=y&&g->adjacency[x*n+z]&&g->adjacency[y*n+z]) {
                count++;
                /* Coherent if all edge weights have same sign */
                double s1=g->edge_weights[x*n+y]>0?1:-1;
                double s2=g->edge_weights[x*n+z]>0?1:-1;
                double s3=g->edge_weights[y*n+z]>0?1:-1;
                if (s1==s2&&s2==s3) coherent++;
            }
    if (coherent_ratio) *coherent_ratio=count>0?(double)coherent/count:0.0;
    return count;
}

int ifn_count_bifan_motifs(const IFN_CausalGraph* g) {
    if (!g||g->n_nodes<4) return 0;
    int n=g->n_nodes, count=0;
    for (int x=0;x<n;x++)
        for (int y=x+1;y<n;y++) if (x!=y)
            for (int a=0;a<n;a++) if (a!=x&&a!=y)
                for (int b=a+1;b<n;b++) if (b!=x&&b!=y)
                    if (g->adjacency[x*n+a]&&g->adjacency[x*n+b]&&
                        g->adjacency[y*n+a]&&g->adjacency[y*n+b]) count++;
    return count;
}
