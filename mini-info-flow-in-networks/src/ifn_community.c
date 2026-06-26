/* ifn_community.c — Community detection for information flow networks.
 *
 * Knowledge points implemented:
 *   1. Edge betweenness centrality — Girvan-Newman foundation
 *   2. Girvan-Newman algorithm — hierarchical community detection
 *   3. Louvain method — greedy modularity optimization
 *   4. Infomap — information-theoretic community detection
 *   5. Community mutual information — structure vs flow alignment
 *   6. NMI — normalized mutual information for partition comparison
 *
 * Ref: Girvan & Newman (2002) PNAS; Blondel et al. (2008) J. Stat. Mech.
 *      Rosvall & Bergstrom (2008) PNAS; Newman (2006) Phys. Rev. E 74, 036104
 */
#include "ifn_community.h"
#include "ifn_causal.h"
#include "ifn_network.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- Edge betweenness (Brandes algorithm) ---- */
void ifn_edge_betweenness(const IFN_CausalGraph* g, double* edge_between) {
    if (!g||!edge_between) return;
    int n=g->n_nodes;
    memset(edge_between,0,n*n*sizeof(double));

    for (int s=0;s<n;s++) {
        /* BFS from source s */
        int* dist=calloc(n,sizeof(int));
        double* sigma=calloc(n,sizeof(double));
        double* delta=calloc(n,sizeof(double));
        int* q=malloc(n*sizeof(int));
        int** pred=malloc(n*sizeof(int*));
        int* pred_count=calloc(n,sizeof(int));
        if (!dist||!sigma||!delta||!q||!pred||!pred_count)
            { free(dist);free(sigma);free(delta);free(q);free(pred);free(pred_count);continue; }

        for (int v=0;v<n;v++) dist[v]=-1;
        dist[s]=0; sigma[s]=1.0;
        int qh=0,qt=0; q[qt++]=s;

        while (qh<qt) {
            int v=q[qh++];
            for (int w=0;w<n;w++) {
                if (!g->adjacency[v*n+w]) continue;
                if (dist[w]==-1) {
                    dist[w]=dist[v]+1;
                    q[qt++]=w;
                }
                if (dist[w]==dist[v]+1) {
                    sigma[w]+=sigma[v];
                    /* Track predecessor */
                    pred[w]=realloc(pred[w],(pred_count[w]+1)*sizeof(int));
                    if (pred[w]) pred[w][pred_count[w]++]=v;
                }
            }
        }

        /* Back-propagate dependency */
        for (int i=qt-1;i>=0;i--) {
            int w=q[i];
            for (int p=0;p<pred_count[w];p++) {
                int v=pred[w][p];
                double c=sigma[v]/sigma[w]*(1.0+delta[w]);
                int idx=v*n+w;
                if (idx>=0&&idx<n*n) edge_between[idx]+=c;
                delta[v]+=c;
            }
        }

        free(dist);free(sigma);free(delta);free(q);
        for (int v=0;v<n;v++) free(pred[v]);
        free(pred); free(pred_count);
    }
}

/* ---- Girvan-Newman algorithm ---- */
void ifn_girvan_newman(const IFN_CausalGraph* g,
    int* community_labels, int* n_communities) {
    if (!g||!community_labels||!n_communities||g->n_nodes<2) return;
    int n=g->n_nodes;

    /* Working copy of graph */
    IFN_CausalGraph* wg=ifn_causal_create(n);
    if (!wg) return;
    for (int i=0;i<n*n;i++) { wg->adjacency[i]=g->adjacency[i]; wg->edge_weights[i]=g->edge_weights[i]; }

    double* ebet=malloc(n*n*sizeof(double));
    if (!ebet) { ifn_causal_free(wg); return; }

    /* Remove edges until disconnection */
    int n_components=n;
    for (int iter=0;iter<n*n&&n_components<10;iter++) {
        ifn_edge_betweenness(wg,ebet);

        /* Find max betweenness edge and remove */
        double max_eb=0.0; int max_i=-1;
        for (int i=0;i<n*n;i++)
            if (wg->adjacency[i]&&ebet[i]>max_eb) { max_eb=ebet[i]; max_i=i; }
        if (max_i<0) break;

        int from=max_i/n, to=max_i%n;
        ifn_causal_remove_edge(wg,from,to);

        /* Count components via DFS */
        int* visited=calloc(n,sizeof(int));
        if (!visited) break;
        n_components=0;
        for (int v=0;v<n;v++) if (!visited[v]) {
            n_components++;
            int* stack=malloc(n*sizeof(int)); int sp=0;
            stack[sp++]=v; visited[v]=n_components;
            while (sp>0) {
                int u=stack[--sp];
                for (int w=0;w<n;w++)
                    if (!visited[w]&&(wg->adjacency[u*n+w]||wg->adjacency[w*n+u]))
                        { visited[w]=n_components; stack[sp++]=w; }
            }
            free(stack);
        }
        for (int v=0;v<n;v++) community_labels[v]=visited[v]-1;
        free(visited);
    }

    *n_communities=n_components;
    free(ebet); ifn_causal_free(wg);
}

/* ---- Louvain method ---- */
double ifn_louvain(const IFN_CausalGraph* g,
    int* community_labels, int max_iter) {
    if (!g||!community_labels||g->n_nodes<1) return 0.0;
    int n=g->n_nodes;

    /* Initialize: each node in own community */
    for (int i=0;i<n;i++) community_labels[i]=i;

    /* Compute total edge weight */
    double m=0.0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) m+=g->adjacency[i*n+j]?1.0:0.0;
    if (m<1) return 0.0;

    /* Node degrees */
    double* k=calloc(n,sizeof(double));
    if (!k) return 0.0;
    for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
            k[i]+=g->adjacency[i*n+j]?1.0:0.0;

    for (int iter=0;iter<max_iter;iter++) {
        int moved=0;
        for (int i=0;i<n;i++) {
            /* Neighbor community weights */
            int n_comms=0;
            int* comms=calloc(n,sizeof(int));
            double* k_in=calloc(n,sizeof(double));
            double* sigma_tot=calloc(n,sizeof(double));
            if (!comms||!k_in||!sigma_tot) { free(comms);free(k_in);free(sigma_tot);continue; }

            for (int j=0;j<n;j++)
                if (g->adjacency[i*n+j]) {
                    int c=community_labels[j];
                    int found=-1;
                    for (int cc=0;cc<n_comms;cc++) if (comms[cc]==c) { found=cc; break; }
                    if (found<0) { comms[n_comms]=c; k_in[n_comms]=1.0; sigma_tot[n_comms]=k[j]; n_comms++; }
                    else k_in[found]+=1.0;
                }

            /* Compute modularity gain for moving i to each community */
            double best_gain=0.0; int best_comm=community_labels[i];
            for (int c=0;c<n_comms;c++) {
                double gain=k_in[c]/m - k[i]*(sigma_tot[c])/(2.0*m*m);
                if (gain>best_gain) { best_gain=gain; best_comm=comms[c]; }
            }
            if (best_gain>1e-12||community_labels[i]!=best_comm) {
                community_labels[i]=best_comm; moved++;
            }
            free(comms); free(k_in); free(sigma_tot);
        }
        if (!moved) break;
    }

    /* Compute modularity Q */
    double Q=0.0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        if (community_labels[i]==community_labels[j])
            Q+=(g->adjacency[i*n+j]?1.0:0.0)-k[i]*k[j]/(2.0*m);
    Q/=2.0*m;
    free(k);
    return Q;
}

/* ---- Infomap (simplified map equation) ---- */
double ifn_infomap(const IFN_CausalGraph* g,
    int* community_labels, int max_iter) {
    if (!g||!community_labels||g->n_nodes<1) return 0.0;
    int n=g->n_nodes;

    /* Initialize with random partition */
    for (int i=0;i<n;i++) community_labels[i]=i%((n>10?n/5:2));

    /* Compute degrees */
    double* deg=calloc(n,sizeof(double));
    if (!deg) return 0.0;
    for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
            deg[i]+=g->adjacency[i*n+j]?1.0:0.0;
    double total_deg=0.0; for (int i=0;i<n;i++) total_deg+=deg[i];

    double L=0.0;
    for (int iter=0;iter<max_iter;iter++) {
        /* Compute map equation components:
         * q = exit probability per module
         * p_i = within-module visit rates */
        double L_new=0.0;
        for (int i=0;i<n;i++) {
            double q=deg[i]/total_deg;
            /* Entropy of within-module code */
            double p_sum=0.0;
            for (int j=0;j<n;j++)
                if (community_labels[j]==community_labels[i])
                    p_sum+=deg[j]/total_deg;
            if (q>1e-12&&p_sum>1e-12)
                L_new+=q*log(p_sum/q);
        }
        if (fabs(L_new-L)<1e-6) { L=L_new; break; }
        L=L_new;

        /* Move nodes to optimize description length */
        for (int i=0;i<n;i++) {
            double best_L=1e100; int best_c=community_labels[i];
            for (int c=0;c<n;c++) if (c!=community_labels[i]) {
                int old_c=community_labels[i];
                community_labels[i]=c;
                double trial_L=0.0;
                for (int v=0;v<n;v++) {
                    double qv=deg[v]/total_deg, ps=0.0;
                    for (int w=0;w<n;w++)
                        if (community_labels[w]==community_labels[v]) ps+=deg[w]/total_deg;
                    if (qv>1e-12&&ps>1e-12) trial_L+=qv*log(ps/qv);
                }
                if (trial_L<best_L) { best_L=trial_L; best_c=c; }
                community_labels[i]=old_c;
            }
            community_labels[i]=best_c;
        }
    }

    free(deg);
    return L>0?L:0.0;
}

/* ---- Community MI ---- */
double ifn_community_mi(const IFN_CausalGraph* g,
    const int* community_labels, int n_communities,
    const double* flow_pattern, int n_edges) {
    if (!g||!community_labels||!flow_pattern||n_communities<1) return 0.0;
    int n=g->n_nodes;
    double* p_comm=calloc(n_communities,sizeof(double));
    double* p_flow=calloc(n_communities,sizeof(double));
    double* joint=calloc(n_communities*n_communities,sizeof(double));
    if (!p_comm||!p_flow||!joint) { free(p_comm);free(p_flow);free(joint);return 0.0; }

    for (int i=0;i<n;i++) p_comm[community_labels[i]]+=1.0/n;
    double total_flow=0.0; for (int i=0;i<n_edges;i++) total_flow+=flow_pattern[i];
    for (int e=0;e<n_edges;e++) {
        int from=g->edges?g->edges[e*2]:0, to=g->edges?g->edges[e*2+1]:0;
        if (from<n&&to<n&&total_flow>1e-12) {
            int c1=community_labels[from], c2=community_labels[to];
            p_flow[c1]+=flow_pattern[e]/total_flow;
            joint[c1*n_communities+c2]+=flow_pattern[e]/total_flow;
        }
    }

    double mi=0.0;
    for (int a=0;a<n_communities;a++)
        for (int b=0;b<n_communities;b++)
            if (joint[a*n_communities+b]>1e-12&&p_flow[a]>1e-12&&p_comm[b]>1e-12)
                mi+=joint[a*n_communities+b]*
                    log(joint[a*n_communities+b]/(p_flow[a]*p_comm[b]));
    free(p_comm);free(p_flow);free(joint);
    return mi>0?mi:0.0;
}

/* ---- Normalized Mutual Information (NMI) between partitions ---- */
double ifn_nmi(const int* labels_a, const int* labels_b, int n) {
    if (!labels_a||!labels_b||n<1) return 0.0;

    int max_a=0,max_b=0;
    for (int i=0;i<n;i++) {
        if (labels_a[i]>max_a) max_a=labels_a[i];
        if (labels_b[i]>max_b) max_b=labels_b[i];
    }
    int na=max_a+1, nb=max_b+1;

    double* p_a=calloc(na,sizeof(double));
    double* p_b=calloc(nb,sizeof(double));
    double* joint=calloc(na*nb,sizeof(double));
    if (!p_a||!p_b||!joint) { free(p_a);free(p_b);free(joint);return 0.0; }

    for (int i=0;i<n;i++) {
        int a=labels_a[i]%na, b=labels_b[i]%nb;
        p_a[a]+=1.0/n; p_b[b]+=1.0/n;
        joint[a*nb+b]+=1.0/n;
    }

    double Ha=0.0,Hb=0.0,MI=0.0;
    for (int a=0;a<na;a++) if (p_a[a]>1e-12) Ha-=p_a[a]*log(p_a[a]);
    for (int b=0;b<nb;b++) if (p_b[b]>1e-12) Hb-=p_b[b]*log(p_b[b]);
    for (int a=0;a<na;a++) for (int b=0;b<nb;b++)
        if (joint[a*nb+b]>1e-12&&p_a[a]>1e-12&&p_b[b]>1e-12)
            MI+=joint[a*nb+b]*log(joint[a*nb+b]/(p_a[a]*p_b[b]));

    double denom=Ha>1e-12&&Hb>1e-12?sqrt(Ha*Hb):0.0;
    free(p_a); free(p_b); free(joint);
    return denom>1e-12?MI/denom:0.0;
}
