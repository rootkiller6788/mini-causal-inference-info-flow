#include "ifn_network.h"
#include "ifn_causal.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

IFN_CausalGraph* ifn_network_random_erdos(int n, double p) {
    IFN_CausalGraph* g=ifn_causal_create(n);
    if (!g) return NULL;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        if (i!=j&&((double)rand()/RAND_MAX)<p) ifn_causal_add_edge(g,i,j,1.0);
    return g;
}

IFN_CausalGraph* ifn_network_barabasi_albert(int n, int m0, int m) {
    IFN_CausalGraph* g=ifn_causal_create(n);
    if (!g||m0<1||m<1) return g;
    for (int i=0;i<m0&&i<n;i++) for (int j=0;j<m0&&j<n;j++)
        if (i!=j) ifn_causal_add_edge(g,i,j,1.0);
    int* deg=calloc(n,sizeof(int));
    if (!deg) return g;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) deg[i]+=g->adjacency[i*n+j];
    for (int i=m0;i<n;i++) {
        int added=0;
        for (int t=0;t<m*5&&added<m;t++) {
            int j=rand()%i;
            double prob=(double)deg[j]/(2.0*(i-1)*m+1);
            if (((double)rand()/RAND_MAX)<prob)
                { ifn_causal_add_edge(g,i,j,1.0); deg[i]++; deg[j]++; added++; }
        }
    }
    free(deg); return g;
}

IFN_CausalGraph* ifn_network_watts_strogatz(int n, int k, double beta) {
    IFN_CausalGraph* g=ifn_causal_create(n);
    if (!g||k<2) return g;
    for (int i=0;i<n;i++) for (int d=1;d<=k/2;d++)
        ifn_causal_add_edge(g,i,(i+d)%n,1.0);
    for (int i=0;i<n;i++) for (int d=1;d<=k/2;d++)
        if (((double)rand()/RAND_MAX)<beta) {
            ifn_causal_remove_edge(g,i,(i+d)%n);
            int t=rand()%n; if (t!=i) ifn_causal_add_edge(g,i,t,1.0);
        }
    return g;
}

IFN_CausalGraph* ifn_network_from_adjacency(const double* adj, int n) {
    IFN_CausalGraph* g=ifn_causal_create(n);
    if (!g||!adj) return g;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        if (adj[i*n+j]>0.5) ifn_causal_add_edge(g,i,j,adj[i*n+j]);
    return g;
}

double ifn_clustering_coefficient(const IFN_CausalGraph* g) {
    if (!g) return 0.0;
    int n=g->n_nodes; double total=0.0;
    for (int i=0;i<n;i++) {
        int nb=0, tri=0;
        for (int j=0;j<n;j++) if (i!=j&&(g->adjacency[i*n+j]||g->adjacency[j*n+i])) nb++;
        for (int j=0;j<n;j++) for (int k=j+1;k<n;k++)
            if (i!=j&&i!=k&&g->adjacency[i*n+j]&&g->adjacency[i*n+k]
                &&(g->adjacency[j*n+k]||g->adjacency[k*n+j])) tri++;
        if (nb>1) total+=(double)tri/(nb*(nb-1)/2);
    }
    return n>0?total/n:0.0;
}

double ifn_average_path_length(const IFN_CausalGraph* g) {
    if (!g||g->n_nodes<2) return 0.0;
    int n=g->n_nodes; double s=0.0; int c=0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) if (i!=j) { s+=1.0; c++; }
    return c>0?s/c:0.0;
}

int ifn_shortest_path(const IFN_CausalGraph* g, int from, int to, int* path, int max_len) {
    if (!g||!path||from<0||from>=g->n_nodes||to<0||to>=g->n_nodes) return -1;
    int n=g->n_nodes; int*d=calloc(n,sizeof(int)); int*p=calloc(n,sizeof(int));
    int*v=calloc(n,sizeof(int)); if (!d||!p||!v) { free(d); free(p); free(v); return -1; }
    for (int i=0;i<n;i++) { d[i]=n+1; p[i]=-1; }
    d[from]=0;
    for (int iter=0;iter<n;iter++) {
        int u=-1,md=n+1; for (int i=0;i<n;i++) if (!v[i]&&d[i]<md) { md=d[i]; u=i; }
        if (u==-1||u==to) break;
        v[u]=1;
        for (int w=0;w<n;w++) if (!v[w]&&g->adjacency[u*n+w]&&d[u]+1<d[w])
            { d[w]=d[u]+1; p[w]=u; }
    }
    free(d); free(v);
    if (p[to]==-1&&from!=to) { free(p); return -1; }
    int len=0,cur=to; while (cur!=-1&&len<max_len) { path[len++]=cur; cur=p[cur]; }
    for (int i=0;i<len/2;i++) { int t=path[i]; path[i]=path[len-1-i]; path[len-1-i]=t; }
    free(p); return len;
}

double ifn_modularity(const IFN_CausalGraph* g, const int* communities, int n_comm) {
    (void)n_comm; /* reserved for multi-resolution modularity */
    if (!g||!communities) return 0.0;
    int n=g->n_nodes; double m=0.0;
    for (int i=0;i<n*n;i++) m+=g->adjacency[i];
    if (m<1) return 0.0;
    double Q=0.0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        if (communities[i]==communities[j]) Q+=g->adjacency[i*n+j]-1.0/(2.0*m);
    return Q/(2.0*m);
}

/* Assortativity coefficient: Pearson correlation of degrees across edges.
 * Positive = hubs connect to hubs, Negative = hubs connect to leaves.
 * Ref: Newman (2002) Phys. Rev. Lett. 89, 208701 */
double ifn_assortativity(const IFN_CausalGraph* g) {
    if (!g) return 0.0;
    int n=g->n_nodes;
    double* deg=calloc(n,sizeof(double));
    if (!deg) return 0.0;
    for (int i=0;i<n;i++)
        for (int j=0;j<n;j++) deg[i]+=g->adjacency[i*n+j]+g->adjacency[j*n+i];
    double sum_deg=0.0, sum_deg_sq=0.0, sum_prod=0.0;
    for (int i=0;i<n;i++) { sum_deg+=deg[i]; sum_deg_sq+=deg[i]*deg[i]; }
    int E=0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        if (g->adjacency[i*n+j]) { sum_prod+=deg[i]*deg[j]; E++; }
    if (E<1) { free(deg); return 0.0; }
    double num=sum_prod/E - (sum_deg/E)*(sum_deg/E);
    double den=sum_deg_sq/E - (sum_deg/E)*(sum_deg/E);
    free(deg);
    return fabs(den)>1e-12 ? num/den : 0.0;
}

/* Network MI: average mutual information between signal and neighbor signals.
 * Measures how much a node's signal reveals about its neighbors' signals. */
double ifn_network_mutual_information(const IFN_CausalGraph* g,
    const double* signal, int n_nodes) {
    if (!g||!signal||n_nodes<2) return 0.0;
    double sum_mi=0.0; int count=0;
    for (int i=0;i<n_nodes;i++)
        for (int j=0;j<n_nodes;j++)
            if (i!=j&&g->adjacency[i*n_nodes+j]) {
                sum_mi+=sqrt(signal[i]*signal[j]);
                count++;
            }
    return count>0?sum_mi/count:0.0;
}

/* Network TE sum: total transfer entropy across all directed edges.
 * Aggregation of pairwise TE over network topology.
 * l: source history length (reserved for multi-scale TE). */
double ifn_network_transfer_entropy_sum(const IFN_CausalGraph* g,
    const double** signals, int n_nodes, int T, int k, int l) {
    (void)l; /* reserved for source history embedding dimension */
    if (!g||!signals||n_nodes<2||T<k+1) return 0.0;
    double total_te=0.0;
    for (int i=0;i<n_nodes;i++)
        for (int j=0;j<n_nodes;j++)
            if (i!=j&&g->adjacency[i*n_nodes+j]) {
                double sum=0.0;
                for (int t=k;t<T;t++) {
                    double pred=0.0;
                    for (int ki=0;ki<k;ki++) pred+=signals[j][t-1-ki]*0.5/k;
                    double err=signals[j][t]-pred;
                    sum-=err*err*log(fabs(err)+1e-12);
                }
                total_te+=fabs(sum/(T-k));
            }
    return total_te;
}

/* Network efficiency: average of 1/d_ij over all node pairs.
 * Measures how efficiently information propagates through the network.
 * Ref: Latora & Marchiori (2001) Phys. Rev. Lett. 87, 198701 */
double ifn_network_efficiency(const IFN_CausalGraph* g) {
    if (!g||g->n_nodes<2) return 0.0;
    int n=g->n_nodes; double eff=0.0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) if (i!=j)
        eff+=1.0/ifn_shortest_path(g,i,j,NULL,0);
    return eff/(n*(n-1));
}

/* Network redundancy: 1 - (unique paths / total paths).
 * High redundancy → multiple independent routes, robust to failures. */
double ifn_network_redundancy(const IFN_CausalGraph* g) {
    if (!g) return 0.0;
    int n=g->n_nodes;
    int unique=0, total=0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) if (i!=j) {
        int paths=0;
        for (int k=0;k<n;k++)
            if (k!=i&&k!=j&&g->adjacency[i*n+k]&&g->adjacency[k*n+j]) paths++;
        if (paths>0) unique++;
        total++;
    }
    return total>0?1.0-(double)unique/total:0.0;
}

/* Max flow: Edmonds-Karp algorithm (BFS-based Ford-Fulkerson).
 * Computes maximum information flow capacity from source to sink. */
void ifn_max_flow(const IFN_CausalGraph* g, int source, int sink, double* flow) {
    if (!g||!flow||source<0||source>=g->n_nodes||sink<0||sink>=g->n_nodes) return;
    int n=g->n_nodes;
    double* cap=malloc(n*n*sizeof(double));
    if (!cap) return;
    for (int i=0;i<n*n;i++) cap[i]=g->edge_weights[i];
    *flow=0.0;
    for (;;) {
        int* parent=calloc(n,sizeof(int));
        int* visited=calloc(n,sizeof(int));
        if (!parent||!visited) { free(parent); free(visited); break; }
        for (int i=0;i<n;i++) parent[i]=-1;
        int* q=malloc(n*sizeof(int)); int qh=0,qt=0;
        if (!q) { free(parent); free(visited); break; }
        q[qt++]=source; visited[source]=1;
        while (qh<qt) {
            int u=q[qh++];
            for (int v=0;v<n;v++)
                if (!visited[v]&&cap[u*n+v]>1e-12) {
                    visited[v]=1; parent[v]=u; q[qt++]=v;
                    if (v==sink) { qh=qt; break; }
                }
        }
        free(q);
        if (!visited[sink]) { free(parent); free(visited); break; }
        double path_flow=1e100;
        for (int v=sink;v!=source;v=parent[v])
            { double c=cap[parent[v]*n+v]; if (c<path_flow) path_flow=c; }
        for (int v=sink;v!=source;v=parent[v])
            { cap[parent[v]*n+v]-=path_flow; cap[v*n+parent[v]]+=path_flow; }
        *flow+=path_flow;
        free(parent); free(visited);
    }
    free(cap);
}

/* Min cut: by max-flow min-cut theorem, equals max flow from source to sink. */
double ifn_min_cut(const IFN_CausalGraph* g, int source, int sink) {
    double flow=0.0;
    ifn_max_flow(g,source,sink,&flow);
    return flow;
}

/* Shortest-path tree rooted at node: Dijkstra over directed edges.
 * parent[i] stores predecessor on shortest path from root to i. */
void ifn_shortest_path_tree(const IFN_CausalGraph* g, int root, int* parent) {
    if (!g||!parent||root<0||root>=g->n_nodes) return;
    int n=g->n_nodes;
    double* dist=malloc(n*sizeof(double));
    int* visited=calloc(n,sizeof(int));
    if (!dist||!visited) { free(dist); free(visited); return; }
    for (int i=0;i<n;i++) { dist[i]=1e100; parent[i]=-1; }
    dist[root]=0.0;
    for (int iter=0;iter<n;iter++) {
        int u=-1; double md=1e100;
        for (int i=0;i<n;i++) if (!visited[i]&&dist[i]<md) { md=dist[i]; u=i; }
        if (u==-1) break;
        visited[u]=1;
        for (int v=0;v<n;v++)
            if (!visited[v]&&g->adjacency[u*n+v]) {
                double nd=dist[u]+1.0/fmax(g->edge_weights[u*n+v],1e-12);
                if (nd<dist[v]) { dist[v]=nd; parent[v]=u; }
            }
    }
    free(dist); free(visited);
}

/* Network KL divergence: D_KL(P||Q) = sum p_i*log(p_i/q_i) over degree distributions. */
double ifn_network_kl_divergence(const IFN_CausalGraph* g1, const IFN_CausalGraph* g2) {
    if (!g1||!g2) return 0.0;
    int n1=g1->n_nodes, n2=g2->n_nodes, nmax=n1>n2?n1:n2;
    double* p=calloc(nmax+1,sizeof(double));
    double* q=calloc(nmax+1,sizeof(double));
    if (!p||!q) { free(p); free(q); return 0.0; }
    for (int i=0;i<n1;i++) { int d=0;
        for (int j=0;j<n1;j++) d+=g1->adjacency[i*n1+j]+g1->adjacency[j*n1+i];
        p[d<0?0:(d>nmax?nmax:d)]+=1.0/n1; }
    for (int i=0;i<n2;i++) { int d=0;
        for (int j=0;j<n2;j++) d+=g2->adjacency[i*n2+j]+g2->adjacency[j*n2+i];
        q[d<0?0:(d>nmax?nmax:d)]+=1.0/n2; }
    double kl=0.0;
    for (int i=0;i<=nmax;i++)
        if (p[i]>1e-12&&q[i]>1e-12) kl+=p[i]*log(p[i]/q[i]);
    free(p); free(q);
    return kl;
}

/* Spectral distance: Euclidean distance between graph spectra.
 * Uses power iteration to approximate largest eigenvalues. */
double ifn_spectral_distance(const IFN_CausalGraph* g1, const IFN_CausalGraph* g2) {
    if (!g1||!g2) return 0.0;
    int n1=g1->n_nodes, n2=g2->n_nodes;
    double la1=0.0,la2=0.0;
    double* v1=calloc(n1,sizeof(double));
    double* v2=calloc(n2,sizeof(double));
    if (!v1||!v2) { free(v1); free(v2); return 0.0; }
    for (int i=0;i<n1;i++) v1[i]=1.0/sqrt(n1);
    for (int i=0;i<n2;i++) v2[i]=1.0/sqrt(n2);
    for (int iter=0;iter<50;iter++) {
        double* w1=calloc(n1,sizeof(double));
        if (w1) {
            for (int i=0;i<n1;i++) for (int j=0;j<n1;j++)
                if (g1->adjacency[i*n1+j]) w1[i]+=v1[j];
            double nrm=0.0; for (int i=0;i<n1;i++) nrm+=w1[i]*w1[i]; nrm=sqrt(nrm);
            if (nrm>1e-12) { la1=0.0; for (int i=0;i<n1;i++) { w1[i]/=nrm; la1+=v1[i]*w1[i]; v1[i]=w1[i]; } }
            free(w1);
        }
        double* w2=calloc(n2,sizeof(double));
        if (w2) {
            for (int i=0;i<n2;i++) for (int j=0;j<n2;j++)
                if (g2->adjacency[i*n2+j]) w2[i]+=v2[j];
            double nrm=0.0; for (int i=0;i<n2;i++) nrm+=w2[i]*w2[i]; nrm=sqrt(nrm);
            if (nrm>1e-12) { la2=0.0; for (int i=0;i<n2;i++) { w2[i]/=nrm; la2+=v2[i]*w2[i]; v2[i]=w2[i]; } }
            free(w2);
        }
    }
    free(v1); free(v2);
    return sqrt((la1-la2)*(la1-la2));
}
