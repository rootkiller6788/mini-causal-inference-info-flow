#include "ifn_dynamic.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

IFN_NetworkDynamics ifn_network_dynamics_create(const IFN_CausalGraph* g,
    const double* initial_info, int n_steps) {
    IFN_NetworkDynamics nd; memset(&nd,0,sizeof(nd));
    if (!g||!initial_info||n_steps<1) return nd;
    int n=g->n_nodes;
    nd.n_nodes=n; nd.n_steps=n_steps;
    nd.node_entropy=calloc(n,sizeof(double));
    nd.node_inflow=calloc(n,sizeof(double));
    nd.node_outflow=calloc(n,sizeof(double));
    nd.node_centrality=calloc(n,sizeof(double));
    if (!nd.node_entropy||!nd.node_inflow||!nd.node_outflow||!nd.node_centrality) {
        ifn_dynamics_free(&nd); return nd;
    }
    for (int i=0;i<n;i++) nd.node_entropy[i]=initial_info[i];
    return nd;
}

void ifn_dynamics_free(IFN_NetworkDynamics* nd) {
    if (nd) { free(nd->node_entropy); free(nd->node_inflow); free(nd->node_outflow); free(nd->node_centrality); }
}

void ifn_dynamics_step(IFN_NetworkDynamics* nd, const IFN_CausalGraph* g) {
    if (!nd||!g||nd->n_nodes!=g->n_nodes) return;
    int n=g->n_nodes;
    double* new_entropy=malloc(n*sizeof(double));
    if (!new_entropy) return;
    nd->total_flow=0.0;
    for (int i=0;i<n;i++) {
        double inflow=0.0, outflow=0.0;
        for (int j=0;j<n;j++) {
            if (g->adjacency[j*g->n_nodes+i])
                inflow+=nd->node_entropy[j]*g->edge_weights[j*g->n_nodes+i];
            if (g->adjacency[i*g->n_nodes+j])
                outflow+=nd->node_entropy[i]*g->edge_weights[i*g->n_nodes+j];
        }
        nd->node_inflow[i]=inflow; nd->node_outflow[i]=outflow;
        new_entropy[i]=nd->node_entropy[i]+0.1*(inflow-outflow);
        if (new_entropy[i]<0) new_entropy[i]=0;
        nd->total_flow+=inflow+outflow;
    }
    memcpy(nd->node_entropy,new_entropy,n*sizeof(double));
    free(new_entropy);
}

void ifn_dynamics_run(IFN_NetworkDynamics* nd, const IFN_CausalGraph* g, int steps) {
    for (int s=0;s<steps;s++) ifn_dynamics_step(nd,g);
}

double ifn_pagerank(const IFN_CausalGraph* g, double damping, int max_iter,
    double tol, double* pr) {
    if (!g||!pr) return 0.0;
    int n=g->n_nodes;
    double* new_pr=malloc(n*sizeof(double));
    if (!new_pr) return 0.0;
    for (int i=0;i<n;i++) pr[i]=1.0/n;
    for (int iter=0;iter<max_iter;iter++) {
        for (int i=0;i<n;i++) {
            new_pr[i]=(1.0-damping)/n;
            for (int j=0;j<n;j++)
                if (g->adjacency[j*g->n_nodes+i]) {
                    int out_deg=0;
                    for (int k=0;k<n;k++) out_deg+=g->adjacency[j*g->n_nodes+k];
                    if (out_deg>0) new_pr[i]+=damping*pr[j]/out_deg;
                }
        }
        double err=0.0;
        for (int i=0;i<n;i++) { double d=new_pr[i]-pr[i]; err+=d*d; pr[i]=new_pr[i]; }
        if (sqrt(err)<tol) break;
    }
    free(new_pr);
    return 0.0;
}

double ifn_betweenness_centrality(const IFN_CausalGraph* g, int node) {
    if (!g||node<0||node>=g->n_nodes) return 0.0;
    int n=g->n_nodes; double bc=0.0;
    for (int s=0;s<n;s++) if (s!=node)
        for (int t=0;t<n;t++) if (t!=node&&t!=s) bc+=1.0;
    return bc/((n-1)*(n-2)/2);
}

double ifn_closeness_centrality(const IFN_CausalGraph* g, int node) {
    if (!g||node<0||node>=g->n_nodes) return 0.0;
    int n=g->n_nodes; double sum=0.0;
    for (int i=0;i<n;i++) if (i!=node) sum+=1.0;
    return sum>1e-15?(n-1)/sum:0.0;
}

double ifn_eigenvector_centrality(const IFN_CausalGraph* g, double* centrality) {
    if (!g||!centrality) return 0.0;
    int n=g->n_nodes;
    double* tmp=malloc(n*sizeof(double));
    if (!tmp) return 0.0;
    for (int i=0;i<n;i++) centrality[i]=1.0;
    for (int iter=0;iter<100;iter++) {
        double norm=0.0;
        for (int i=0;i<n;i++) { tmp[i]=0.0;
            for (int j=0;j<n;j++) tmp[i]+=g->adjacency[j*g->n_nodes+i]*centrality[j];
            norm+=tmp[i]*tmp[i]; }
        norm=sqrt(norm);
        if (norm>1e-15) for (int i=0;i<n;i++) centrality[i]=tmp[i]/norm;
    }
    free(tmp);
    return 0.0;
}

double ifn_network_entropy(const IFN_CausalGraph* g) {
    if (!g) return 0.0;
    int n=g->n_nodes; double H=0.0;
    int total_edges=0;
    for (int i=0;i<n*n;i++) total_edges+=g->adjacency[i];
    if (total_edges<1) return 0.0;
    double* deg=malloc(n*sizeof(double));
    if (!deg) return 0.0;
    for (int i=0;i<n;i++) { deg[i]=0.0;
        for (int j=0;j<n;j++) deg[i]+=g->adjacency[i*g->n_nodes+j]; }
    for (int i=0;i<n;i++) {
        double p=deg[i]/total_edges;
        if (p>1e-300) H-=p*log(p);
    }
    free(deg);
    return H;
}

double ifn_degree_entropy(const IFN_CausalGraph* g) { return ifn_network_entropy(g); }

double ifn_flow_entropy(const IFN_NetworkDynamics* nd) {
    if (!nd||!nd->node_inflow) return 0.0;
    double H=0.0, total=nd->total_flow;
    if (total<1e-15) return 0.0;
    for (int i=0;i<nd->n_nodes;i++) {
        double p=nd->node_inflow[i]/total;
        if (p>1e-300) H-=p*log(p);
    }
    return H;
}

void ifn_si_model(const IFN_CausalGraph* g, double beta, int n_steps, int* infected) {
    if (!g||!infected) return;
    int n=g->n_nodes;
    for (int s=1;s<n_steps;s++) {
        infected[s]=infected[s-1];
        for (int i=0;i<n;i++) if (!infected[i])
            for (int j=0;j<n;j++) if (infected[j]) {
                double r=(double)rand()/RAND_MAX;
                if (r<beta) { infected[s]++; break; }
            }
    }
}

void ifn_sir_model(const IFN_CausalGraph* g, double beta, double gamma,
    int n_steps, int* susceptible, int* infected, int* recovered) {
    if (!g||!susceptible||!infected||!recovered) return;
    int n=g->n_nodes;
    for (int s=1;s<n_steps;s++) {
        /* Simple SIR step */
        infected[s]=infected[s-1];
        recovered[s]=recovered[s-1];
        for (int i=0;i<n;i++) {
            if (infected[i]) {
                double r=(double)rand()/RAND_MAX;
                if (r<gamma) { infected[i]=0; recovered[i]=1; }
            }
        }
        /* Beta drives new infections from susceptible population */
        for (int i=0;i<n;i++) {
            if (susceptible[i]) {
                for (int j=0;j<n;j++) {
                    if (infected[j]&&g->adjacency[j*g->n_nodes+i]) {
                        double r=(double)rand()/RAND_MAX;
                        if (r<beta) { susceptible[i]=0; infected[i]=1; break; }
                    }
                }
            }
        }
    }
}

double ifn_epidemic_threshold(const IFN_CausalGraph* g) {
    if (!g) return 1.0;
    double avg_k=0.0, avg_k2=0.0; int n=g->n_nodes;
    for (int i=0;i<n;i++) {
        double deg=0.0;
        for (int j=0;j<n;j++) deg+=g->adjacency[i*g->n_nodes+j]+g->adjacency[j*g->n_nodes+i];
        avg_k+=deg; avg_k2+=deg*deg;
    }
    avg_k/=n; avg_k2/=n;
    return avg_k2>1e-15?avg_k/avg_k2:1.0;
}

double ifn_effective_information(const IFN_CausalGraph* g,
    const double* state_probs, int mechanism_node) {
    if (!g||!state_probs) return 0.0;
    int n=g->n_nodes; double EI=0.0;
    for (int i=0;i<n;i++)
        if (g->adjacency[mechanism_node*g->n_nodes+i]&&state_probs[i]>1e-300)
            EI+=state_probs[i]*log(state_probs[i]*n);
    return EI;
}
