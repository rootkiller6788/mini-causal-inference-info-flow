#include "ifn_causal.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

IFN_CausalGraph* ifn_causal_create(int n_nodes) {
    if (n_nodes<1) return NULL;
    IFN_CausalGraph* g=calloc(1,sizeof(IFN_CausalGraph));
    if (!g) return NULL;
    g->n_nodes=n_nodes;
    g->adjacency=calloc(n_nodes*n_nodes,sizeof(int));
    g->moral_graph=calloc(n_nodes*n_nodes,sizeof(int));
    g->edge_weights=calloc(n_nodes*n_nodes,sizeof(double));
    if (!g->adjacency||!g->moral_graph||!g->edge_weights) { ifn_causal_free(g); return NULL; }
    return g;
}

void ifn_causal_free(IFN_CausalGraph* g) {
    if (g) { free(g->adjacency); free(g->moral_graph); free(g->edge_weights); free(g->edges); free(g); }
}

void ifn_causal_add_edge(IFN_CausalGraph* g, int from, int to, double weight) {
    if (!g||from<0||from>=g->n_nodes||to<0||to>=g->n_nodes||from==to) return;
    g->adjacency[from*g->n_nodes+to]=1;
    g->edge_weights[from*g->n_nodes+to]=weight;
    g->n_edges++;
    g->edges=realloc(g->edges,g->n_edges*2*sizeof(int));
    if (g->edges) { g->edges[(g->n_edges-1)*2]=from; g->edges[(g->n_edges-1)*2+1]=to; }
}

void ifn_causal_remove_edge(IFN_CausalGraph* g, int from, int to) {
    if (!g||!g->adjacency[from*g->n_nodes+to]) return;
    g->adjacency[from*g->n_nodes+to]=0;
    g->edge_weights[from*g->n_nodes+to]=0.0;
}

bool ifn_causal_has_edge(const IFN_CausalGraph* g, int from, int to) {
    return g&&from>=0&&from<g->n_nodes&&to>=0&&to<g->n_nodes&&g->adjacency[from*g->n_nodes+to];
}

static bool has_cycle_dfs(const IFN_CausalGraph* g, int v, int* visited, int* rec_stack) {
    if (rec_stack[v]) return true;
    if (visited[v]) return false;
    visited[v]=1; rec_stack[v]=1;
    for (int u=0;u<g->n_nodes;u++)
        if (g->adjacency[v*g->n_nodes+u])
            if (has_cycle_dfs(g,u,visited,rec_stack)) return true;
    rec_stack[v]=0;
    return false;
}

bool ifn_is_dag(const IFN_CausalGraph* g) {
    if (!g) return false;
    int* visited=calloc(g->n_nodes,sizeof(int));
    int* rec_stack=calloc(g->n_nodes,sizeof(int));
    if (!visited||!rec_stack) { free(visited); free(rec_stack); return false; }
    for (int v=0;v<g->n_nodes;v++)
        if (!visited[v]&&has_cycle_dfs(g,v,visited,rec_stack)) { free(visited); free(rec_stack); return false; }
    free(visited); free(rec_stack);
    return true;
}

static void topo_dfs(const IFN_CausalGraph* g, int v, int* visited, int* order, int* idx) {
    visited[v]=1;
    for (int u=0;u<g->n_nodes;u++)
        if (g->adjacency[v*g->n_nodes+u]&&!visited[u]) topo_dfs(g,u,visited,order,idx);
    order[(*idx)++]=v;
}

void ifn_topological_sort(const IFN_CausalGraph* g, int* order) {
    if (!g||!order) return;
    int* visited=calloc(g->n_nodes,sizeof(int));
    if (!visited) return;
    int idx=0;
    for (int v=0;v<g->n_nodes;v++)
        if (!visited[v]) topo_dfs(g,v,visited,order,&idx);
    for (int i=0;i<idx/2;i++) { int t=order[i]; order[i]=order[idx-1-i]; order[idx-1-i]=t; }
    free(visited);
}

void ifn_transitive_closure(const IFN_CausalGraph* g, int* closure) {
    if (!g||!closure) return;
    int n=g->n_nodes;
    memcpy(closure,g->adjacency,n*n*sizeof(int));
    for (int k=0;k<n;k++)
        for (int i=0;i<n;i++)
            for (int j=0;j<n;j++)
                if (closure[i*n+k]&&closure[k*n+j]) closure[i*n+j]=1;
}

bool ifn_d_separated(const IFN_CausalGraph* g, int x, int y, const int* z, int n_z) {
    if (!g||x<0||x>=g->n_nodes||y<0||y>=g->n_nodes) return false;
    int* closure=malloc(g->n_nodes*g->n_nodes*sizeof(int));
    if (!closure) return false;
    ifn_moralize(g,g->moral_graph);
    memcpy(closure,g->moral_graph,g->n_nodes*g->n_nodes*sizeof(int));
    for (int k=0;k<g->n_nodes;k++) {
        bool in_z=false;
        for (int i=0;i<n_z;i++) if (k==z[i]) in_z=true;
        if (in_z) for (int i=0;i<g->n_nodes;i++) closure[i*g->n_nodes+k]=closure[k*g->n_nodes+i]=0;
    }
    for (int k=0;k<g->n_nodes;k++) {
        bool in_z=false;
        for (int i=0;i<n_z;i++) if (k==z[i]) in_z=true;
        if (!in_z) for (int i=0;i<g->n_nodes;i++) for (int j=0;j<g->n_nodes;j++)
            if (closure[i*g->n_nodes+k]&&closure[k*g->n_nodes+j]) closure[i*g->n_nodes+j]=1;
    }
    bool sep=!closure[x*g->n_nodes+y];
    free(closure);
    return sep;
}

void ifn_moralize(const IFN_CausalGraph* g, int* moral) {
    if (!g||!moral) return;
    int n=g->n_nodes;
    memcpy(moral,g->adjacency,n*n*sizeof(int));
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        if (g->adjacency[i*n+j]) moral[j*n+i]=1;
    for (int i=0;i<n;i++)
        for (int p1=0;p1<n;p1++)
            for (int p2=p1+1;p2<n;p2++)
                if (g->adjacency[p1*n+i]&&g->adjacency[p2*n+i]) {
                    moral[p1*n+p2]=moral[p2*n+p1]=1;
                }
}

double ifn_causal_flow(const IFN_CausalGraph* g, int source, int target) {
    if (!g||source<0||source>=g->n_nodes||target<0||target>=g->n_nodes) return 0.0;
    double flow=0.0;
    for (int k=0;k<g->n_nodes;k++)
        if (ifn_causal_has_edge(g,k,target)) flow+=g->edge_weights[k*g->n_nodes+target];
    return flow;
}

void ifn_causal_path_effect(const IFN_CausalGraph* g, int source, int target,
    const int* path, int path_len, double* direct, double* indirect) {
    if (!g||!path||path_len<2||source<0||source>=g->n_nodes||target<0||target>=g->n_nodes) return;
    *direct=g->edge_weights[path[0]*g->n_nodes+path[1]];
    *indirect=*direct;
    for (int i=1;i<path_len-1;i++)
        *indirect*=g->edge_weights[path[i]*g->n_nodes+path[i+1]];
}

void ifn_do_intervention(IFN_CausalGraph* g, int node, double value) {
    if (!g||node<0||node>=g->n_nodes) return;
    for (int i=0;i<g->n_nodes;i++) {
        g->adjacency[i*g->n_nodes+node]=0;
        g->edge_weights[i*g->n_nodes+node]=0.0;
    }
    g->edge_weights[node*g->n_nodes+node]=value;
}

double ifn_causal_effect(const IFN_CausalGraph* g, int cause, int effect,
    const double* data, int n) {
    if (!g||!data||n<1||cause<0||cause>=n||effect<0||effect>=g->n_nodes) return 0.0;
    double sum=0.0;
    for (int i=0;i<n;i++) if (i==cause) sum+=data[i];
    return sum/n;
}
