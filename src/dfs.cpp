#include <iostream>
#ifdef __linux__
#include <cstring>
#endif

// Function to perform Depth-First Search (DFS) using recursion
void dfs(int node, int n, const int* matrix, int* cluster, int currentCluster) {
    // Assign the current cluster ID to the node
    cluster[node] = currentCluster;

    // Explore all neighbors of the node
    for (int i = 0; i < n; ++i) {
        // Check if there's an edge and the node has not been assigned to a cluster
        if (matrix[node * n + i] == 1 && cluster[i] == -1) {
            dfs(i, n, matrix, cluster, currentCluster);
        }
    }
}

void findConnectedComponents(int Nmon, const int* matrix, int* cluster) {
    int currentCluster = 0;  // Cluster ID counter

    // Iterate through all nodes
    for (int i = 0; i < Nmon; ++i) {
        if (cluster[i] == -1)
        {  // If node hasn't been assigned to a cluster
            dfs(i, Nmon, matrix, cluster, currentCluster);  // Perform DFS
            currentCluster++;  // Move to the next cluster ID
        }
    }
}