#include <iostream>
#include <cstring>

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

int main_() {
    const int Nmon = 6;  // Number of nodes in the graph

    // Connectivity matrix (1D array, flattened version of a 2D matrix)
    // Example graph:
    // 0 - 1 - 2, 3 - 4, 5
    int matrix[Nmon * Nmon] =
    {
        0, 1, 1, 0, 0, 0,  // Node 0 connections
        1, 0, 1, 0, 0, 0,  // Node 1 connections
        1, 1, 0, 0, 0, 0,  // Node 2 connections
        0, 0, 0, 0, 0, 0,  // Node 3 connections
        0, 0, 0, 0, 0, 1,  // Node 4 connections
        0, 0, 0, 0, 1, 0   // Node 5 connections
    };

    // Array to store cluster IDs for each node
    int * clusterID = new int[Nmon];
    memset(clusterID, -1, Nmon*sizeof(int));
    // Find all connected components
    findConnectedComponents(Nmon, matrix, clusterID);

    // Print the cluster assignments
    std::cout << "Node : Cluster ID" << std::endl;
    for (int i = 0; i < Nmon; ++i)
    {
        std::cout << i << " : " << clusterID[i] << std::endl;
    }

    return 0;
}