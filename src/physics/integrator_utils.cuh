/**
 * @file integrator_utils.cuh
 * @brief Implements several helper/utility functions used in the integrator.
 */

#pragma once

/**
 * @brief A makro that calculates the monomer pair indices from the threadID.
 * 
 * This makro determines the layout of monomers pairs in the contact matrices!
 * The layout in the matrices is: M = {ij} = {00, 10, ..., N0, 01, 11, ..., N1, ..., N-1N, NN}
 */
#define CALC_MONOMER_INDICES(threadID, i, j, matrix_i, matrix_j, Nmon)  \
    i = threadID % Nmon;                                                \
    j = threadID / Nmon;                                                \
    matrix_i = i + j * Nmon;                                            \
    matrix_j = j + i * Nmon;

    
/**
 * @brief Calculates the contact surface radius between two monomers.
 * 
 * This function uses Newtons method to calculate the contact radius a
 * by solving the equation (see Wada et al. - 2007; eq (4))
 * normal_displacement / equilibrium_displacement = 3 * (a / a_0)^2 - 2 * sqrt(a / a_0)
 * 
 * This equation does not have a solution when normal_displacement < -critical_displacement.
 * This case is handled by returning a fixed value (the value of a at the critical diplacement) instead.
 * 
 * @param normal_displacement: The current normal displacement.
 * @param a_0: The equilibrium contact surface radius.
 * @param reduced_radius: The reduced radius of the two monomers.
 * 
 * @returns a: The contact surface radius.
 */
__device__ double get_contact_radius(
    const double delta_N,
    const double a_0,
    const double R
) {
    double delta_N_0 = a_0 * a_0 / (3. * R);

    // Substitute the values in the iteration process (a / a_0) =: x, (delta_N / delta_N_0) =: y.
    // The equation that is to be solved becomes: 0 = 3 * x^2 - 2 * sqrt(x) - y.
    double y = delta_N / delta_N_0;
    
    // There is no solution to the equation when the critical displacement is exeeded.
    // Instead the value at the critical displacement is returned.
    double critical_displacement = - pow(9. / 16., 2. / 3.) * delta_N_0;
    if (delta_N <= critical_displacement) {
        return pow(1. / 6., 2. / 3.) * a_0;
    }

    // An initial guess of a = a_0 is used. 
    // This value also ensures that the algorithm converges to the correct solution of the equation.
    double x_n = 1.; 

    // Use Newtons method to determine the solution.
    for (int n = 0; n < 20; n++) {
        // Recursiveley adjust the guess using the update rule x_n+1 = x_n - f(x_n) / f'(x_n).
        // TODO: Check for optimization opportunities. This piece of code is executed very often.
        x_n = x_n - (3. * x_n * x_n - 2. * sqrt(x_n) - y) / (6. * x_n - 1. / sqrt(x_n));
    }

    // Return the resubstituted root of the equation.
    return x_n * a_0;
}

/**
 * @brief Calculates the reduced radius of two monomers.
 * 
 * @param r_i: The radius of monomer i.
 * @param r_j: The radius of monomer j.
 * @returns The reduced radius.
 */
__host__ __device__ double get_R(const double r_i, const double r_j) {
    return (r_i * r_j) / (r_i + r_j);
}

/**
 * @brief Calculates the shear modulus of a monomer.
 * 
 * @param E_i: Youngs modulus of the monomer.
 * @param nu_i: Poissons ratio of the monomer.
 * @returns The shear modulus.
 */
__host__ __device__ double get_G_i(const double E_i, const double nu_i) {
    return E_i / (2. * (1. + nu_i));
}

/**
 * @brief Calculates the combined shear modulus of two monomers.
 * 
 * @param E_i: The shear modulus of monomer i.
 * @param E_j: The shear modulus of monomer j.
 * @param nu_j: Poissons ratio of monomer i.
 * @param nu_j: Poissons ratio of monomer j.
 * @returns The combined shear modulus.
 */
__host__ __device__ double get_E_s(const double E_i, const double E_j, const double nu_i, const double nu_j) {
    return 1. / (((1 - nu_i * nu_i) / E_i) + ((1 - nu_j * nu_j) / E_j));
}

/**
 * @brief Calculates the combined Youngs modulus of two monomers.
 * 
 * @param G_i: Youngs modulus of monomer i.
 * @param G_j: Youngs modulus of monomer j.
 * @param nu_j: Poissons ratio of monomer i.
 * @param nu_j: Poissons ratio of monomer j.
 * @returns The combined Youngs modulus.
 */
__host__ __device__ double get_G_s(const double G_i, const double G_j, const double nu_i, const double nu_j) {
    return 1. / ((1. - nu_i * nu_i) / G_i + (1. - nu_j * nu_j) / G_j);
}

/**
 * @brief Calculates the combined surface energy of two monomers.
 * 
 * @param gamma_i: The surface energy of monomer i.
 * @param gamma_j: The surface energy of monomer j.
 * @returns The combined surface energy.
 */
__host__ __device__ double get_gamma(const double gamma_i, const double gamma_j) {
    return gamma_i + gamma_j - 2.0 / (1.0 / gamma_i + 1.0 / gamma_j);
}

/**
 * @brief Calculates the equilibrium contact surface radius of two monomers.
 * 
 * @param gamma: The combined surface energy of the two monomers.
 * @param R: The reduced radius of two monomers.
 * @param E_s: The combined shear modulus of the two monomers.
 * @returns The equilibrium contact surface radius.
 */
__host__ __device__ double get_a_0(const double gamma, const double R, const double E_s) {
    return pow(9 * PI * gamma * R * R / E_s, 1.0 / 3.0);
}

/**
 * @brief Calculates the critical normal displacement of two monomers.
 * 
 * @param a_0: The equilibrium contact surface radius of the two monomers.
 * @param R: The reduced radius of the two monomers.
 * @returns The critical normal displacement.
 */
__host__ __device__ double get_delta_N_crit(const double a_0, const double R) {
    return 0.5 * a_0 * a_0 / (R * pow(6.0, 1.0 / 3.0));
}

/**
 * @brief Calculates the critical sliding displacement of two monomers.
 * 
 * @param nu_i: Poissons ratio of monomer i.
 * @param nu_j: Poissons ratio of monomer j.
 * @param a_0: The equilibrium contact surface radius of the two monomers.
 * @returns The critical sliding displacement.
 */
__host__ __device__ double get_delta_S_crit(const double nu_i, const double nu_j, const double a_0) {
    return (2.0 - 0.5 * (nu_i + nu_j)) * a_0 / (16.0 * PI);
}

/**
 * @brief Calculates the normal potential between two monomers.
 * 
 * @param F_c: The normal force at separation of the two monomers.
 * @param delta_N_crit: The critical normal displacement of the two monomers.
 * @param a: The contact surface radius of the two monomers.
 * @param a_0: The equilibrium contact surface radius of the two monomers.
 * @returns The normal potential.
 */
__host__ __device__ double get_U_N(const double F_c, const double delta_N_crit, const double a, const double a_0) {
    return F_c * delta_N_crit * (0.84661389438303971 + 4. * pow(6., (1. / 3.)) * ((4. / 5.) * pow(a / a_0, 5.) - (4. / 3.) * pow(a / a_0, (7. / 2.)) + (1. / 3.) * pow(a / a_0, 2.)));
}

/**
 * @brief Calculates the sliding potential between two monomers.
 * 
 * @param k_s: The strenght of the sliding interaction between the two monomers.
 * @param sliding_displacement: The sliding displacement of the two monomers.
 * @returns The sliding potential.
 */
__host__ __device__ double get_U_S(const double k_s, const double3 sliding_displacement) {
    return 0.5 * k_s * vec_lenght_sq(sliding_displacement);
}

/**
 * @brief Calculates the rolling potential between two monomers.
 * 
 * @param k_r: The strenght of the rolling interaction between the two monomers.
 * @param rolling_displacement: The rolling displacement of the two monomers.
 * @returns The rolling potential.
 */
__host__ __device__ double get_U_R(const double k_r, const double3 rolling_displacement) {
    return  0.5 * k_r * vec_lenght_sq(rolling_displacement);
}

/**
 * @brief Calculates the twisting potential between two monomers.
 * 
 * @param k_t: The strenght of the twisting interaction between the two monomers.
 * @param twisting_displacement: The twisting displacement of the two monomers. This is conceptually different from Wada (2007). Only the intergated part of eq. (24) is included here.
 * @returns The twisting potential.
 */
__host__ __device__ double get_U_T(const double k_t, const double twisting_displacement) {
    return  0.5 * k_t * twisting_displacement * twisting_displacement;
}

/**
 * @brief A function to recursively search for nodes connected to the current node.
 * 
 * This function uses a depth-first search (DFS) algorithm to explore the graph represented by the connection matrix.
 * 
 * @param node: The current node being explored.
 * @param n: The total number of nodes.
 * @param matrix: The connection matrix representing the graph.
 * @param clusters: The array storing the cluster IDs for each node.
 * @param currentCluster: The ID of the current cluster being assigned.
 */
void dfs(const int node, const int n, const double3* matrix, int* clusters, const int currentCluster) {
    // Assign the current cluster ID to the node.
    clusters[node] = currentCluster;

    // Check all nodes to see if they are connected to the current node.
    for (int i = 0; i < n; i++) {
        // Check if there's an edge and the node has not been assigned to a cluster
        if (vec_lenght_sq(matrix[node * n + i]) != 0. && clusters[i] == -1) {
            dfs(i, n, matrix, clusters, currentCluster);
        }
    }
}

/**
 * @brief A function to find distinct dust aggregates using the contact pointer matrix as an adjacency matrix.
 * 
 * This function iterates through all nodes in the graph and uses a depth-first search (DFS) algorithm to find connected monomers.
 * 
 * @param Nmon: The number of nodes in the graph.
 * @param contact_pointer: The matrix of current contact pointers.
 * @param cluster: The array storing the cluster IDs for each node.
 */
void findMonomerClusters(const int Nmon, const double3* contact_pointer, int* cluster) {
    int currentCluster = 0;  // Cluster ID counter

    // Iterate through all nodes
    for (int i = 0; i < Nmon; ++i) {
        if (cluster[i] == -1)
        {  // If node hasn't been assigned to a cluster
            dfs(i, Nmon, contact_pointer, cluster, currentCluster);
            currentCluster++;  // Move to the next cluster ID
        }
    }
}