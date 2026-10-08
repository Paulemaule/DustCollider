/**
 * @file integrator_utils.cuh
 * @brief Implements several helper/utility functions used in the physics calculations.
 */

#pragma once

#include <stack>

#include "utils/constant.cuh"
#include "utils/vector.cuh"

/**
 * @brief A macro that calculates the monomer pair indices from the threadID.
 * 
 * This macro determines the layout of monomers pairs in the contact matrices!
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
__host__ __device__ double get_contact_radius(
    const double delta_N,
    const double a_0,
    const double R
) {
    double delta_N_0 = a_0 * a_0 / (3. * R);

    // Substitute the values in the iteration process (a / a_0) =: x, (delta_N / delta_N_0) =: y.
    // The equation that is to be solved becomes: 0 = 3 * x^2 - 2 * sqrt(x) - y.
    double y = delta_N / delta_N_0;
    
    // There is no solution to the equation when the critical displacement is exceeded.
    // Instead the value at the critical displacement is returned.
    double critical_displacement = - pow(9. / 16., 1. / 3.) * delta_N_0;
    if (delta_N <= critical_displacement) {
        return pow(1. / 6., 2. / 3.) * a_0;
    }

    // An initial guess of a = a_0 is used. 
    // This value also ensures that the algorithm converges to the correct solution of the equation.
    double x_n = 1.; 

    // Use Newtons method to determine the solution.
    for (int n = 0; n < 20; n++) {
        // Recursively adjust the guess using the update rule x_n+1 = x_n - f(x_n) / f'(x_n).
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
 * @brief Calculates the combined (reduced) Youngs modulus E* of two monomers.
 *
 * 1 / E* = (1 - nu_i^2) / E_i + (1 - nu_j^2) / E_j
 *
 * @param E_i: Youngs modulus of monomer i.
 * @param E_j: Youngs modulus of monomer j.
 * @param nu_i: Poissons ratio of monomer i.
 * @param nu_j: Poissons ratio of monomer j.
 * @returns The combined Youngs modulus E*.
 */
__host__ __device__ double get_E_s(const double E_i, const double E_j, const double nu_i, const double nu_j) {
    return 1. / (((1 - nu_i * nu_i) / E_i) + ((1 - nu_j * nu_j) / E_j));
}

/**
 * @brief Calculates the combined (reduced) shear modulus G* of two monomers.
 *
 * 1 / G* = (2 - nu_i) / G_i + (2 - nu_j) / G_j
 *
 * @param G_i: Shear modulus of monomer i.
 * @param G_j: Shear modulus of monomer j.
 * @param nu_i: Poissons ratio of monomer i.
 * @param nu_j: Poissons ratio of monomer j.
 * @returns The combined shear modulus G*.
 */
__host__ __device__ double get_G_s(const double G_i, const double G_j, const double nu_i, const double nu_j) {
    return 1. / ((2. - nu_i) / G_i + (2. - nu_j) / G_j);
}

/**
 * @brief Calculates the reduced shear modulus of two monomers.
 *
 * G = G_i * G_j / (G_i + G_j)
 *
 * @param G_i: Shear modulus of monomer i.
 * @param G_j: Shear modulus of monomer j.
 * @returns The reduced shear modulus.
 */
__host__ __device__ double get_G(const double G_i, const double G_j) {
    return G_i * G_j / (G_i + G_j);
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
 * @brief Calculates the critical force at the separation of two monomers.
 *
 * @param gamma: The combined surface energy of the two monomers.
 * @param R: The reduced radius of the two monomers.
 * @returns The critical force.
 */
__host__ __device__ double get_F_c(const double gamma, const double R) {
    return 3. * PI * gamma * R;
}

/**
 * @brief Calculates the strength of the sliding force and torque between two monomers.
 *
 * @param G_s: The combined shear modulus of the two monomers.
 * @param a_0: The equilibrium contact surface radius of the two monomers.
 * @returns The strength of the sliding interaction.
 */
__host__ __device__ double get_k_s(const double G_s, const double a_0) {
    return 8. * G_s * a_0;
}

/**
 * @brief Calculates the strength of the rolling torque between two monomers.
 *
 * @param F_c: The critical force at the separation of the two monomers.
 * @param R: The reduced radius of the two monomers.
 * @returns The strength of the rolling interaction.
 */
__host__ __device__ double get_k_r(const double F_c, const double R) {
    return 4. * F_c / R;
}

/**
 * @brief Calculates the strength of the twisting torque between two monomers.
 *
 * @param G: The reduced shear modulus of the two monomers.
 * @param a_0: The equilibrium contact surface radius of the two monomers.
 * @returns The strength of the twisting interaction.
 */
__host__ __device__ double get_k_t(const double G, const double a_0) {
    return 16. * G * a_0 * a_0 * a_0 / 3.;
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
 * @param k_s: The strength of the sliding interaction between the two monomers.
 * @param sliding_displacement: The sliding displacement of the two monomers.
 * @returns The sliding potential.
 */
__host__ __device__ double get_U_S(const double k_s, const double3 sliding_displacement) {
    return 0.5 * k_s * vec_length_sq(sliding_displacement);
}

/**
 * @brief Calculates the rolling potential between two monomers.
 * 
 * @param k_r: The strength of the rolling interaction between the two monomers.
 * @param rolling_displacement: The rolling displacement of the two monomers.
 * @returns The rolling potential.
 */
__host__ __device__ double get_U_R(const double k_r, const double3 rolling_displacement) {
    return  0.5 * k_r * vec_length_sq(rolling_displacement);
}

/**
 * @brief Calculates the twisting potential between two monomers.
 * 
 * @param k_t: The strength of the twisting interaction between the two monomers.
 * @param twisting_displacement: The twisting displacement of the two monomers. This is conceptually different from Wada (2007). Only the intergated part of eq. (24) is included here.
 * @returns The twisting potential.
 */
__host__ __device__ double get_U_T(const double k_t, const double twisting_displacement) {
    return  0.5 * k_t * twisting_displacement * twisting_displacement;
}

/**
 * @brief Calculates the normal displacement (the compression) of two monomers.
 *
 * @param position_i: The position of monomer i.
 * @param position_j: The position of monomer j.
 * @param r_i: The radius of monomer i.
 * @param r_j: The radius of monomer j.
 * @returns The normal displacement, positive when the monomers overlap.
 */
__host__ __device__ double get_normal_displacement(const double3 position_i, const double3 position_j, const double r_i, const double r_j) {
    return r_i + r_j - vec_dist_len(position_i, position_j);
}

/**
 * @brief Calculates the helper vector of the sliding displacement of two monomers (see Wada et al. 2007).
 *
 * zeta_0 = r_i * n_i - r_j * n_j + (r_i + r_j) * n
 *
 * @param pointer_i: The contact pointer of monomer i in the lab frame, pointing from its center to the contact.
 * @param pointer_j: The contact pointer of monomer j in the lab frame, pointing from its center to the contact.
 * @param normal: The unit vector from monomer j to monomer i, (x_i - x_j) / |x_i - x_j|.
 * @param r_i: The radius of monomer i.
 * @param r_j: The radius of monomer j.
 * @returns The helper vector zeta_0.
 */
__host__ __device__ double3 get_contact_displacement(const double3 pointer_i, const double3 pointer_j, const double3 normal, const double r_i, const double r_j) {
    double3 res;
    res.x = r_i * pointer_i.x - r_j * pointer_j.x + (r_i + r_j) * normal.x;
    res.y = r_i * pointer_i.y - r_j * pointer_j.y + (r_i + r_j) * normal.y;
    res.z = r_i * pointer_i.z - r_j * pointer_j.z + (r_i + r_j) * normal.z;

    return res;
}

/**
 * @brief Calculates the sliding displacement of two monomers, the tangential part of zeta_0.
 *
 * @param contact_displacement: The helper vector zeta_0, see get_contact_displacement.
 * @param normal: The unit vector from monomer j to monomer i, (x_i - x_j) / |x_i - x_j|.
 * @returns The sliding displacement.
 */
__host__ __device__ double3 get_sliding_displacement(const double3 contact_displacement, const double3 normal) {
    double3 res;
    res.x = contact_displacement.x - vec_dot(contact_displacement, normal) * normal.x;
    res.y = contact_displacement.y - vec_dot(contact_displacement, normal) * normal.y;
    res.z = contact_displacement.z - vec_dot(contact_displacement, normal) * normal.z;

    return res;
}

/**
 * @brief Calculates the rolling displacement of two monomers.
 *
 * xi = R * (n_i + n_j)
 *
 * @param pointer_i: The contact pointer of monomer i in the lab frame, pointing from its center to the contact.
 * @param pointer_j: The contact pointer of monomer j in the lab frame, pointing from its center to the contact.
 * @param R: The reduced radius of the two monomers.
 * @returns The rolling displacement.
 */
__host__ __device__ double3 get_rolling_displacement(const double3 pointer_i, const double3 pointer_j, const double R) {
    double3 res;
    res.x = R * (pointer_i.x + pointer_j.x);
    res.y = R * (pointer_i.y + pointer_j.y);
    res.z = R * (pointer_i.z + pointer_j.z);

    return res;
}

/**
 * @brief Assigns a cluster ID to all monomers reachable from `seed` via the contact graph.
 *
 * Uses an iterative DFS to avoid stack overflow on large or deeply connected aggregates.
 *
 * @param seed           Starting monomer index.
 * @param Nmon           Total number of monomers.
 * @param contact_pointer  Nmon×Nmon contact-pointer matrix; entry [i*Nmon+j] is non-zero iff i and j are in contact.
 * @param cluster        Cluster-ID array; unvisited entries must be -1 on entry.
 * @param cluster_id     ID to assign to all reachable monomers.
 */
static void dfs_iterative(const int seed, const int Nmon, const double3* contact_pointer, int* cluster, const int cluster_id) {
    std::stack<int> stack;
    stack.push(seed);
    cluster[seed] = cluster_id;

    while (!stack.empty()) {
        const int node = stack.top();
        stack.pop();

        for (int i = 0; i < Nmon; i++) {
            if (cluster[i] == -1 && vec_length_sq(contact_pointer[node * Nmon + i]) != 0.0) {
                cluster[i] = cluster_id;
                stack.push(i);
            }
        }
    }
}

/**
 * @brief Labels each monomer with the ID of the connected cluster it belongs to.
 *
 * Treats the contact_pointer matrix as an adjacency matrix and performs one iterative
 * DFS per unvisited monomer. Monomers in the same connected component receive the
 * same cluster ID (0-indexed, assigned in discovery order).
 *
 * @param Nmon             Total number of monomers.
 * @param contact_pointer  Nmon×Nmon contact-pointer matrix (host pointer).
 * @param cluster          Output array of length Nmon; must be pre-filled with -1.
 */
void findMonomerClusters(const int Nmon, const double3* contact_pointer, int* cluster) {
    int current_id = 0;

    for (int i = 0; i < Nmon; i++) {
        if (cluster[i] == -1) {
            dfs_iterative(i, Nmon, contact_pointer, cluster, current_id);
            current_id++;
        }
    }
}