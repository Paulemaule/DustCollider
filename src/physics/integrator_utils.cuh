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
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_E_s(const double E_i, const double E_j, const double nu_i, const double nu_j) {
    return 1. / (((1 - nu_i * nu_i) / E_i) + ((1 - nu_j * nu_j) / E_j));
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_G_s(const double G_i, const double G_j, const double nu_i, const double nu_j) {
    return 1. / ((1. - nu_i * nu_i) / G_i + (1. - nu_j * nu_j) / G_j);
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_gamma(const double gamma_i, const double gamma_j) {
    return gamma_i + gamma_j - 2.0 / (1.0 / gamma_i + 1.0 / gamma_j);
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_a_0(const double gamma, const double R, const double E_s) {
    return pow(9 * PI * gamma * R * R / E_s, 1.0 / 3.0);
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_delta_N_crit(const double a_0, const double R) {
    return 0.5 * a_0 * a_0 / (R * pow(6.0, 1.0 / 3.0));
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_delta_S_crit(const double nu_i, const double nu_j, const double a_0) {
    return (2.0 - 0.5 * (nu_i + nu_j)) * a_0 / (16.0 * PI);
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_U_N(const double F_c, const double delta_c, const double a, const double a_0) {
    return F_c * delta_c * (0.84661389438303971 + 4. * pow(6., (1. / 3.)) * ((4. / 5.) * pow(a / a_0, 5.) - (4. / 3.) * pow(a / a_0, (7. / 2.)) + (1. / 3.) * pow(a / a_0, 2.)));
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_U_S(const double k_s, const double3 sliding_displacement) {
    return 0.5 * k_s * vec_lenght_sq(sliding_displacement);
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_U_R(const double k_r, const double3 rolling_displacement) {
    return  0.5 * k_r * vec_lenght_sq(rolling_displacement);
}

/**
 * @brief Calculates the
 * 
 * @param
 * @returns The 
 */
__host__ __device__ double get_U_T(const double k_t, const double3 twisting_displacement) {
    return  0.5 * k_t * vec_lenght_sq(twisting_displacement);
}