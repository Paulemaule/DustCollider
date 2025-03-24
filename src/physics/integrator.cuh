/**
 * @file integrator.cuh
 * @brief Implements a synchronized leapfrog algorithm for time integration.
 * 
 * This file defines CUDA kernels for performing a synchronized leapfrog integration 
 * scheme, specifically using the synchronized Predictor-Evaluate-Corrector (PEC) method. 
 * 
 * The integration scheme updates the state according to this formula:
 * x(t+dt) : Predictor
 * x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
 * 
 * a(t+dt) : Evaluator
 * a(t+dt) = f(x(t+dt))
 * 
 * v(t+dt) : Corrector
 * v(t+dt) = v(t) + 0.5 * (a(t) + a(t+dt)) * dt
 * 
 * The integration scheme is implemented in the corresponding CUDA Kernels.
 */

#pragma once

#include "utils/vector.cuh"
#include "utils/typedefs.cuh"
#include "utils/constant.cuh"

#include "physics/integrator_utils.cuh"

/**
 * @brief Implements the predictor of the synchronized leapfrog algorithm for the positions of the monomers.
 * 
 * Calculates the positions of the monomers for the next timestep according to the formula
 * x(t+dt) = x(t) + dt * v(t) + 0.5 * dt * dt * a(t)
 * 
 * @param position_curr: The current position x(t) of the monomers.
 * @param velocity_curr: The current velocities v(t) of the monomers.
 * @param force_curr: The current forces F(t) acting on the monomers.
 * @param position_next: The location where the next positions x(t+dt) will be written to.
 * @param mass: The masses of the monomers.
 * @param timestep: The timestep of the simulation.
 * @param Nmon: The number of monomers.
 */
__global__ void predictor(
    const double3*              position_curr,
    const double3*              velocity_curr,
    const double3*              force_curr,

    double3*                    position_next,
    
    const double*               mass,
    const double                timestep,
    const double                Nmon
) {
    // Retrieve the ID of the current thread.
    int threadID = blockDim.x * blockIdx.x + threadIdx.x; // The ID of the current thread.

    // Terminate threads that do not correspond to a monomer.
    if (threadID >= Nmon) return;

    double inv_mass = 1.0 / mass[threadID];

    position_next[threadID].x = position_curr[threadID].x + timestep * velocity_curr[threadID].x + 0.5 * timestep * timestep * inv_mass * force_curr[threadID].x;
    position_next[threadID].y = position_curr[threadID].y + timestep * velocity_curr[threadID].y + 0.5 * timestep * timestep * inv_mass * force_curr[threadID].y;
    position_next[threadID].z = position_curr[threadID].z + timestep * velocity_curr[threadID].z + 0.5 * timestep * timestep * inv_mass * force_curr[threadID].z;
}

/**
 * @brief Implements the predictor of the synchronized leapfrog algorithm for the rotation of the monomers.
 * 
 * Calculates the rotations of the monomers contacts for the next timestep according to the formula
 * x(t+dt) = x(t) + dt * v(t) + 0.5 * dt * dt * a(t)
 */
__global__ void predictor_pointer(
    const double4*              rotation_curr,
    const double*               twisting_curr,
    const double3*              omega_curr,
    const double3*              torque_curr,
    
    double4*                    rotation_next,
    double*                     twisting_next,

    const double*               moment,
    const double                timestep,
    const int                   Nmon
) {
    // Retrieve the ID of the current thread
    int threadID = blockDim.x * blockIdx.x + threadIdx.x; // The ID of the current thread.
    
    // Terminate threads that do not correspond to a monomer pair.
    if (threadID >= (Nmon * Nmon)) return;

    // Index of the monomer in the arrays.
    int i, j;
    // Index of the monomer in the pair matrices.
    int matrix_i, matrix_j;
    
    CALC_MONOMER_INDICES(threadID, i, j, matrix_i, matrix_j, Nmon);

    // Skip threads that dont correspond to a monomer pair.
    if (i == j) return;

    // Integrate contact rotation matrix
    double4 rot = rotation_curr[matrix_i];
    double3 omega = omega_curr[i];

    double moment_i = moment[i];
    double3 omega_dot = torque_curr[i];
    omega_dot.x /= moment_i;
    omega_dot.y /= moment_i;
    omega_dot.z /= moment_i;

    // CHECK: Is this formula correct?
    double4 e_dot, e_ddot;
    e_dot.w = - 0.5 * (rot.x * omega.x + rot.y * omega.y + rot.z * omega.z);
    e_dot.x = 0.5 * (rot.w * omega.x - rot.y * omega.z + rot.z * omega.y);
    e_dot.y = 0.5 * (rot.w * omega.y - rot.z * omega.x + rot.x * omega.z);
    e_dot.z = 0.5 * (rot.w * omega.z - rot.x * omega.y + rot.y * omega.x);

    double temp = 0.5 * e_dot.w;

    e_ddot.w = - 0.25 * (rot.w * vec_lenght_sq(omega) + 2.0 * (rot.x * omega_dot.x + rot.y * omega_dot.y + rot.z * omega_dot.z));
    e_ddot.x = temp * omega.x + 0.5 * (rot.w * omega_dot.x - rot.y * omega_dot.z + rot.z * omega_dot.y);
    e_ddot.y = temp * omega.y + 0.5 * (rot.w * omega_dot.y - rot.z * omega_dot.x + rot.x * omega_dot.z);
    e_ddot.z = temp * omega.z + 0.5 * (rot.w * omega_dot.z - rot.x * omega_dot.y + rot.y * omega_dot.x);

    rot.w = rot.w + timestep * e_dot.w + 0.5 * timestep * timestep * e_ddot.w;
    rot.x = rot.x + timestep * e_dot.x + 0.5 * timestep * timestep * e_ddot.x;
    rot.y = rot.y + timestep * e_dot.y + 0.5 * timestep * timestep * e_ddot.y;
    rot.z = rot.z + timestep * e_dot.z + 0.5 * timestep * timestep * e_ddot.z;

    // Renormalize the rotation quaternion.
    quat_normalize(rot);

    // Apply the changes.
    rotation_next[matrix_i].w = rot.w;
    rotation_next[matrix_i].x = rot.x;
    rotation_next[matrix_i].y = rot.y;
    rotation_next[matrix_i].z = rot.z;

    // Integrate contact twisting matrix
    // FIXME: Finish this implementation.    
    twisting_next[matrix_i] = twisting_curr[matrix_i];
};

/**
 * @brief Impements the corrector of the synchronized leapfrog algorithm.
 * 
 * Calculates the velocities and angular velocities of the monomers for the next timestep according to the formula
 * v(t+dt) = v(t) + 0.5 * dt * (a(t) + a(t+dt))
 * 
 * @param velocity_curr: The current velocity v(t) of the monomers.
 * @param omega_curr: The current angular velocities w(t) of the monomers.
 * @param force_curr: The current forces F(t) acting on the monomers.
 * @param force_next: The next forces F(t+dt) acting on the monomers.
 * @param torque_curr: The current torques M(t) acting on the monomers.
 * @param torque_next: The next torques M(t+dt) acting on the monomers.
 * @param velocity_next: The location where the next velocities v(t+dt) will be written to.
 * @param omega_next: The location where the next angular velocities w(t+dt) will be written to.
 * @param mass: The masses of the monomers.
 * @param moment: The angular momenta of the monomers.
 * @param timestep: The timestep of the simulation.
 * @param Nmon: The number of monomers.
 */
__global__ void corrector(
    const double3*              velocity_curr,
    const double3*              omega_curr,
    const double3*              force_curr,
    const double3*              force_next,
    const double3*              torque_curr,
    const double3*              torque_next,

    double3*                    velocity_next,
    double3*                    omega_next,

    const double*               mass,
    const double*               moment,
    const double                timestep,
    const int                   Nmon
) {
    // Retrieve the ID of the current thread
    int threadID = blockDim.x * blockIdx.x + threadIdx.x; // The ID of the current thread.

    // Terminate threads that do not correspond to a monomer.
    if (threadID >= Nmon) return;

    double inv_mass = 1.0 / mass[threadID];
    double inv_moment = 1.0 / moment[threadID];

    // Update the monomer velocity
    velocity_next[threadID].x = velocity_curr[threadID].x + 0.5 * timestep * inv_mass * (force_curr[threadID].x + force_next[threadID].x);
    velocity_next[threadID].y = velocity_curr[threadID].y + 0.5 * timestep * inv_mass * (force_curr[threadID].y + force_next[threadID].y);
    velocity_next[threadID].z = velocity_curr[threadID].z + 0.5 * timestep * inv_mass * (force_curr[threadID].z + force_next[threadID].z);

    // Update the monomer angular velocity
    omega_next[threadID].x = omega_curr[threadID].x + 0.5 * inv_moment * timestep * (torque_curr[threadID].x + torque_next[threadID].x);
    omega_next[threadID].y = omega_curr[threadID].y + 0.5 * inv_moment * timestep * (torque_curr[threadID].y + torque_next[threadID].y);
    omega_next[threadID].z = omega_curr[threadID].z + 0.5 * inv_moment * timestep * (torque_curr[threadID].z + torque_next[threadID].z);
}

// TODO: The collaborator list is incomplete.
/**
 * @brief Implements the evaluation step of the synchronized leapfrog algorithm.
 * 
 * Calculates the forces and torques acting on the monomers according to the interaction model 
 * developed by Johnson, Kendall and Roberts (the JKR model) with further extensions by 
 * Dominik and Nübold (2002), Wada et al. (2007) and (?)
 * 
 * The JKR model is used to calculate the inter monomer forces.
 * For this purpose the contact pointer approach described in Dominik and Nübuld (2002) as well as Wada et al (2007) is implemented.
 * The damping force proposed by (?) is included for (?).
 */
__global__ void evaluate(
    const double3*              position_next,
    const double3*              pointer_curr,
    const double4*              rotation_next,
    const double*               twisting_next,
    const double*               compression_old,

    double3*                    force_next,
    double3*                    torque_next,
    double4*                    potential_energy,
    double4*                    inelastic_counter,

    const double*               mass,
    const double*               radius,
    const double*               youngs_modulus,
    const double*               poisson_number,
    const double*               surface_energy,
    const double*               crit_rolling_displacement,
    const double*               viscous_damping_timescale,
    const double                timestep,
    const int                   Nmon
) {
    // Retrieve the ID of the current thread
    int threadID = blockDim.x * blockIdx.x + threadIdx.x; // The ID of the current thread.
    
    // Terminate threads that do not correspond to a monomer pair.
    if (threadID >= (Nmon * Nmon)) return;

    // Index of the monomer in the arrays.
    int i, j;
    // Index of the monomer in the pair matrices.
    int matrix_i, matrix_j;
    
    CALC_MONOMER_INDICES(threadID, i, j, matrix_i, matrix_j, Nmon);

    // This branching is unnecessary because this case is filtered out later when the contacts are checked.
    // The case is then skipped because the pair i-i is allways marked as unconnected.
    //if (i == j) return; 

    // Initialize containers for the force and torque.
    double3 force = { 0.0, 0.0, 0.0 };              // A container for the total force resulting from the pair interaction.
    double3 torque = { 0.0, 0.0, 0.0 };             // A container for the total torque resulting from the pair interaction.

    // DETERMINE PAIR PROPERTIES
    double r_i = radius[i];                         // The radius of monomer i.
    double r_j = radius[j];                         // The radius of monomer j.

    double E_i = youngs_modulus[i];                 // Youngs modulus of monomer i.
    double E_j = youngs_modulus[j];                 // Youngs modulus of monomer j.

    double nu_i = poisson_number[i];                // Poisson number of monomer i.
    double nu_j = poisson_number[j];                // Poisson number of monomer j.

    double gamma_i = surface_energy[i];             // The surface energy of monomer i.
    double gamma_j = surface_energy[j];             // The surface energy of monomer j.

    double3 pointer_i = pointer_curr[matrix_i];     // The contact pointer, pointing from the center of monomer i to the contact location with monomer j.
    double3 pointer_j = pointer_curr[matrix_j];     // The contact pointer, pointing from the center of monomer j to the contact location with monomer i.

    double3 position_i = position_next[i];          // The position of monomer i.
    double3 position_j = position_next[j];          // The position of monomer j.

    double G_i = get_G_i(E_i, nu_i);          // The shear modulus of monomer i.
    double G_j = get_G_i(E_j, nu_j);          // The shear modulus of monomer j.

    // The reduced radius of the monomer pair.
    double R = get_R(r_i, r_j);

    // The combined Youngs modulus of the monomer pair.
    double E_s = get_E_s(E_i, E_j, nu_i, nu_j);

    // The combined shear modulus of the monomer pair.
    double G_s = get_G_s(G_i, G_j, nu_i, nu_j);

    // The reduced shear modulus of the monomer pair.
    double G = G_i * G_j / (G_i + G_j);

    // The surface energy of the monomer pair.
    double gamma = get_gamma(gamma_i, gamma_j);

    // The equilibrium contact radius of the monomer pair.
    double a_0 = get_a_0(gamma, R, E_s);

    // The force and torque strenghts.
    double F_c = 3. * PI * gamma * R;                   // The critical force at monomer separation.
    double k_s = 8. * G_s * a_0;                        // The strenght of the sliding force and torque.
    double k_r = 4. * F_c / R;                          // The strenght of the rolling torque.
    double k_t = 16. * G * a_0 * a_0 * a_0 / 3.;        // The strenght of the twisting torque.

    // The critical displacements.
    double delta_N_crit = get_delta_N_crit(a_0, R);

    // The viscous damping timescale of the pair.
    double t_vis = 0.5 * (viscous_damping_timescale[i] + viscous_damping_timescale[j]);

    // Calculate contact effects only if the monomers are in contact.
    if (vec_lenght(pointer_i) != 0. && vec_lenght(pointer_j) != 0.) {
        // Corotate the contact pointers
        pointer_i = quat_apply_inverse(rotation_next[matrix_i], pointer_i);
        pointer_j = quat_apply_inverse(rotation_next[matrix_j], pointer_j);

        // The pointer from the current position of monomer i to the current position of monomer j.
        double3 pointer_pos = vec_get_normal(position_i, position_j);

        // Caculate the displacements
        double normal_displacement;             // The displacement in the normal-dof of the contact.
        double3 sliding_displacement;           // The displacement in the sliding-dof of the contact.
        double3 rolling_displacement;           // The displacement in the rolling-dof of the contact.
        double twisting_displacement;           // The displacement in the twisting-dof of the contact.

        normal_displacement = r_i + r_j - vec_dist_len(position_i, position_j);
        
        double3 contact_displacement;           // A helper variable in the displacement calculation (zeta_0 in Wada07).
        contact_displacement.x = r_i * pointer_i.x - r_j * pointer_j.x + (r_i + r_j) * pointer_pos.x;
        contact_displacement.y = r_i * pointer_i.y - r_j * pointer_j.y + (r_i + r_j) * pointer_pos.y;
        contact_displacement.z = r_i * pointer_i.z - r_j * pointer_j.z + (r_i + r_j) * pointer_pos.z;

        sliding_displacement.x = contact_displacement.x - vec_dot(contact_displacement, pointer_pos) * pointer_pos.x;
        sliding_displacement.y = contact_displacement.y - vec_dot(contact_displacement, pointer_pos) * pointer_pos.y;
        sliding_displacement.z = contact_displacement.z - vec_dot(contact_displacement, pointer_pos) * pointer_pos.z;

        rolling_displacement.x = R * (pointer_i.x + pointer_j.x);
        rolling_displacement.y = R * (pointer_i.y + pointer_j.y);
        rolling_displacement.z = R * (pointer_i.z + pointer_j.z);

        twisting_displacement = twisting_next[matrix_i];

        // FORCE AND TORQUE CALCULATIONS
        // Normal
        double a = get_contact_radius(normal_displacement, a_0, R);
        double normal_force = 4 * F_c * (pow((a / a_0), 3.) - pow((a / a_0), 1.5));

        force.x += normal_force * pointer_pos.x;
        force.y += normal_force * pointer_pos.y;
        force.z += normal_force * pointer_pos.z;

        // Damping
        double vis_damping_strenght = 2.0 * t_vis / (nu_i * nu_j) * E_s;
        double delta_N_dot = (normal_displacement - compression_old[matrix_i]) / timestep;
        double damping_force = vis_damping_strenght * a * delta_N_dot;

        force.x += damping_force * pointer_pos.x;
        force.y += damping_force * pointer_pos.y;
        force.z += damping_force * pointer_pos.z;

        // Sliding
        double particle_distance = vec_dist_len(position_i, position_j);
        double sliding_force = - k_s * (r_i + r_j - vec_dot(contact_displacement, pointer_pos)) / particle_distance;

        force.x += sliding_force * sliding_displacement.x;
        force.y += sliding_force * sliding_displacement.y;
        force.z += sliding_force * sliding_displacement.z;

        double3 tmp_s = vec_cross(pointer_i, sliding_displacement);

        torque.x += - r_i * k_s * tmp_s.x;
        torque.y += - r_i * k_s * tmp_s.y;
        torque.z += - r_i * k_s * tmp_s.z;

        // Rolling
        double3 tmp_r = vec_cross(pointer_i, rolling_displacement);

        torque.x += - R * k_r * tmp_r.x;
        torque.y += - R * k_r * tmp_r.y;
        torque.z += - R * k_r * tmp_r.z;

        // Twisting
        //torque.x += - k_t * twisting_displacement.x;
        //torque.y += - k_t * twisting_displacement.y;
        //torque.z += - k_t * twisting_displacement.z;

        // Add the forces and torques from the pair interaction onto the total for monomer i.
        atomicAdd(&force_next[i].x, force.x);
        atomicAdd(&force_next[i].y, force.y);
        atomicAdd(&force_next[i].z, force.z);
        
        atomicAdd(&torque_next[i].x, torque.x);
        atomicAdd(&torque_next[i].y, torque.y);
        atomicAdd(&torque_next[i].z, torque.z);

        // Add energy dissipation.
        atomicAdd(&inelastic_counter->w, 0.5 * damping_force * (normal_displacement - compression_old[matrix_i]));

        // Add the potential energies.
        atomicAdd(&potential_energy->w, 0.5 * get_U_N(F_c, delta_N_crit, a, a_0));
        atomicAdd(&potential_energy->x, 0.5 * get_U_S(k_s, sliding_displacement));
        atomicAdd(&potential_energy->y, 0.5 * get_U_R(k_r, rolling_displacement));
        atomicAdd(&potential_energy->z, 0.5 * get_U_T(k_t, { 0., 0., 0. }));
    }
}

/**
 * @brief Checks for inelastic motion and updates the contact pointers accordingly.
 * 
 * This function checks for critical displacements in all degrees of freedom.
 * When critical displacement is encountered the contact pointers are adjusted.
 */
__global__ void updatePointers(
    const double3*              position_next,
    const double3*              pointer_curr,
    const double4*              rotation_curr, // FIXME: Check if this leads to problems. Conceptually rotation_next should be used as it is the current contact rotation, but this would lead to memory races.
    const double*               compression_curr, // TODO: This is never used in this function

    double3*                    pointer_next,
    double4*                    rotation_next,
    double*                     twisting_next,
    double*                     compression_next,
    double4*                    inelastic_counter,

    const double*               radius,
    const double*               youngs_modulus,
    const double*               poisson_number,
    const double*               surface_energy,
    const double*               crit_rolling_displacement,
    const int                   Nmon
) {
    // Retrieve the ID of the current thread
    int threadID = blockDim.x * blockIdx.x + threadIdx.x; // The ID of the current thread.
    
    // Terminate threads that do not correspond to a monomer pair.
    if (threadID >= (Nmon * Nmon)) return;

    // Index of the monomer in the arrays.
    int i, j;
    // Index of the monomer in the pair matrices.
    int matrix_i, matrix_j;
    
    CALC_MONOMER_INDICES(threadID, i, j, matrix_i, matrix_j, Nmon);

    if (i == j) return;                             // There is no self interaction
    //if (i > j) return;                              // The contact pointer updates need to be symmetrical

    // DETERMINE PAIR PROPERTIES
    double r_i = radius[i];                         // The radius of monomer i.
    double r_j = radius[j];                         // The radius of monomer j.

    double E_i = youngs_modulus[i];                 // Youngs modulus of monomer i.
    double E_j = youngs_modulus[j];                 // Youngs modulus of monomer j.

    double nu_i = poisson_number[i];                // Poisson number of monomer i.
    double nu_j = poisson_number[j];                // Poisson number of monomer j.

    double gamma_i = surface_energy[i];             // The surface energy of monomer i.
    double gamma_j = surface_energy[j];             // The surface energy of monomer j.

    double3 pointer_i = pointer_curr[matrix_i];     // The contact pointer, pointing from the center of monomer i to the contact location with monomer j.
    double3 pointer_j = pointer_curr[matrix_j];     // The contact pointer, pointing from the center of monomer j to the contact location with monomer i.

    double3 position_i = position_next[i];          // The position of monomer i.
    double3 position_j = position_next[j];          // The position of monomer j.

    double G_i = get_G_i(E_i, nu_i);                // The shear modulus of monomer i.
    double G_j = get_G_i(E_j, nu_j);                // The shear modulus of monomer j.

    // The reduced radius of the monomer pair.
    double R = get_R(r_i, r_j);

    // The combined Youngs modulus of the monomer pair.
    double E_s = get_E_s(E_i, E_j, nu_i, nu_j);

    // The combined shear modulus of the monomer pair.
    double G_s = get_G_s(G_i, G_j, nu_i, nu_j);

    // The reduced shear modulus of the monomer pair.
    double G = G_i * G_j / (G_i + G_j);

    // The surface energy of the monomer pair.
    double gamma = get_gamma(gamma_i, gamma_j);

    // The equilibrium contact radius of the monomer pair.
    double a_0 = get_a_0(gamma, R, E_s);

    // The force and torque strenghts.
    double F_c = 3. * PI * gamma * R;                           // The critical force at monomer separation.
    double k_s = 8. * G_s * a_0;                                // The strenght of the sliding force and torque.
    double k_r = 4. * F_c / R;                                  // The strenght of the rolling torque.
    double k_t = 16. * G * a_0 * a_0 * a_0 / 3.;                // The strenght of the twisting torque.

    // The critical displacements.
    double delta_N_crit = get_delta_N_crit(a_0, R);             // The critical normal displacement of the monomer pair.
    // FIXME: In wada 2007 it is mentioned that this formula is not valid for materials with strong intermolecular forces.
    double delta_S_crit = get_delta_S_crit(nu_i, nu_j, a_0);    // The critical sliding displacement of the monomer pair.
    double delta_R_crit = 0.5 * (crit_rolling_displacement[i] + crit_rolling_displacement[j]);  // The critical rolling displacement of the monomer pair.
    double delta_T_crit = 1. / (16. * PI);                      // The critical twisting displacement of the monomer pair.

    if (vec_lenght_sq(pointer_i) != 0. && vec_lenght_sq(pointer_j) != 0.) {
        // Corotate the contact pointers.
        double4 rotation_i = rotation_curr[matrix_i];
        double4 rotation_j = rotation_curr[matrix_j];
        pointer_i = quat_apply_inverse(rotation_i, pointer_i);
        pointer_j = quat_apply_inverse(rotation_j, pointer_j);

        double3 pointer_pos = vec_get_normal(position_i, position_j); // The current pointer from monomer j to monomer i.

        // Caculate the displacements
        double normal_displacement;             // The displacement in the normal-dof of the contact.
        double3 sliding_displacement;           // The displacement in the sliding-dof of the contact.
        double3 rolling_displacement;           // The displacement in the rolling-dof of the contact.
        double twisting_displacement;           // The displacement in the twisting-dof of the contact.

        normal_displacement = r_i + r_j - vec_dist_len(position_i, position_j);
        
        double3 contact_displacement;           // A helper variable in the displacement calculation.
        contact_displacement.x = r_i * pointer_i.x - r_j * pointer_j.x + (r_i + r_j) * pointer_pos.x;
        contact_displacement.y = r_i * pointer_i.y - r_j * pointer_j.y + (r_i + r_j) * pointer_pos.y;
        contact_displacement.z = r_i * pointer_i.z - r_j * pointer_j.z + (r_i + r_j) * pointer_pos.z;

        // ASK: This differs significantly from Stefans implementation, why?
        sliding_displacement.x = contact_displacement.x - vec_dot(contact_displacement, pointer_pos) * pointer_pos.x;
        sliding_displacement.y = contact_displacement.y - vec_dot(contact_displacement, pointer_pos) * pointer_pos.y;
        sliding_displacement.z = contact_displacement.z - vec_dot(contact_displacement, pointer_pos) * pointer_pos.z;

        rolling_displacement.x = R * (pointer_i.x + pointer_j.x);
        rolling_displacement.y = R * (pointer_i.y + pointer_j.y);
        rolling_displacement.z = R * (pointer_i.z + pointer_j.z);

        twisting_displacement = twisting_next[matrix_i];

        // Check for critical displacements
        if (- normal_displacement > delta_N_crit) {
            // The monomers loose contact
            pointer_next[matrix_i]      = { 0., 0., 0. };
            rotation_next[matrix_i]     = { 0., 0., 0., 0. };
            twisting_next[matrix_i]     = 0.;
            compression_next[matrix_i]  = 0.;

            // All potential energy stored in the connection is lost.
            atomicAdd(&inelastic_counter->w, 0.5 * get_U_N(F_c, delta_N_crit, get_contact_radius(normal_displacement, a_0, R), a_0));
            atomicAdd(&inelastic_counter->x, 0.5 * get_U_S(k_s, sliding_displacement));
            atomicAdd(&inelastic_counter->y, 0.5 * get_U_R(k_r, rolling_displacement));
            atomicAdd(&inelastic_counter->z, 0.5 * get_U_T(k_t, { 0., 0., 0. }));

            return; 
        } else {
            // The monomers stay in contact
            pointer_next[matrix_i] = pointer_curr[matrix_i];
            compression_next[matrix_i] = normal_displacement;
            // rotation_next, twisting_next do not need to be updated here, as they are allready being updated in the pointer_corrector
        }

        // FIXME: This is not perfect. The correction to rolling, can, if inelestic sliding also occurs, be inaccurate because rolling displacement has been calculated with an inaccurate pointer.
        double sliding_displacement_abs = vec_lenght(sliding_displacement);
        double rolling_displacement_abs = vec_lenght(rolling_displacement);

        if (sliding_displacement_abs > delta_S_crit) {
            // Inelastic sliding motion.
            double3 correction;
            correction.x = sliding_displacement.x * (1. - delta_S_crit / sliding_displacement_abs);
            correction.y = sliding_displacement.y * (1. - delta_S_crit / sliding_displacement_abs);
            correction.z = sliding_displacement.z * (1. - delta_S_crit / sliding_displacement_abs);

            // FIXME: Ensure that the calculation of the angle theta_1 (wada07) is correct.
            double correction_factor = vec_dot(pointer_i, correction) / vec_lenght(correction);
            correction_factor = 1. / (1. - correction_factor * correction_factor);
            correction_factor = correction_factor / (2. * r_i);

            // Calculate the corrected pointer
            pointer_i.x -= correction.x * correction_factor;
            pointer_i.y -= correction.y * correction_factor;
            pointer_i.z -= correction.z * correction_factor;

            // Re-normalize the pointer.
            vec_normalize(pointer_i);
            
            // Track dissipated energy.
            atomicAdd(&inelastic_counter->x, 0.5 * k_s * delta_S_crit * (sliding_displacement_abs - delta_S_crit));
        }

        if (rolling_displacement_abs > delta_R_crit) {
            // Inelastic rolling motion.
            double3 correction;
            correction.x = rolling_displacement.x * (1. - delta_R_crit / rolling_displacement_abs);
            correction.y = rolling_displacement.y * (1. - delta_R_crit / rolling_displacement_abs);
            correction.z = rolling_displacement.z * (1. - delta_R_crit / rolling_displacement_abs);

            double correction_factor = vec_dot(pointer_i, correction) / vec_lenght(correction);
            correction_factor = 1. / (1. - correction_factor * correction_factor);
            correction_factor = correction_factor / (2. * r_i);
            
            // Calculate the corrected pointer
            pointer_i.x -= correction.x * correction_factor;
            pointer_i.y -= correction.y * correction_factor;
            pointer_i.z -= correction.z * correction_factor;
 
            // Re-normalize the pointer.
            vec_normalize(pointer_i);

            // Track dissipated energy.
            atomicAdd(&inelastic_counter->y, 0.5 * k_r * delta_R_crit * (rolling_displacement_abs - delta_R_crit));
        }

        // If there were any corrections, apply them to the contact pointer.
        if (sliding_displacement_abs > delta_S_crit || rolling_displacement_abs > delta_R_crit) {
            // Corotate the corrected pointer to calculate the pointer in the monomer fixed system.
            pointer_next[matrix_i] = quat_apply(rotation_i, pointer_i);
        }

        if (twisting_displacement * twisting_displacement > delta_T_crit * delta_T_crit) {
            // Inelastic twisting motion.
            int sign = 1 - (2. * signbit(twisting_displacement)); // Extract the sign of the twisting displacement
            twisting_next[matrix_i] = sign * delta_T_crit;
            
            // Track dissipated energy.
            atomicAdd(&inelastic_counter->z, 0.5 * k_t * delta_T_crit * (twisting_displacement - delta_T_crit));
        }
    } else {
        double normal_displacement;             // The displacement in the normal-dof of the contact.
        normal_displacement = r_i + r_j - vec_dist_len(position_i, position_j);

        // Check if the monomer pair is currently making contact
        if (normal_displacement >= 0.) {
            // The monomers are touching -> initialize the contact pointer.
            pointer_next[matrix_i] = vec_get_normal(position_j, position_i);
            rotation_next[matrix_i] = { 0., 0., 0., 1. };
            compression_next[matrix_i] = normal_displacement;

            atomicAdd(&inelastic_counter->w, - 0.5 * get_U_N(F_c, delta_N_crit, get_contact_radius(normal_displacement, a_0, R), a_0));
        } else {
            pointer_next[matrix_i] = { 0., 0., 0. };
            rotation_next[matrix_i] = { 0., 0., 0., 0. };
            twisting_next[matrix_i] = 0.;
            compression_next[matrix_i] = normal_displacement;
        }
    }
}