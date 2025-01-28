#pragma once

#include "utils/constant.cuh"
#include "utils/vector.cuh"
#include "utils/typedefs.cuh"

__device__ void updateNeighbourhoodRelations();
__device__ void updateContacts();
__device__ void updateParticleInteraction();

/**
 * The evaluate step in PECE scheme
 */
__global__ void evaluate() {
    updateNeighbourhoodRelations();
    updateContacts();
    updateParticleInteraction();
}

/**
 * @brief Checks if monomers are making or breaking contact.
 * 
 * This device funtion checks if, for any given monomer pair, the pair is making or breaking physical contact.
 * The function then adjusts the contact matrices accordingly.
 */
__device__ void updateNeighbourhoodRelations(
    const double3* pos,
    double3* matrix_con,
    double3* matrix_norm,
    double4* matrix_rot,
    double* matrix_comp,
    double* matrix_twist,
    const double* amon,
    const material* mat,
    const int* matIDs,
    const int Nmon
) {
    // Find thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    // Skip threads that dont correspond to a monomer pair.
    if (threadID < Nmon * Nmon) return;

    // Find the individual monomer indices.
    int i = threadID % Nmon;
    int j = threadID / Nmon;

    int matrix_i = i * Nmon + j;
    //int matrix_j = 1 * Nmon * Nmon + i * Nmon + j;

    if (i == j) return;

    // Calculate pair properties // TODO: Investigate if this can be done during simulation initialization instead of during each evaluation step...
    double3 pos_A = pos[i];
    double3 pos_B = pos[j];

    int mat_id_A = matIDs[i];
    int mat_id_B = matIDs[j];

    double amon_A = amon[i];
    double amon_B = amon[j];
    double R = (amon_A * amon_B) / (amon_A + amon_B);

    double nu_A = mat[mat_id_A].nu;
    double nu_B = mat[mat_id_B].nu;
    
    double E_A = mat[mat_id_A].E;
    double E_B = mat[mat_id_B].E;
    double Es = ((1 - nu_A * nu_A) / E_A) + ((1 - nu_B * nu_B) / E_B);
    Es = 1.0 / Es;

    double gamma_A = mat[mat_id_A].gamma;
    double gamma_B = mat[mat_id_B].gamma;
    double gamma = gamma_A + gamma_B - 2.0 / (1.0 / gamma_A + 1.0 / gamma_B);

    double a0 = pow(9 * PI * gamma * R * R / Es, 1.0 / 3.0);
    double delta_c = 0.5 * a0 * a0 / (R * pow(6.0, 1.0 / 3.0));
    
    double breaking_dist = amon_A + amon_B + delta_c;
    double contact_dist = amon_A + amon_B;

    // Calculate current pair state
    double distance = vec_dist_len(pos_A, pos_B);

    // TODO: Check if these branches lead to warp divergence!
    if (distance < contact_dist) {
        if (matrix_comp[matrix_i] == -1.0) {
            double3 n = vec_get_normal(pos_A, pos_B);

            // Initialize compression lenght
            matrix_comp[matrix_i] = amon_A + amon_B - distance;

            // Initialize contact pointers
            matrix_con[matrix_i].x = -n.x;
            matrix_con[matrix_i].y = -n.y;
            matrix_con[matrix_i].z = -n.z;

            // Initialize contact normal
            matrix_norm[matrix_i].x = n.x;
            matrix_norm[matrix_i].y = n.y;
            matrix_norm[matrix_i].z = n.z;

            // Initialize pointer rotation
            matrix_rot[matrix_i].w = 1;
            matrix_rot[matrix_i].x = 0;
            matrix_rot[matrix_i].y = 0;
            matrix_rot[matrix_i].z = 0;

            // Initialize twisting angle
            matrix_twist[matrix_i] = 0;
        }
    }

    if (distance > breaking_dist) {
        matrix_comp[matrix_i] = -1.0;

        matrix_con[matrix_i].x = 0;
        matrix_con[matrix_i].y = 0;
        matrix_con[matrix_i].z = 0;

        matrix_norm[matrix_i].x = 0;
        matrix_norm[matrix_i].y = 0;
        matrix_norm[matrix_i].z = 0;

        matrix_rot[matrix_i].w = 0;
        matrix_rot[matrix_i].x = 0;
        matrix_rot[matrix_i].y = 0;
        matrix_rot[matrix_i].z = 0;

        matrix_twist[matrix_i] = 0;
    }
} 

/**
 * @brief Calculates the current contact pointer.
 * 
 * This device function calculates the current contact pointer from the rotation matrix and initial pointer.
 */
__device__ void updateNormal(
    const double3 n_A, // This is getting rotated
    double3* matrix_con, // This is where the rotated result is getting stored
    const double4* matrix_rot,
    const int Nmon
) {
    // Find thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    // Skip threads that dont correspond to a monomer pair.
    if (threadID < Nmon * Nmon) return; // TODO: Make sure this is necessary, will cause warp divergence

    // Find the individual monomer indices.
    int i = threadID % Nmon;
    int j = threadID / Nmon;

    // Find the index of the monomer pair in the contact matrix
    int matrix_i = i * Nmon + j;
    
    double4 rot_A = matrix_rot[matrix_i]; // TODO: Rename

    double3 initN_A;
    initN_A.x = 2.0 * ((0.5 - rot_A.y * rot_A.y - rot_A.z * rot_A.z) * n_A.x + (rot_A.x * rot_A.y + rot_A.z * rot_A.w) * n_A.y + (rot_A.x * rot_A.z - rot_A.y * rot_A.w) * n_A.z);
    initN_A.y = 2.0 * ((rot_A.x * rot_A.y - rot_A.z * rot_A.w) * n_A.x + (0.5 - rot_A.x * rot_A.x - rot_A.z * rot_A.z) * n_A.y + (rot_A.y * rot_A.z + rot_A.x * rot_A.w) * n_A.z);
    initN_A.z = 2.0 * ((rot_A.x * rot_A.z + rot_A.y * rot_A.w) * n_A.x + (rot_A.y * rot_A.z - rot_A.x * rot_A.w) * n_A.y + (0.5 - rot_A.x * rot_A.x - rot_A.y * rot_A.y) * n_A.z);

    matrix_con[matrix_i].x = initN_A.x;
    matrix_con[matrix_i].y = initN_A.y;
    matrix_con[matrix_i].z = initN_A.z;
}

/**
 * @brief 
 * 
 * 
 */
__device__ void updateContacts(
    const double3* omega_curr,
    const double3* omega_tot_curr,
    const double3* torque_curr,
    double3* mag_curr,
    double4* matrix_rot_curr,
    const double3* matrix_comp_curr,
    const double* moment,
    const double timestep,
    const int Nmon
) {
    // FIXME: The indices 'threadID', 'i' and 'j' got confused here i think. Make sure they are used correctly!!
    // Find thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    // Skip threads that dont correspond to a monomer pair.
    if (threadID < Nmon * Nmon) return; // TODO: Make sure this is necessary, will cause warp divergence

    // Find the individual monomer indices.
    int i = threadID % Nmon;
    int j = threadID / Nmon;

    // Find the index of the monomer pair in the contact matrix
    int matrix_i = i * Nmon + j;

    double3 omega, omega_dot;
    double4 e_dot, e_ddot;
    double tmp = 0;

    double inv_moment = 1.0 / moment[i];

    omega_dot.x = inv_moment * torque_curr[i].x;

    // FIXME: This whole bit needs to be placed in a seperate kernel. At this point this code gets executed Nmon times instead of one because it is in a Nmon*Nmon Kernel!
    // I could/should (?) treat this the same way I treat the contact pointers with a quaternion that tracks the rotation of the magnetization that gets calculated each timestep from the initial magnetization.
    // COROTATE THE MAGNETIZATION
    omega.x = omega_tot_curr[i].x;
    omega.y = omega_tot_curr[i].y;
    omega.z = omega_tot_curr[i].z;

    double len_mag = vec_lenght(mag_curr[i]);
    double len_omega = vec_lenght(omega);

    if (len_mag * len_omega > 0) {
        double4 q_mag;

        double3 tmp_mag;
        tmp_mag.x = mag_curr[i].x;
        tmp_mag.y = mag_curr[i].y;
        tmp_mag.z = mag_curr[i].z;
        vec_normalize(tmp_mag);

        q_mag.w = 0;
        q_mag.x = tmp_mag.x;
        q_mag.y = tmp_mag.y;
        q_mag.z = tmp_mag.z;

        e_dot.w = -0.5 * (q_mag.x * omega.x + q_mag.y * omega.y + q_mag.z * omega.z);
        e_dot.x = 0.5 * (q_mag.w * omega.x - q_mag.y * omega.z + q_mag.z * omega.y);
        e_dot.y = 0.5 * (q_mag.w * omega.y - q_mag.z * omega.x + q_mag.x * omega.z);
        e_dot.z = 0.5 * (q_mag.w * omega.z - q_mag.x * omega.y + q_mag.y * omega.x);

        tmp = 0.5 * e_dot.w;

        e_ddot.w = -0.25 * (q_mag.w * vec_lenght_sq(omega) + 2.0 * (q_mag.x * omega_dot.x + q_mag.y * omega_dot.y + q_mag.z * omega_dot.z));
        e_ddot.x = tmp * omega.x + 0.5 * (q_mag.w * omega_dot.x - q_mag.y * omega_dot.z + q_mag.z * omega_dot.y);
        e_ddot.y = tmp * omega.y + 0.5 * (q_mag.w * omega_dot.y - q_mag.z * omega_dot.x + q_mag.x * omega_dot.z);
        e_ddot.z = tmp * omega.z + 0.5 * (q_mag.w * omega_dot.z - q_mag.x * omega_dot.y + q_mag.y * omega_dot.x);

        // Integration step
        q_mag.w += timestep * e_dot.w + 0.5 * timestep * timestep * e_ddot.w;
        q_mag.x += timestep * e_dot.x + 0.5 * timestep * timestep * e_ddot.x;
        q_mag.y += timestep * e_dot.y + 0.5 * timestep * timestep * e_ddot.y;
        q_mag.z += timestep * e_dot.z + 0.5 * timestep * timestep * e_ddot.z;
    
        tmp_mag.x = q_mag.x;
        tmp_mag.y = q_mag.y;
        tmp_mag.z = q_mag.z;

        double len_tmp = vec_lenght(tmp_mag);

        if (len_tmp > 0) {
            mag_curr[i].x = len_mag * tmp_mag.x / len_tmp;
            mag_curr[i].y = len_mag * tmp_mag.y / len_tmp;
            mag_curr[i].z = len_mag * tmp_mag.z / len_tmp;
        }
    }

    // INTEGRATE THE CONTACT POINTER QUATERNION
    double4 rot = matrix_rot_curr[matrix_i];

    e_dot.w = -0.5 * (rot.x * omega.x + rot.y * omega.y + rot.z * omega.z);
    e_dot.x = 0.5 * (rot.w * omega.x - rot.y * omega.z + rot.z * omega.y);
    e_dot.y = 0.5 * (rot.w * omega.y - rot.z * omega.x + rot.x * omega.z);
    e_dot.z = 0.5 * (rot.w * omega.z - rot.x * omega.y + rot.y * omega.x);

    tmp = 0.5 * e_dot.w;

    e_ddot.w = -0.25 * (rot.w * vec_lenght_sq(omega) + 2.0 * (rot.x * omega_dot.x + rot.y * omega_dot.y + rot.z * omega_dot.z));
    e_ddot.x = tmp * omega.x + 0.5 * (rot.w * omega_dot.x - rot.y * omega_dot.z + rot.z * omega_dot.y);
    e_ddot.y = tmp * omega.y + 0.5 * (rot.w * omega_dot.y - rot.z * omega_dot.x + rot.x * omega_dot.z);
    e_ddot.z = tmp * omega.z + 0.5 * (rot.w * omega_dot.z - rot.x * omega_dot.y + rot.y * omega_dot.x);

    // Integration step
    rot.w += timestep * e_dot.w + 0.5 * timestep * timestep * e_ddot.w;
    rot.x += timestep * e_dot.x + 0.5 * timestep * timestep * e_ddot.x;
    rot.y += timestep * e_dot.y + 0.5 * timestep * timestep * e_ddot.y;
    rot.z += timestep * e_dot.z + 0.5 * timestep * timestep * e_ddot.z;

    quat_normalize(rot);

    // Store the new rotation matrices.
    matrix_rot_curr[matrix_i].w = rot.w;
    matrix_rot_curr[matrix_i].x = rot.x;
    matrix_rot_curr[matrix_i].y = rot.y;
    matrix_rot_curr[matrix_i].z = rot.z;
}

/**
 * 
 */
__device__ double gpu_getJKRContactRadius(
    const double compression_length,
    const double r0,
    const double R
) {
    double c1_contact_radius = sqrt(r0);
    double c2_contact_radius = 0.25 * 2.0 / 3.0 * pow(r0, 1.5);

    // contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
    // use equilibrium contact radius as starting value
    double k = compression_length * R / 3.0;
    double x_pow3;
    double x_new;
    double x_old = c1_contact_radius;

    // use Newton-Raphson method to find root
    for (int i = 0; i < 20; ++i)
    {
        x_pow3 = x_old * x_old * x_old;
        x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius);

        if (std::abs(x_new - x_old) / x_new < 1.e-14)
            break;

        x_old = x_new;
    }

    return x_new * x_new;
}

/**
 * 
 */
__device__ void updateParticleInteraction(
    const double3* pos_next,
    double3* force_next,
    const double3* torque_curr,
    double3* torque_next,
    double3* dMdt_next,
    const double3* matrix_con_curr,
    double3* matrix_norm,
    const double3* omega_curr,
    const double3* omega_tot_curr,
    const double3* mag_curr,
    const double4* matrix_rot_curr,
    double* matrix_comp_curr,
    const double* matrix_twist_curr,
    const double* amon,
    const double* moment,
    const material* mat,
    const int* matIDs,
    const double3 B_ext,
    const double timestep,
    const int Nmon
) {
    // Find thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    // Skip threads that dont correspond to a monomer pair.
    if (threadID < Nmon * Nmon) return; // TODO: Make sure this is necessary, will cause warp divergence

    // Find the individual monomer indices.
    int i = threadID % Nmon;
    int j = threadID / Nmon;

    // Find the index of the monomer pair in the contact matrix
    int matrix_i = i * Nmon + j;
    int matrix_j = Nmon * Nmon + i * Nmon + j; // FIXME: I dont understand this formula...

    // Initialize some quantities
    vec_set(force_next[i], 0.0);
    vec_set(torque_next[i], 0.0);
    vec_set(dMdt_next[i], 0.0);

    double3 pos_A = pos_next[i];
    double3 pos_B = pos_next[j];

    int matID_A = matIDs[i];
    int matID_B = matIDs[j];

    double amon_A = amon[i];
    double amon_B = amon[j];

    // Calculate the current contact pointer and monomer distance.
    double3 n_c = vec_diff(pos_A, pos_B);
    double particle_distance = vec_lenght(n_c);
    vec_normalize(n_c);

    // TODO: What does this do?? It initializes force_new to a small value proportional to n_c. And why do you need force_temp here?
    vec3D force_tmp;
    force_tmp.x = 1e-12 * n_c.x;
    force_tmp.y = 1e-12 * n_c.y;
    force_tmp.z = 1e-12 * n_c.z;

    // FIXME: Investigate this!
    force_next[i].x += force_tmp.x;
    force_next[i].y += force_tmp.y;
    force_next[i].z += force_tmp.z;

    // FIXME: Magnetitic effect go here!

    if (matrix_comp_curr[matrix_i] == -1.) return;

    // FIXME: The inelastic motion branches are bound to cause tons of warp divergence...
    // Tracks if an update to the contact pointers is necessary (-> inelastic motion)
    bool contact_update_necessary = false;

    // Calculate required quantities for force calculations.
    double3 omega_A = omega_curr[i];
    double3 omega_B = omega_curr[j];

    double R = (amon_A * amon_B) / (amon_A + amon_B);

    double moment_A = moment[i];
    double moment_B = moment[j];

    double nu_A = mat[matID_A].nu;
    double nu_B = mat[matID_B].nu;

    double E_A = mat[matID_A].E;
    double E_B = mat[matID_B].E;
    double Es = (1 - nu_A * nu_A) / E_A + (1 - nu_B * nu_B) / E_B;
    Es = 1. / Es;

    double T_vis_A = mat[matID_A].tvis;
    double T_vis_B = mat[matID_B].tvis;
    double T_vis = 0.5 * (T_vis_A + T_vis_B);

    double G_A = 0.5 * E_A / (1. + nu_A);
    double G_B = 0.5 * E_B / (1. + nu_B);
    double Gs = (1. - nu_A * nu_A) / G_A + (1. - nu_B * nu_B) / G_B;
    Gs = 1. / Gs;

    double gamma_A = mat[matID_A].gamma;
    double gamma_B = mat[matID_B].gamma;
    double gamma = gamma_A + gamma_B - 2. / (1. / gamma_A + 1. / gamma_B);

    double a0 = pow(9 * PI * gamma * R * R / Es, 1. / 3.);
    double delta_c = 0.5 * a0 * a0 / (R * pow(6.0, 1.0 / 3.0));

    // Calculate current contact pointers.
    double3 n_A = matrix_con_curr[matrix_i];
    double3 n_B = matrix_con_curr[matrix_j];
    double3 delta_n = vec_diff(n_A, n_B);

    // ###############      FORCE CALCULATIONS FROM HERE      ###############
    double compression_length = amon_A + amon_B - particle_distance;
    double contact_radius = gpu_getJKRContactRadius(compression_length, a0, R);
            
    // Calculate the NORMAL FORCE.
    double Fc = 3 * PI * gamma * R;
    double force_elastic = 4.0 * Fc * (pow(contact_radius / a0, 3.0) - pow(contact_radius / a0, 3.0 / 2.0));

    force_next[i].x += force_elastic * n_c.x;
    force_next[i].y += force_elastic * n_c.y;
    force_next[i].z += force_elastic * n_c.z;

    // TODO: Research this.
    // Calculate the DAMPING FORCE.
    double old_compression_length = matrix_comp_curr[matrix_i];
    double vis_damp_const = 2.0 * T_vis / (nu_A * nu_B) * Es;
    double delta_dot = (compression_length - old_compression_length) / timestep;
    double force_damp = vis_damp_const * contact_radius * delta_dot;

    force_next[i].x += force_damp * n_c.x;
    force_next[i].y += force_damp * n_c.y;
    force_next[i].z += force_damp * n_c.z;

    // Determine sliding and rolling displacement.
    double dot = vec_dot(delta_n, n_c);

    // TODO: Is amon_A correct here?
    double3 displacement;
    displacement.x = amon_A * (delta_n.x - dot * n_c.x);
    displacement.y = amon_A * (delta_n.y - dot * n_c.y);
    displacement.z = amon_A * (delta_n.z - dot * n_c.z);

    double displacement_norm = vec_lenght(displacement);

    double crit_sliding_displacement_modifier = 1.0;
    double crit_sliding_displacement = crit_sliding_displacement_modifier * (2.0 - nu_A) / (16.0 * PI) * a0;

    // FIXME: This process needs to be redesigned, the contact pointers need to be updated symmetrically (in 2 threads at once), but each thread can only update the primary monomer it is assigned to.
    // When the inelastic sliding regime is regime is reached.
    if (displacement_norm > crit_sliding_displacement) {
        // Determine correction of contact pointers
        double3 displacement_correction;
        displacement_correction.x = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.x;
        displacement_correction.y = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.y;
        displacement_correction.z = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.z;

        // Calculate correction factor (see Wada et al. 2007 appendix for details)
        double inv_norm = 1.0 / vec_lenght_sq(displacement_correction);

        dot = vec_dot(n_A, displacement_correction);
        double alpha_A = 1.0 / (1.0 - dot * dot * inv_norm);

        dot = vec_dot(n_A, displacement_correction);
        double alpha_B = 1.0 / (1.0 - dot * dot * inv_norm);

        // Apply contact pointer corrections to the current contact pointers.
        double particle_radius_inv = 1.0 / amon_A;
        n_A.x -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.x;
        n_A.y -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.y;
        n_A.z -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.z;

        n_B.x += 0.5 * particle_radius_inv * alpha_B * displacement_correction.x;
        n_B.y += 0.5 * particle_radius_inv * alpha_B * displacement_correction.y;
        n_B.z += 0.5 * particle_radius_inv * alpha_B * displacement_correction.z;

        vec_normalize(n_A);
        vec_normalize(n_B);

        // Track that the contact pointers need to be updates.
        contact_update_necessary = true;
    }

    double crit_rolling_displacement = 0.5 * (mat[matID_A].xi + mat[matID_B].xi);

    displacement.x = R * (n_A.x + n_B.x);
    displacement.y = R * (n_A.y + n_B.y);
    displacement.z = R * (n_A.z + n_B.z);
    displacement_norm = vec_lenght(displacement);

    // When the inelastic rolling regime is reached.
    if (displacement_norm > crit_rolling_displacement) {
        // Determine correction of contact pointers
        double3 displacement_correction;
        displacement_correction.x = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.x;
        displacement_correction.y = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.y;
        displacement_correction.z = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.z;

        // Calculate correction factor (see Wada et al. _B007 appendix for details)
        double inv_norm = 1.0 / vec_lenght_sq(displacement_correction);

        dot = vec_dot(n_A, displacement_correction);
        double alpha_A = 1.0 / (1.0 - dot * dot * inv_norm);

        dot = vec_dot(n_B, displacement_correction);
        double alpha_B = 1.0 / (1.0 - dot * dot * inv_norm);

        // Apply the correction to the current contact pointers.
        double particle_radius_inv = 1.0 / amon_A;
        n_A.x -= particle_radius_inv * alpha_A * displacement_correction.x;
        n_A.y -= particle_radius_inv * alpha_A * displacement_correction.y;
        n_A.z -= particle_radius_inv * alpha_A * displacement_correction.z;

        n_B.x -= particle_radius_inv * alpha_B * displacement_correction.x;
        n_B.y -= particle_radius_inv * alpha_B * displacement_correction.y;
        n_B.z -= particle_radius_inv * alpha_B * displacement_correction.z;

        vec_normalize(n_A);
        vec_normalize(n_B);

        // Track that the contact pointers need to be updated.
        contact_update_necessary = true;
    }

    // Update the contact pointers when neccessary.
    if (contact_update_necessary)
    {
        //updateNormal(n_A, n_B, matrix_con, matrix_rot, i, j, Nmon); // FIXME: This needs to be redesigned!
    }

    // Calculate the SLIDING FORCE
    double sliding_modifier = 1.0; // TODO: Make this a macro or maybe just remove it.
    double k_s = sliding_modifier * 8.0 * Gs * a0;

    // Sliding displacement
    double tmp_A = amon_A * vec_dot(n_A, n_c);
    double tmp_B = amon_B * vec_dot(n_B, n_c);

    // FIXME: Check this in Wada07, but: This formula is incorrect it should read r_i * n_i - r_j * n_j + (r_i + r_j) * n_c
    double3 displacement_zeta;
    displacement_zeta.x = amon_A * n_A.x - amon_B * n_B.x - (tmp_A - tmp_B) * n_c.x;
    displacement_zeta.y = amon_A * n_A.y - amon_B * n_B.y - (tmp_A - tmp_B) * n_c.y;
    displacement_zeta.z = amon_A * n_A.z - amon_B * n_B.z - (tmp_A - tmp_B) * n_c.z;

    // FIXME: The following equations are using zeta and not zeta_0, see Wada07
    double3 tmp;
    tmp.x = amon_B * n_B.x - amon_A * n_A.x;
    tmp.y = amon_B * n_B.y - amon_A * n_A.y;
    tmp.z = amon_B * n_B.z - amon_A * n_A.z;

    // FIXME: I dont understand what this does, but it cant be right?
    double force_sliding = -k_s * vec_dot(displacement_zeta, tmp) / particle_distance;

    // FIXME: This is a hack, might resolve itself once the formulaes are correct.
    // Clamps the sliding force to a maximum value
    if (abs(force_sliding) > 1.0e-10)
        force_sliding = 1e-10;

    force_next[i].x += force_sliding * n_c.x;
    force_next[i].y += force_sliding * n_c.y;
    force_next[i].z += force_sliding * n_c.z;

    // Calculate the SLIDING TORQUE
    double3 torque_sliding;
    tmp = vec_cross(n_A, displacement_zeta); // FIXME: This uses zeta_0 and not zeta. See Wada07

    // TODO: Unnecessary, right?
    torque_sliding.x = - amon_A * k_s * tmp.x;
    torque_sliding.y = - amon_A * k_s * tmp.y;
    torque_sliding.z = - amon_A * k_s * tmp.z;

    torque_next[i].x += torque_sliding.x;
    torque_next[i].y += torque_sliding.y;
    torque_next[i].z += torque_sliding.z;

    // Calculate the ROLLING TORQUE
    double3 displacement_xi;
    double3 torque_rolling;

    double rolling_modifier = 1.0; // TODO: Change this into a macro or just remove it.
    double k_r = rolling_modifier * 4.0 * Fc / R; // Remove /R from here because it gets devided out later anyway.

    // The rolling displacement.
    displacement_xi.x = R * (n_A.x + n_B.x);
    displacement_xi.y = R * (n_A.y + n_B.y);
    displacement_xi.z = R * (n_A.z + n_B.z);

    tmp = vec_cross(n_A, displacement_xi);
    torque_rolling.x = -k_r * R * tmp.x;
    torque_rolling.y = -k_r * R * tmp.y;
    torque_rolling.z = -k_r * R * tmp.z;

    torque_next[i].x += torque_rolling.x;
    torque_next[i].y += torque_rolling.y;
    torque_next[i].z += torque_rolling.z;

    // Calculate the TWISTING FORCE
    double3 delta_omega_old, delta_omega_new, twisting_torque;
    double crit_twisting_displacement = 1.0 / (16.0 * PI); // TODO: Check this, should this not be PI / 16.0?

    double twisting_modifier = 1.0; // TODO: Change this into a macro or just remove it.
    double k_t = twisting_modifier * 16.0 / 3.0 * Gs * a0 * a0 * a0;
    
    // TODO: Check if this is true, when exactly does matrix_twist get modified?
    double twisting_displacement = matrix_twist_curr[matrix_i];
    double moment_inv_A = 1.0 / moment_A;
    double moment_inv_B = 1.0 / moment_B;

    // Store current contact normal.
    double3 n_c_old = matrix_norm[i * Nmon + j];

    // Difference in angular momenta, ie change in twisting displacement.
    delta_omega_old.x = omega_A.x - omega_B.x;
    delta_omega_old.y = omega_A.y - omega_B.y;
    delta_omega_old.z = omega_A.z - omega_B.z;

    // update twisting displacement - use second order integration: omega^n+1 = omega^n + 0.5 (<delta_omega^n, n_c^n> + <delta_omega^n+1, n_c^n+1>)
    delta_omega_new.x = delta_omega_old.x + timestep * (moment_inv_A * torque_curr[i].x - moment_inv_B * torque_curr[j].x);
    delta_omega_new.y = delta_omega_old.y + timestep * (moment_inv_A * torque_curr[i].y - moment_inv_B * torque_curr[j].y);
    delta_omega_new.z = delta_omega_old.z + timestep * (moment_inv_A * torque_curr[i].z - moment_inv_B * torque_curr[j].z);

    twisting_displacement += 0.5 * timestep * (vec_dot(delta_omega_old, n_c_old) + vec_dot(delta_omega_new, n_c));

    // Clamps the twisting displacement // TODO: the clamping is only in the positive direction...
    if (twisting_displacement > crit_twisting_displacement)
        twisting_displacement = crit_twisting_displacement;

    twisting_torque.x = k_t * twisting_displacement * n_c.x;
    twisting_torque.y = k_t * twisting_displacement * n_c.y;
    twisting_torque.z = k_t * twisting_displacement * n_c.z;

    torque_next[i].x -= twisting_torque.x;
    torque_next[i].y -= twisting_torque.y;
    torque_next[i].z -= twisting_torque.z;/**/

    // Update the contact normal and compression lenghts.
    matrix_norm[matrix_i].x = n_c.x;
    matrix_norm[matrix_i].y = n_c.y;
    matrix_norm[matrix_i].z = n_c.z;

    matrix_comp_curr[matrix_i] = compression_length;
}