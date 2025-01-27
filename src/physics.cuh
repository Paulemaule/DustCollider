#pragma once

/**
 * @brief Depth-First search for single cluster recognition.
 * 
 * Searches the entire graph starting from node to find all connected nodes and assigns the current cluster ID to it.
 * 
 * @param node: The current node that is being looked at.
 * @param n: The number of monomers.
 * @param *matrix: The connection graph.
 * @param *cluster: The cluster memberships.
 * @param currentCluster: The ID of the current cluster.
 */
inline void cpu_dfs(int node, int n, const double* matrix, int* cluster, int currentCluster)
{
    // Assign the current cluster ID to the node
    cluster[node] = currentCluster;

    // Explore all neighbors of the node
    for (int i = 0; i < n; ++i)
    {
        // Check if there's an edge and the node has not been assigned to a cluster
        if (matrix[node * n + i] != -1 && cluster[i] == -1)
            cpu_dfs(i, n, matrix, cluster, currentCluster);
    }
}

/**
 * @brief Determines clustermembership of all monomers.
 * 
 * @param Nmon: The number of monomers.
 * @param *matrix: The connection graph.
 * @param *cluster: The cluster memberships.
 */
inline void cpu_findConnectedComponents(int Nmon, const double* matrix, int* cluster)
{
    int currentCluster = 0;  // Cluster ID counter

    // Iterate through all nodes
    for (int i = 0; i < Nmon; i++)
    {
        if (cluster[i] == -1)
        {  // If node hasn't been assigned to a cluster
            cpu_dfs(i, Nmon, matrix, cluster, currentCluster);  // Perform DFS
            currentCluster++;  // Move to the next cluster ID
        }
    }
}

/**
 * @brief Updates the contact matrices.
 * 
 * This functions checks if monomers are close enough to make contact or far enough to break contact and initializes/zeros the contact matrices accordingly.
 * 
 * @param *pos: Array of positions.
 * @param *matrix_con: Array of contact pointers.
 * @param *matrix_norm: Array of contact normals.
 * @param *matrix_rot: Array of contact rotations.
 * @param *matrix_comp: Array of contact compression lenghts.
 * @param *matrix_twist: Array of contact twisting angles.
 * @param *amon: Array of monomer radii.
 * @param *matIDs: Array of material IDs.
 * @param *Nmon: Number of monomers.
 */
inline void cpu_updateNeighbourhoodRelations(vec3D* pos, vec3D* matrix_con, vec3D* matrix_norm, quat* matrix_rot, double* matrix_comp, double* matrix_twist, double* amon, material* mat, int* matIDs, int Nmon)
{
    // Iterate over monomer pairs
    for (int i = 0; i < Nmon; i++)
    {
        for (int j = 1; j < Nmon; j++)
        {
            if (i == j)
                continue;

            // Determine and calculate necessary values.
            vec3D pos_A = pos[i];
            vec3D pos_B = pos[j];

            int mat_id_A = matIDs[i];
            int mat_id_B = matIDs[j];

            double a_mon_A = amon[i];
            double a_mon_B = amon[j];
            double R = (a_mon_A * a_mon_B) / (a_mon_A + a_mon_B); // Reduced radius

            double nu_A = mat[mat_id_A].nu;
            double nu_B = mat[mat_id_B].nu;

            double E_A = mat[mat_id_A].E;
            double E_B = mat[mat_id_B].E;

            double Es = (1 - nu_A * nu_A) / E_A + (1 - nu_B * nu_B) / E_B; // "Reduced" Young's modulus
            Es = 1. / Es;

            double gamma_A = mat[mat_id_A].gamma;
            double gamma_B = mat[mat_id_B].gamma;
            double gamma = gamma_A + gamma_B - 2. / (1. / gamma_A + 1. / gamma_B); // Contact surface energy

            double a0 = pow(9 * PI * gamma * R * R / Es, 1. / 3.);
            double delta_c = 0.5 * a0 * a0 / (R * pow(6.0, 1.0 / 3.0));
            double breaking_dist = a_mon_A + a_mon_B + delta_c;

            double contact_distance = (a_mon_A + a_mon_B);
            double distance = cpu_vec3D_distance(pos_A, pos_B);

            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

            // Check if monomer surfaces are in contact.
            if (distance < contact_distance)
            {
                // Check if monomers are allready marked as in contact
                if (matrix_comp[i * Nmon + j] == -1.)
                {
                    vec3D n = cpu_vec3D_get_normal(pos_A, pos_B);

                    // Initialize contact pointers
                    matrix_con[index_A].x = -n.x;
                    matrix_con[index_A].y = -n.y;
                    matrix_con[index_A].z = -n.z;
                    matrix_con[index_B].x = n.x; // FIXME: Potential memory race!
                    matrix_con[index_B].y = n.y;
                    matrix_con[index_B].z = n.z;

                    // Initialize contact pointer rotation
                    matrix_rot[index_A].e0 = 1;
                    matrix_rot[index_A].e1 = 0;
                    matrix_rot[index_A].e2 = 0;
                    matrix_rot[index_A].e3 = 0;
                    matrix_rot[index_B].e0 = 1; // FIXME: Potential memroy race!
                    matrix_rot[index_B].e1 = 0;
                    matrix_rot[index_B].e2 = 0;
                    matrix_rot[index_B].e3 = 0;

                    // Initialize contact normal
                    matrix_norm[i * Nmon + j].x = n.x;
                    matrix_norm[i * Nmon + j].x = n.y;
                    matrix_norm[i * Nmon + j].x = n.x;

                    double compression_length = a_mon_A + a_mon_B - distance;

                    matrix_comp[i * Nmon + j] = compression_length;
                    matrix_twist[i * Nmon + j] = 0;
                }
            }

            if (distance > breaking_dist)
            {
                // Empty the contact pointers.
                matrix_con[index_A].x = 0;
                matrix_con[index_A].y = 0;
                matrix_con[index_A].z = 0;
                matrix_con[index_B].x = 0;
                matrix_con[index_B].y = 0;
                matrix_con[index_B].z = 0;

                // TODO: What about matrix_rot ?

                // Mark the monomers as disconnected.
                matrix_comp[i * Nmon + j] = -1.;
                matrix_twist[i * Nmon + j] = 0;
            }
        }
    }
}

/**
 * @brief Updates the contact matrices.
 * 
 * This functions checks if monomers are close enough to make contact or far enough to break contact and initializes/zeros the contact matrices accordingly.
 * 
 * @param *pos: Array of positions.
 * @param *matrix_con: Array of contact pointers.
 * @param *matrix_rot: Array of contact pointer rotations.
 * @param *matrix_comp: Array of contact compression lenghts.
 * @param *matrix_twist: Array of contact twisting angles.
 * @param *amon: Array of monomer radii.
 * @param *mat: Array of material parameters.
 * @param *matIDs: Array of material IDs.
 * @param Nmon: Number of monomers.
 */
__global__ void gpu_updateNeighbourhoodRelations(vec3D* pos, vec3D* matrix_con, vec3D* matrix_norm, quat* matrix_rot, double* matrix_comp, double* matrix_twist, double* amon, material* mat, int* matIDs, int Nmon)
{
    // Calculate the index of the current thread, this index is associated with a single monomer.
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // Skip all threads that dont correnspond to a monomer.
    if (i < Nmon)
    {
        // Iterate over all other monomers.
        for (int j = 1; j < Nmon; j++)
        {
            if (i == j)
                continue;
            
            // Determine necessary values.
            vec3D pos_A = pos[i];
            vec3D pos_B = pos[j]; // FIXME: Memory race issues?

            int mat_id_A = matIDs[i];
            int mat_id_B = matIDs[j];

            double a_mon_A = amon[i];
            double a_mon_B = amon[j];
            double R = (a_mon_A * a_mon_B) / (a_mon_A + a_mon_B);

            double nu_A = mat[mat_id_A].nu;
            double nu_B = mat[mat_id_B].nu;

            double E_A = mat[mat_id_A].E;
            double E_B = mat[mat_id_B].E;

            double Es = (1 - nu_A * nu_A) / E_A + (1 - nu_B * nu_B) / E_B;
            Es = 1. / Es;

            double gamma_A = mat[mat_id_A].gamma;
            double gamma_B = mat[mat_id_B].gamma;
            double gamma = gamma_A + gamma_B - 2. / (1. / gamma_A + 1. / gamma_B);

            double a0 = pow(9 * PI * gamma * R * R / Es, 1. / 3.);
            double delta_c = 0.5 * a0 * a0 / (R * pow(6.0, 1.0 / 3.0));
            double breaking_dist = a_mon_A + a_mon_B + delta_c;

            double contact_distance = (a_mon_A + a_mon_B);
            double distance = gpu_vec3D_distance(pos_A, pos_B);

            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

            // Check if the monomers are close enough to be in contact.
            if (distance < contact_distance)
            {
                // Check if the monomers are allready marked as in contact.
                if (matrix_comp[i * Nmon + j] == -1.)
                {
                    vec3D n = gpu_vec3D_get_normal(pos_A, pos_B);

                    // Initialize contact pointers
                    matrix_con[index_A].x = -n.x;
                    matrix_con[index_A].y = -n.y;
                    matrix_con[index_A].z = -n.z;

                    // FIXME: Memory race issues!
                    matrix_con[index_B].x = n.x;
                    matrix_con[index_B].y = n.y;
                    matrix_con[index_B].z = n.z;

                    // Initialize contact pointer rotation.
                    matrix_rot[index_A].e0 = 1;
                    matrix_rot[index_A].e1 = 0;
                    matrix_rot[index_A].e2 = 0;
                    matrix_rot[index_A].e3 = 0;

                    // FIXME: Memory race issues!
                    matrix_rot[index_B].e0 = 1;
                    matrix_rot[index_B].e1 = 0;
                    matrix_rot[index_B].e2 = 0;
                    matrix_rot[index_B].e3 = 0;

                    // Initialize contact normal.
                    matrix_norm[i * Nmon + j].x = n.x;
                    matrix_norm[i * Nmon + j].x = n.y;
                    matrix_norm[i * Nmon + j].x = n.x;

                    // Initialize compression lenght.
                    double compression_length = a_mon_A + a_mon_B - distance; // TODO: Redundand line.
                    matrix_comp[i * Nmon + j] = compression_length;

                    // Initialize twisting angle.
                    matrix_twist[i * Nmon + j] = 0;
                }
            }

            // Check if monomers are far enough from each other to break contact.
            if (distance > breaking_dist)
            {
                // Zero the contact pointers.
                matrix_con[index_A].x = 0;
                matrix_con[index_A].y = 0;
                matrix_con[index_A].z = 0;

                matrix_con[index_B].x = 0;
                matrix_con[index_B].y = 0;
                matrix_con[index_B].z = 0;

                // TODO: What about contact normal and contact pointer rotation?

                // Mark the monomer pair as disconnected.
                matrix_comp[i * Nmon + j] = -1.;

                // Zero the twisting angle.
                matrix_twist[i * Nmon + j] = 0;
            }
        }
    }
}

/**
 * @brief Update the contact pointers (???).
 * 
 * @param &n_A: The current contact pointer of monomer A.
 * @param &n_B: The current contact pointer of monomer B.
 * @param *matrix_con: An array of contact pointers.
 * @param *matrix_rot: An array of contact pointer rotations.
 * @param i: Index of the first monomer.
 * @param j: Index of the second monomer.
 * @param Nmon: Number of monomers.
 */
inline void cpu_updateNormal(vec3D& n_A, vec3D& n_B, vec3D* matrix_con, quat* matrix_rot, int i, int j, int Nmon)
{
    // Determine the indices of the contact pointers and rotations in the matrices.
    int index_A = 0 * Nmon * Nmon + i * Nmon + j;
    int index_B = 1 * Nmon * Nmon + i * Nmon + j;

    quat rot_A = matrix_rot[index_A];
    quat rot_B = matrix_rot[index_B];

    vec3D init_n_A, init_n_B;

    // Perform quaternion rotation on the contact pointers according to the contact pointer rotation matrix.
    init_n_A.x = 2.0 * ((0.5 - rot_A.e2 * rot_A.e2 - rot_A.e3 * rot_A.e3) * n_A.x + (rot_A.e1 * rot_A.e2 + rot_A.e3 * rot_A.e0) * n_A.y + (rot_A.e1 * rot_A.e3 - rot_A.e2 * rot_A.e0) * n_A.z);
    init_n_A.y = 2.0 * ((rot_A.e1 * rot_A.e2 - rot_A.e3 * rot_A.e0) * n_A.x + (0.5 - rot_A.e1 * rot_A.e1 - rot_A.e3 * rot_A.e3) * n_A.y + (rot_A.e2 * rot_A.e3 + rot_A.e1 * rot_A.e0) * n_A.z);
    init_n_A.z = 2.0 * ((rot_A.e1 * rot_A.e3 + rot_A.e2 * rot_A.e0) * n_A.x + (rot_A.e2 * rot_A.e3 - rot_A.e1 * rot_A.e0) * n_A.y + (0.5 - rot_A.e1 * rot_A.e1 - rot_A.e2 * rot_A.e2) * n_A.z);

    init_n_B.x = 2.0 * ((0.5 - rot_B.e2 * rot_B.e2 - rot_B.e3 * rot_B.e3) * n_B.x + (rot_B.e1 * rot_B.e2 + rot_B.e3 * rot_B.e0) * n_B.y + (rot_B.e1 * rot_B.e3 - rot_B.e2 * rot_B.e0) * n_B.z);
    init_n_B.y = 2.0 * ((rot_B.e1 * rot_B.e2 - rot_B.e3 * rot_B.e0) * n_B.x + (0.5 - rot_B.e1 * rot_B.e1 - rot_B.e3 * rot_B.e3) * n_B.y + (rot_B.e2 * rot_B.e3 + rot_B.e1 * rot_B.e0) * n_B.z);
    init_n_B.z = 2.0 * ((rot_B.e1 * rot_B.e3 + rot_B.e2 * rot_B.e0) * n_B.x + (rot_B.e2 * rot_B.e3 - rot_B.e1 * rot_B.e0) * n_B.y + (0.5 - rot_B.e1 * rot_B.e1 - rot_B.e2 * rot_B.e2) * n_B.z);

    // Store the rotated contact pointers.
    matrix_con[index_A].x = init_n_A.x;
    matrix_con[index_A].y = init_n_A.y;
    matrix_con[index_A].z = init_n_A.z;

    matrix_con[index_B].x = init_n_B.x;
    matrix_con[index_B].y = init_n_B.y;
    matrix_con[index_B].z = init_n_B.z;
}

/**
 * @brief Update the contact pointers (???).
 * 
 * @param &n_A: The current contact pointer of monomer A.
 * @param &n_B: The current contact pointer of monomer B.
 * @param *matrix_con: An array of contact pointers.
 * @param *matrix_rot: An array of contact pointer rotations.
 * @param i: Index of the first monomer.
 * @param j: Index of the second monomer.
 * @param Nmon: Number of monomers.
 */
__device__ void gpu_updateNormal(vec3D& n_A, vec3D& n_B, vec3D* matrix_con, quat* matrix_rot, int i, int j, int Nmon)
{
    // Determine the indices of the contact pointers and rotations in the matrices.
    int index_A = 0 * Nmon * Nmon + i * Nmon + j;
    int index_B = 1 * Nmon * Nmon + i * Nmon + j;

    quat rot_A = matrix_rot[index_A];
    quat rot_B = matrix_rot[index_B];

    // Perform quaternion based rotation on the contact pointers.
    vec3D init_n_A, init_n_B;

    init_n_A.x = 2.0 * ((0.5 - rot_A.e2 * rot_A.e2 - rot_A.e3 * rot_A.e3) * n_A.x + (rot_A.e1 * rot_A.e2 + rot_A.e3 * rot_A.e0) * n_A.y + (rot_A.e1 * rot_A.e3 - rot_A.e2 * rot_A.e0) * n_A.z);
    init_n_A.y = 2.0 * ((rot_A.e1 * rot_A.e2 - rot_A.e3 * rot_A.e0) * n_A.x + (0.5 - rot_A.e1 * rot_A.e1 - rot_A.e3 * rot_A.e3) * n_A.y + (rot_A.e2 * rot_A.e3 + rot_A.e1 * rot_A.e0) * n_A.z);
    init_n_A.z = 2.0 * ((rot_A.e1 * rot_A.e3 + rot_A.e2 * rot_A.e0) * n_A.x + (rot_A.e2 * rot_A.e3 - rot_A.e1 * rot_A.e0) * n_A.y + (0.5 - rot_A.e1 * rot_A.e1 - rot_A.e2 * rot_A.e2) * n_A.z);

    init_n_B.x = 2.0 * ((0.5 - rot_B.e2 * rot_B.e2 - rot_B.e3 * rot_B.e3) * n_B.x + (rot_B.e1 * rot_B.e2 + rot_B.e3 * rot_B.e0) * n_B.y + (rot_B.e1 * rot_B.e3 - rot_B.e2 * rot_B.e0) * n_B.z);
    init_n_B.y = 2.0 * ((rot_B.e1 * rot_B.e2 - rot_B.e3 * rot_B.e0) * n_B.x + (0.5 - rot_B.e1 * rot_B.e1 - rot_B.e3 * rot_B.e3) * n_B.y + (rot_B.e2 * rot_B.e3 + rot_B.e1 * rot_B.e0) * n_B.z);
    init_n_B.z = 2.0 * ((rot_B.e1 * rot_B.e3 + rot_B.e2 * rot_B.e0) * n_B.x + (rot_B.e2 * rot_B.e3 - rot_B.e1 * rot_B.e0) * n_B.y + (0.5 - rot_B.e1 * rot_B.e1 - rot_B.e2 * rot_B.e2) * n_B.z);

    // Store the results.
    matrix_con[index_A].x = init_n_A.x;
    matrix_con[index_A].y = init_n_A.y;
    matrix_con[index_A].z = init_n_A.z;

    // FIXME: Potential memory race!
    matrix_con[index_B].x = init_n_B.x;
    matrix_con[index_B].y = init_n_B.y;
    matrix_con[index_B].z = init_n_B.z;
}

/**
 * @brief Updates the monomer positions based on current velocity and force.
 * 
 * @param *pos_old: The old positions array.
 * @param *pos_new: The new positions array.
 * @param *force_old: The old forces array.
 * @param *vel: The array of monomer velocities.
 * @param *double: An array of monomer masses.
 * @param time_step: The timestep of the simulation.
 * @param Nmon: The number of monomers.
 */
inline void cpu_predictor(vec3D* pos_old, vec3D* pos_new, vec3D* force_old, vec3D* vel, double* mass, double time_step, int Nmon)
{
    // Iterate over the monomers.
    for (int i = 0; i < Nmon; i++)
    {
        double mass_inv = 1. / mass[i];

        // TODO: Understand why force plays a role here...
        // Update the monomer positions.
        pos_new[i].x = pos_old[i].x + time_step * vel[i].x + 0.5 * mass_inv * time_step * time_step * force_old[i].x;
        pos_new[i].y = pos_old[i].y + time_step * vel[i].y + 0.5 * mass_inv * time_step * time_step * force_old[i].y;
        pos_new[i].z = pos_old[i].z + time_step * vel[i].z + 0.5 * mass_inv * time_step * time_step * force_old[i].z;
    }
}

/**
 * @brief Updates a monomers position based on its current velocity and force.
 * 
 * @param *pos_old: The old positions array.
 * @param *pos_new: The new positions array.
 * @param *force_old: The old forces array.
 * @param *vel: The array of monomer velocities.
 * @param *mass: An array of monomer masses.
 * @param time_step: The timestep of the simulation.
 * @param Nmon: The number of monomers.
 */
__global__ void gpu_predictor(vec3D* pos_old, vec3D* pos_new, vec3D* force_old, vec3D* vel, double* mass, double time_step, int Nmon)
{
    // Determines the index of the current thread, which is associated with a single monomer.
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // Skip all threads that dont correspond to a monomer.
    if (i < Nmon)
    {
        double mass_inv = 1. / mass[i];

        // Update the position of the monomer.
        pos_new[i].x = pos_old[i].x + time_step * vel[i].x + 0.5 * mass_inv * time_step * time_step * force_old[i].x;
        pos_new[i].y = pos_old[i].y + time_step * vel[i].y + 0.5 * mass_inv * time_step * time_step * force_old[i].y;
        pos_new[i].z = pos_old[i].z + time_step * vel[i].z + 0.5 * mass_inv * time_step * time_step * force_old[i].z;
    }
}

/**
 * @brief Updates amonomers velocity, angular velocity and magnetization.
 * 
 * @param *force_old: Array of old forces.
 * @param *force_new: Array of new forces.
 * @param *torque_old: Array of old torques.
 * @param *torque_new: Array of new torques.
 * @param *dMdt_old: Array of old change in magnetization.
 * @param *dMdt_new: Array of new change in magnetization.
 * @param *vel: Array of velocities.
 * @param *omega: Array of angular velocities.
 * @param *omega_tot: Array of total angular velocities due to rotation and curved trajectories.
 * @param *mag: Array of magnetizations.
 * @param *mass: Array of monomer masses.
 * @param *moment: Array of monomer moments of inertia.
 * @param *mat: Array of materials.
 * @param *int: Array of material IDs of materials.
 * @param time_step: The timestep of the simulation.
 * @param Nmon: The number of monomers.
 */
inline void cpu_corrector(vec3D* force_old, vec3D* force_new, vec3D* torque_old, vec3D* torque_new, vec3D* dMdt_old, vec3D* dMdt_new, vec3D* vel, vec3D* omega, vec3D* omega_tot, vec3D* mag, double* mass, double* moment, material* mat, int* matIDs, double time_step, int Nmon)
{
    // Iterate over the monomers.
    for (int i = 0; i < Nmon; i++)
    {
        double mass_inv = 1. / mass[i];
        double moment_of_inertia_inv = 1. / moment[i];

        // Calculate the acceleration of the monomer
        vec3D acc;

        acc.x = 0.5 * mass_inv * (force_new[i].x + force_old[i].x);
        acc.y = 0.5 * mass_inv * (force_new[i].y + force_old[i].y);
        acc.z = 0.5 * mass_inv * (force_new[i].z + force_old[i].z);

        // Update the velocities.
        vel[i].x += acc.x * time_step;
        vel[i].y += acc.y * time_step;
        vel[i].z += acc.z * time_step;

        // Update the angular momenta.
        omega[i].x += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].x + torque_old[i].x);
        omega[i].y += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].y + torque_old[i].y);
        omega[i].z += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].z + torque_old[i].z);

        // Determine the curvature of the trajectory.
        double v_sq = cpu_vec3D_length_sq(vel[i]);
        vec3D cross = cpu_vec3D_cross(vel[i], acc);

        // Update the total angular momentum.
        if (v_sq > 0)
        {
            omega_tot[i].x = omega[i].x + cross.x / v_sq;
            omega_tot[i].y = omega[i].y + cross.y / v_sq;
            omega_tot[i].z = omega[i].z + cross.z / v_sq;
        }
        else
        {
            omega_tot[i].x = omega[i].x;
            omega_tot[i].y = omega[i].y;
            omega_tot[i].z = omega[i].z;
        }

        // Update the magnetization.
        int mat_id = matIDs[i];
        double chi = mat[mat_id].chi;
        double Msat = mat[mat_id].Msat;

        if (abs(chi) > 0) //magnetic material
        {
            // Update the magnetization
            mag[i].x += 0.5 * time_step * (dMdt_new[i].x + dMdt_old[i].x);
            mag[i].y += 0.5 * time_step * (dMdt_new[i].y + dMdt_old[i].y);
            mag[i].z += 0.5 * time_step * (dMdt_new[i].z + dMdt_old[i].z);

            // Ensure the magnetization abides by saturation magnetization.
            double len_mag = cpu_vec3D_length(mag[i]);

            if (len_mag > 0)
            {
                if (abs(chi) > LIMIT_FER) // Ferromagnetic materials.
                {
                    // Rescale magnetization to saturation magnetization.
                    mag[i].x = Msat * mag[i].x / len_mag;
                    mag[i].y = Msat * mag[i].y / len_mag;
                    mag[i].z = Msat * mag[i].z / len_mag;
                }
                else // Non ferromagnetic materials.
                {
                    if (len_mag > Msat)
                    {
                        // Rescale magnetization to saturation magnetization if it exceeds it.
                        mag[i].x = Msat * mag[i].x / len_mag;
                        mag[i].y = Msat * mag[i].y / len_mag;
                        mag[i].z = Msat * mag[i].z / len_mag;
                    }
                }
            }
        }
    }
}

/**
 * @brief Updates velocity, angular velocity and magnetization of a single monomer.
 * 
 * @param *force_old: Array of old forces.
 * @param *force_new: Array of new forces.
 * @param *torque_old: Array of old torques.
 * @param *torque_new: Array of new torques.
 * @param *dMdt_old: Array of old change in magnetization.
 * @param *dMdt_new: Array of new change in magnetization.
 * @param *vel: Array of velocities.
 * @param *omega: Array of angular velocities.
 * @param *omega_tot: Array of total angular velocities due to rotation and curved trajectories.
 * @param *mag: Array of magnetizations.
 * @param *mass: Array of monomer masses.
 * @param *moment: Array of monomer moments of inertia.
 * @param *mat: Array of materials.
 * @param *int: Array of material IDs of materials.
 * @param time_step: The timestep of the simulation.
 * @param Nmon: The number of monomers.
 */
__global__ void gpu_corrector(vec3D* force_old, vec3D* force_new, vec3D* torque_old, vec3D* torque_new, vec3D* dMdt_old, vec3D* dMdt_new, vec3D* vel, vec3D* omega, vec3D* omega_tot, vec3D* mag, double* mass, double* moment, material* mat, int* matIDs, double time_step, int Nmon)
{
    // Determine thread index which corrensponds to a single monomer.
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // Skips indices which do not correnspond to a monomer.
    if (i < Nmon)
    {
        double mass_inv = 1. / mass[i];
        double moment_of_inertia_inv = 1. / moment[i];

        // Calculate the acceleration.
        vec3D acc;

        acc.x = 0.5 * mass_inv * (force_new[i].x + force_old[i].x);
        acc.y = 0.5 * mass_inv * (force_new[i].y + force_old[i].y);
        acc.z = 0.5 * mass_inv * (force_new[i].z + force_old[i].z);

        // Update the velocity.
        vel[i].x += acc.x * time_step;
        vel[i].y += acc.y * time_step;
        vel[i].z += acc.z * time_step;

        // Update the angular momentum.
        omega[i].x += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].x + torque_old[i].x);
        omega[i].y += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].y + torque_old[i].y);
        omega[i].z += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].z + torque_old[i].z);

        // Calculate the trajectories curvature.
        double v_sq = gpu_vec3D_length_sq(vel[i]);
        vec3D cross = gpu_vec3D_cross(vel[i], acc);

        // Calculate the total angular momentum.
        if (v_sq > 0)
        {
            omega_tot[i].x = omega[i].x + cross.x / v_sq;
            omega_tot[i].y = omega[i].y + cross.y / v_sq;
            omega_tot[i].z = omega[i].z + cross.z / v_sq;
        }
        else
        {
            omega_tot[i].x = omega[i].x;
            omega_tot[i].y = omega[i].y;
            omega_tot[i].z = omega[i].z;
        }

        // Update the magnetization.
        int mat_id = matIDs[i];
        double chi = mat[mat_id].chi;
        double Msat = mat[mat_id].Msat;

        if (abs(chi) > 0)
        {
            // Update the magnetization.
            mag[i].x += 0.5 * time_step * (dMdt_new[i].x + dMdt_old[i].x);
            mag[i].y += 0.5 * time_step * (dMdt_new[i].y + dMdt_old[i].y);
            mag[i].z += 0.5 * time_step * (dMdt_new[i].z + dMdt_old[i].z);

            double len_mag = gpu_vec3D_length(mag[i]);

            if (len_mag > 0)
            {
                if (abs(chi) > LIMIT_FER) // Ferromagnetic materials.
                {
                    // Rescale the magnetization to the saturation magnetization.
                    mag[i].x = Msat * mag[i].x / len_mag;
                    mag[i].y = Msat * mag[i].y / len_mag;
                    mag[i].z = Msat * mag[i].z / len_mag;
                }
                else // Non ferromagnetic materials.
                {
                    // Rescale the magnetization to saturation magnetization if it is exceeded.
                    if (len_mag > Msat)
                    {
                        mag[i].x = Msat * mag[i].x / len_mag;
                        mag[i].y = Msat * mag[i].y / len_mag;
                        mag[i].z = Msat * mag[i].z / len_mag;
                    }
                }
            }
        }
    }
}

/**
 * @brief Exchanges the old pointers for position, force, torque and magnetization change with the new.
 * 
 * @param *&pos_old: The old positions array.
 * @param *&pos_new: The new positions array.
 * @param *&force_old: The old forces array.
 * @param *&force_new: The new forces array.
 * @param *&torque_old: The old torques array.
 * @param *&torque_new: The new torques array.
 * @param *&dMdt_old: The old changes in magnetization array.
 * @param *&dMdt_new: The new changes in magnetization array.
 */
inline void switch_pointer(vec3D*& pos_old, vec3D*& pos_new, vec3D*& force_old, vec3D*& force_new, vec3D*& torque_old, vec3D*& torque_new, vec3D*& dMdt_old, vec3D*& dMdt_new)
{
    vec3D* temp;

    temp = pos_old;
    pos_old = pos_new;
    pos_new = temp;

    temp = force_old;
    force_old = force_new;
    force_new = temp;

    temp = torque_old;
    torque_old = torque_new;
    torque_new = temp;

    temp = dMdt_old;
    dMdt_old = dMdt_new;
    dMdt_new = temp;
}

/**
 * @brief Exchanges the old pointers for position, force, torque and magnetization change with the new.
 * 
 * @param *&pos_old: The old positions array.
 * @param *&pos_new: The new positions array.
 * @param *&force_old: The old forces array.
 * @param *&force_new: The new forces array.
 * @param *&torque_old: The old torques array.
 * @param *&torque_new: The new torques array.
 * @param *&dMdt_old: The old changes in magnetization array.
 * @param *&dMdt_new: The new changes in magnetization array.
 */
__device__ void gpu_switch_pointer(vec3D*& pos_old, vec3D*& pos_new, vec3D*& force_old, vec3D*& force_new, vec3D*& torque_old, vec3D*& torque_new, vec3D*& dMdt_old, vec3D*& dMdt_new)
{
    vec3D* temp;

    temp = pos_old;
    pos_old = pos_new;
    pos_new = temp;

    temp = force_old;
    force_old = force_new;
    force_new = temp;

    temp = torque_old;
    torque_old = torque_new;
    torque_new = temp;

    temp = dMdt_old;
    dMdt_old = dMdt_new;
    dMdt_new = temp;
}

/**
 * @brief Corotates the magnetization and calculates the new rotation matrix for the contact pointers.
 * 
 * @param *omega: Array of the angular momenta.
 * @param *omega_tot: Array of the total angular momenta (self rotation + curved trajectory).
 * @param *torque: Array of the torques.
 * @param *mag: Array of the magnetizations.
 * @param *matrix_rot: Array of the contact pointer rotations.
 * @param *matrix_comp: Array of the contact compression lenghts.
 * @param *moment: Array of the moments of inertia.
 * @param Nmon: Number of monomers.
 * @param time_step: The timestep of the simulation.
 */
inline void cpu_updateContacts(vec3D* omega, vec3D* omega_tot, vec3D* torque, vec3D* mag, quat* matrix_rot, double* matrix_comp, double* moment, int Nmon, double time_step)
{
    // Iterate over the monomers.
    for (int i = 0; i < Nmon; i++)
    {
        vec3D omega_A;
        vec3D omega_A_dot;
        double moment_A_inv = 1.0 / moment[i];

        quat e_dot;
        quat e_ddot;
        double temp = 0;

        // Calculate change in angular momentum based on total torque
        omega_A_dot.x = moment_A_inv * torque[i].x;
        omega_A_dot.y = moment_A_inv * torque[i].y;
        omega_A_dot.z = moment_A_inv * torque[i].z;

        // FIXME: This is unphysical, the rotation of the magnetization is tied to self rotation NOT trajectory curvature...
        omega_A.x = omega_tot[i].x;
        omega_A.y = omega_tot[i].y;
        omega_A.z = omega_tot[i].z;

        double len_mag = cpu_vec3D_length(mag[i]);
        double len_omega_A = cpu_vec3D_length(omega_A);
        
        // When both magnetization and total anuglar momentum are not zero.
        if (len_mag * len_omega_A > 0)
        {
            quat q_mag;

            // Calculate a normalized magnetization.
            vec3D tmp_mag;

            tmp_mag.x = mag[i].x;
            tmp_mag.y = mag[i].y;
            tmp_mag.z = mag[i].z;

            cpu_vec3D_normalize(tmp_mag);

            // Corotate the magnetization based on the total angular momentum (// TODO: ?).
            q_mag.e0 = 0;
            q_mag.e1 = tmp_mag.x;
            q_mag.e2 = tmp_mag.y;
            q_mag.e3 = tmp_mag.z;

            e_dot.e0 = -0.5 * (q_mag.e1 * omega_A.x + q_mag.e2 * omega_A.y + q_mag.e3 * omega_A.z);
            e_dot.e1 = 0.5 * (q_mag.e0 * omega_A.x - q_mag.e2 * omega_A.z + q_mag.e3 * omega_A.y);
            e_dot.e2 = 0.5 * (q_mag.e0 * omega_A.y - q_mag.e3 * omega_A.x + q_mag.e1 * omega_A.z);
            e_dot.e3 = 0.5 * (q_mag.e0 * omega_A.z - q_mag.e1 * omega_A.y + q_mag.e2 * omega_A.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (q_mag.e0 * cpu_vec3D_length_sq(omega_A) + 2.0 * (q_mag.e1 * omega_A_dot.x + q_mag.e2 * omega_A_dot.y + q_mag.e3 * omega_A_dot.z));
            e_ddot.e1 = temp * omega_A.x + 0.5 * (q_mag.e0 * omega_A_dot.x - q_mag.e2 * omega_A_dot.z + q_mag.e3 * omega_A_dot.y);
            e_ddot.e2 = temp * omega_A.y + 0.5 * (q_mag.e0 * omega_A_dot.y - q_mag.e3 * omega_A_dot.x + q_mag.e1 * omega_A_dot.z);
            e_ddot.e3 = temp * omega_A.z + 0.5 * (q_mag.e0 * omega_A_dot.z - q_mag.e1 * omega_A_dot.y + q_mag.e2 * omega_A_dot.x);

            q_mag.e0 += time_step * e_dot.e0 + 0.5 * time_step * time_step * e_ddot.e0;
            q_mag.e1 += time_step * e_dot.e1 + 0.5 * time_step * time_step * e_ddot.e1;
            q_mag.e2 += time_step * e_dot.e2 + 0.5 * time_step * time_step * e_ddot.e2;
            q_mag.e3 += time_step * e_dot.e3 + 0.5 * time_step * time_step * e_ddot.e3;

            tmp_mag.x = q_mag.e1;
            tmp_mag.y = q_mag.e2;
            tmp_mag.z = q_mag.e3;

            double len_tmp = cpu_vec3D_length(tmp_mag);

            // FIXME: Enable the corotation of magnetization.
            /*if (len_tmp > 0)
            {
                mag[i].x = len_mag * tmp_mag.x / len_tmp;
                mag[i].y = len_mag * tmp_mag.y / len_tmp;
                mag[i].z = len_mag * tmp_mag.z / len_tmp;
            }*/
        }

        // Iterate over all monomer pairs.
        for (int j = 1; j < Nmon; j++)
        {
            if (i == j)
                continue;

            // Skip unconnected particle pairs.
            if (matrix_comp[i * Nmon + j] == -1.)
                continue;

            vec3D omega_B;
            vec3D omega_B_dot;

            // Calculate the change in angular momentum due to torque for monomers A and B.
            omega_A_dot.x = moment_A_inv * torque[i].x;
            omega_A_dot.y = moment_A_inv * torque[i].y;
            omega_A_dot.z = moment_A_inv * torque[i].z;

            omega_A.x = omega[i].x;
            omega_A.y = omega[i].y;
            omega_A.z = omega[i].z;

            double moment_B_inv = 1.0 / moment[j];

            omega_B_dot.x = moment_B_inv * torque[j].x;
            omega_B_dot.y = moment_B_inv * torque[j].y;
            omega_B_dot.z = moment_B_inv * torque[j].z;

            omega_B.x = omega[j].x;
            omega_B.y = omega[j].y;
            omega_B.z = omega[j].z;

            // Determine the index of the monomers in the contact matrices.
            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

            quat rot_A = matrix_rot[index_A];
            quat rot_B = matrix_rot[index_B];

            // Determine contact pointer rotation quaternion for both monomers.
            e_dot.e0 = -0.5 * (rot_A.e1 * omega_A.x + rot_A.e2 * omega_A.y + rot_A.e3 * omega_A.z);
            e_dot.e1 = 0.5 * (rot_A.e0 * omega_A.x - rot_A.e2 * omega_A.z + rot_A.e3 * omega_A.y);
            e_dot.e2 = 0.5 * (rot_A.e0 * omega_A.y - rot_A.e3 * omega_A.x + rot_A.e1 * omega_A.z);
            e_dot.e3 = 0.5 * (rot_A.e0 * omega_A.z - rot_A.e1 * omega_A.y + rot_A.e2 * omega_A.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (rot_A.e0 * cpu_vec3D_length_sq(omega_A) + 2.0 * (rot_A.e1 * omega_A_dot.x + rot_A.e2 * omega_A_dot.y + rot_A.e3 * omega_A_dot.z));
            e_ddot.e1 = temp * omega_A.x + 0.5 * (rot_A.e0 * omega_A_dot.x - rot_A.e2 * omega_A_dot.z + rot_A.e3 * omega_A_dot.y);
            e_ddot.e2 = temp * omega_A.y + 0.5 * (rot_A.e0 * omega_A_dot.y - rot_A.e3 * omega_A_dot.x + rot_A.e1 * omega_A_dot.z);
            e_ddot.e3 = temp * omega_A.z + 0.5 * (rot_A.e0 * omega_A_dot.z - rot_A.e1 * omega_A_dot.y + rot_A.e2 * omega_A_dot.x);

            rot_A.e0 += time_step * e_dot.e0 + 0.5 * time_step * time_step * e_ddot.e0;
            rot_A.e1 += time_step * e_dot.e1 + 0.5 * time_step * time_step * e_ddot.e1;
            rot_A.e2 += time_step * e_dot.e2 + 0.5 * time_step * time_step * e_ddot.e2;
            rot_A.e3 += time_step * e_dot.e3 + 0.5 * time_step * time_step * e_ddot.e3;

            e_dot.e0 = -0.5 * (rot_B.e1 * omega_B.x + rot_B.e2 * omega_B.y + rot_B.e3 * omega_B.z);
            e_dot.e1 = 0.5 * (rot_B.e0 * omega_B.x - rot_B.e2 * omega_B.z + rot_B.e3 * omega_B.y);
            e_dot.e2 = 0.5 * (rot_B.e0 * omega_B.y - rot_B.e3 * omega_B.x + rot_B.e1 * omega_B.z);
            e_dot.e3 = 0.5 * (rot_B.e0 * omega_B.z - rot_B.e1 * omega_B.y + rot_B.e2 * omega_B.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (rot_B.e0 * cpu_vec3D_length_sq(omega_B) + 2.0 * (rot_B.e1 * omega_B_dot.x + rot_B.e2 * omega_B_dot.y + rot_B.e3 * omega_B_dot.z));
            e_ddot.e1 = temp * omega_B.x + 0.5 * (rot_B.e0 * omega_B_dot.x - rot_B.e2 * omega_B_dot.z + rot_B.e3 * omega_B_dot.y);
            e_ddot.e2 = temp * omega_B.y + 0.5 * (rot_B.e0 * omega_B_dot.y - rot_B.e3 * omega_B_dot.x + rot_B.e1 * omega_B_dot.z);
            e_ddot.e3 = temp * omega_B.z + 0.5 * (rot_B.e0 * omega_B_dot.z - rot_B.e1 * omega_B_dot.y + rot_B.e2 * omega_B_dot.x);

            rot_B.e0 += time_step * e_dot.e0 + 0.5 * time_step * time_step * e_ddot.e0;
            rot_B.e1 += time_step * e_dot.e1 + 0.5 * time_step * time_step * e_ddot.e1;
            rot_B.e2 += time_step * e_dot.e2 + 0.5 * time_step * time_step * e_ddot.e2;
            rot_B.e3 += time_step * e_dot.e3 + 0.5 * time_step * time_step * e_ddot.e3;

            cpu_quat_normalize(rot_A);
            cpu_quat_normalize(rot_B);

            // Store the new rotation matrices.
            matrix_rot[index_A].e0 = rot_A.e0;
            matrix_rot[index_A].e1 = rot_A.e1;
            matrix_rot[index_A].e2 = rot_A.e2;
            matrix_rot[index_A].e3 = rot_A.e3;

            matrix_rot[index_B].e0 = rot_B.e0;
            matrix_rot[index_B].e1 = rot_B.e1;
            matrix_rot[index_B].e2 = rot_B.e2;
            matrix_rot[index_B].e3 = rot_B.e3;
        }
    }
}

/**
 * @brief Corotates the magnetization and calculates the new rotation matrix for the contact pointers.
 * 
 * @param *omega: Array of the angular momenta.
 * @param *omega_tot: Array of the total angular momenta (self rotation + curved trajectory).
 * @param *torque: Array of the torques.
 * @param *mag: Array of the magnetizations.
 * @param *matrix_rot: Array of the contact pointer rotations.
 * @param *matrix_comp: Array of the contact compression lenghts.
 * @param *moment: Array of the moments of inertia.
 * @param Nmon: Number of monomers.
 * @param time_step: The timestep of the simulation.
 */
__global__ void gpu_updateContacts(vec3D* omega, vec3D* omega_tot, vec3D* torque, vec3D* mag, quat* matrix_rot, double* matrix_comp, double* moment, int Nmon, double time_step)
{
    // Determine the thread index which corresponds to a single monomer.
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // Skip threads that do not correspond to a monomer.
    if (i < Nmon)
    {
        vec3D omega_A;
        vec3D omega_A_dot;
        double moment_A_inv = 1.0 / moment[i];

        quat e_dot;
        quat e_ddot;
        double temp = 0;

        omega_A_dot.x = moment_A_inv * torque[i].x;
        omega_A_dot.y = moment_A_inv * torque[i].y;
        omega_A_dot.z = moment_A_inv * torque[i].z;

        // FIXME: The corotation of magnetization is based on omega not omega_tot!
        omega_A.x = omega_tot[i].x;
        omega_A.y = omega_tot[i].y;
        omega_A.z = omega_tot[i].z;

        double len_mag = gpu_vec3D_length(mag[i]);
        double len_omega_A = gpu_vec3D_length(omega_A);

        // When both magnetization and total angular momentum are not zero.
        if (len_mag * len_omega_A > 0)
        {
            quat q_mag;

            // Calculate the normalized magnetization
            vec3D tmp_mag;

            tmp_mag.x = mag[i].x;
            tmp_mag.y = mag[i].y;
            tmp_mag.z = mag[i].z;

            gpu_vec3D_normalize(tmp_mag);

            // Corotate the magnetization based on the total angular momentum.
            q_mag.e0 = 0;
            q_mag.e1 = tmp_mag.x;
            q_mag.e2 = tmp_mag.y;
            q_mag.e3 = tmp_mag.z;

            e_dot.e0 = -0.5 * (q_mag.e1 * omega_A.x + q_mag.e2 * omega_A.y + q_mag.e3 * omega_A.z);
            e_dot.e1 = 0.5 * (q_mag.e0 * omega_A.x - q_mag.e2 * omega_A.z + q_mag.e3 * omega_A.y);
            e_dot.e2 = 0.5 * (q_mag.e0 * omega_A.y - q_mag.e3 * omega_A.x + q_mag.e1 * omega_A.z);
            e_dot.e3 = 0.5 * (q_mag.e0 * omega_A.z - q_mag.e1 * omega_A.y + q_mag.e2 * omega_A.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (q_mag.e0 * gpu_vec3D_length_sq(omega_A) + 2.0 * (q_mag.e1 * omega_A_dot.x + q_mag.e2 * omega_A_dot.y + q_mag.e3 * omega_A_dot.z));
            e_ddot.e1 = temp * omega_A.x + 0.5 * (q_mag.e0 * omega_A_dot.x - q_mag.e2 * omega_A_dot.z + q_mag.e3 * omega_A_dot.y);
            e_ddot.e2 = temp * omega_A.y + 0.5 * (q_mag.e0 * omega_A_dot.y - q_mag.e3 * omega_A_dot.x + q_mag.e1 * omega_A_dot.z);
            e_ddot.e3 = temp * omega_A.z + 0.5 * (q_mag.e0 * omega_A_dot.z - q_mag.e1 * omega_A_dot.y + q_mag.e2 * omega_A_dot.x);

            q_mag.e0 += time_step * e_dot.e0 + 0.5 * time_step * time_step * e_ddot.e0;
            q_mag.e1 += time_step * e_dot.e1 + 0.5 * time_step * time_step * e_ddot.e1;
            q_mag.e2 += time_step * e_dot.e2 + 0.5 * time_step * time_step * e_ddot.e2;
            q_mag.e3 += time_step * e_dot.e3 + 0.5 * time_step * time_step * e_ddot.e3;

            tmp_mag.x = q_mag.e1;
            tmp_mag.y = q_mag.e2;
            tmp_mag.z = q_mag.e3;

            double len_tmp = gpu_vec3D_length(tmp_mag);

            // FIXME: Activate the corotation of magnetization.
            /*if (len_tmp > 0)
            {
                mag[i].x = len_mag * tmp_mag.x / len_tmp;
                mag[i].y = len_mag * tmp_mag.y / len_tmp;
                mag[i].z = len_mag * tmp_mag.z / len_tmp;
            }*/
        }

        // Iterate over all monomer pairs.
        for (int j = 1; j < Nmon; j++)
        {
            if (i == j)
                continue;

            // Skip monomers that are not in contact.
            if (matrix_comp[i * Nmon + j] == -1.)
                continue;

            vec3D omega_B;
            vec3D omega_B_dot;

            // Calculate the change in angular momentum due to torque.
            omega_A_dot.x = moment_A_inv * torque[i].x;
            omega_A_dot.y = moment_A_inv * torque[i].y;
            omega_A_dot.z = moment_A_inv * torque[i].z;

            omega_A.x = omega[i].x;
            omega_A.y = omega[i].y;
            omega_A.z = omega[i].z;

            double moment_B_inv = 1.0 / moment[j];

            omega_B_dot.x = moment_B_inv * torque[j].x;
            omega_B_dot.y = moment_B_inv * torque[j].y;
            omega_B_dot.z = moment_B_inv * torque[j].z;

            omega_B.x = omega[j].x;
            omega_B.y = omega[j].y;
            omega_B.z = omega[j].z;

            // Determine the indices of the monomers in the contact matrices.
            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

            quat rot_A = matrix_rot[index_A];
            quat rot_B = matrix_rot[index_B];

            // Determine contact pointer rotation quaternions.
            e_dot.e0 = -0.5 * (rot_A.e1 * omega_A.x + rot_A.e2 * omega_A.y + rot_A.e3 * omega_A.z);
            e_dot.e1 = 0.5 * (rot_A.e0 * omega_A.x - rot_A.e2 * omega_A.z + rot_A.e3 * omega_A.y);
            e_dot.e2 = 0.5 * (rot_A.e0 * omega_A.y - rot_A.e3 * omega_A.x + rot_A.e1 * omega_A.z);
            e_dot.e3 = 0.5 * (rot_A.e0 * omega_A.z - rot_A.e1 * omega_A.y + rot_A.e2 * omega_A.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (rot_A.e0 * gpu_vec3D_length_sq(omega_A) + 2.0 * (rot_A.e1 * omega_A_dot.x + rot_A.e2 * omega_A_dot.y + rot_A.e3 * omega_A_dot.z));
            e_ddot.e1 = temp * omega_A.x + 0.5 * (rot_A.e0 * omega_A_dot.x - rot_A.e2 * omega_A_dot.z + rot_A.e3 * omega_A_dot.y);
            e_ddot.e2 = temp * omega_A.y + 0.5 * (rot_A.e0 * omega_A_dot.y - rot_A.e3 * omega_A_dot.x + rot_A.e1 * omega_A_dot.z);
            e_ddot.e3 = temp * omega_A.z + 0.5 * (rot_A.e0 * omega_A_dot.z - rot_A.e1 * omega_A_dot.y + rot_A.e2 * omega_A_dot.x);

            rot_A.e0 += time_step * e_dot.e0 + 0.5 * time_step * time_step * e_ddot.e0;
            rot_A.e1 += time_step * e_dot.e1 + 0.5 * time_step * time_step * e_ddot.e1;
            rot_A.e2 += time_step * e_dot.e2 + 0.5 * time_step * time_step * e_ddot.e2;
            rot_A.e3 += time_step * e_dot.e3 + 0.5 * time_step * time_step * e_ddot.e3;

            e_dot.e0 = -0.5 * (rot_B.e1 * omega_B.x + rot_B.e2 * omega_B.y + rot_B.e3 * omega_B.z);
            e_dot.e1 = 0.5 * (rot_B.e0 * omega_B.x - rot_B.e2 * omega_B.z + rot_B.e3 * omega_B.y);
            e_dot.e2 = 0.5 * (rot_B.e0 * omega_B.y - rot_B.e3 * omega_B.x + rot_B.e1 * omega_B.z);
            e_dot.e3 = 0.5 * (rot_B.e0 * omega_B.z - rot_B.e1 * omega_B.y + rot_B.e2 * omega_B.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (rot_B.e0 * gpu_vec3D_length_sq(omega_B) + 2.0 * (rot_B.e1 * omega_B_dot.x + rot_B.e2 * omega_B_dot.y + rot_B.e3 * omega_B_dot.z));
            e_ddot.e1 = temp * omega_B.x + 0.5 * (rot_B.e0 * omega_B_dot.x - rot_B.e2 * omega_B_dot.z + rot_B.e3 * omega_B_dot.y);
            e_ddot.e2 = temp * omega_B.y + 0.5 * (rot_B.e0 * omega_B_dot.y - rot_B.e3 * omega_B_dot.x + rot_B.e1 * omega_B_dot.z);
            e_ddot.e3 = temp * omega_B.z + 0.5 * (rot_B.e0 * omega_B_dot.z - rot_B.e1 * omega_B_dot.y + rot_B.e2 * omega_B_dot.x);

            rot_B.e0 += time_step * e_dot.e0 + 0.5 * time_step * time_step * e_ddot.e0;
            rot_B.e1 += time_step * e_dot.e1 + 0.5 * time_step * time_step * e_ddot.e1;
            rot_B.e2 += time_step * e_dot.e2 + 0.5 * time_step * time_step * e_ddot.e2;
            rot_B.e3 += time_step * e_dot.e3 + 0.5 * time_step * time_step * e_ddot.e3;

            gpu_quat_normalize(rot_A);
            gpu_quat_normalize(rot_B);

            // Store the new rotation matrices.
            matrix_rot[index_A].e0 = rot_A.e0;
            matrix_rot[index_A].e1 = rot_A.e1;
            matrix_rot[index_A].e2 = rot_A.e2;
            matrix_rot[index_A].e3 = rot_A.e3;

            // FIXME: Potential memory race!
            matrix_rot[index_B].e0 = rot_B.e0;
            matrix_rot[index_B].e1 = rot_B.e1;
            matrix_rot[index_B].e2 = rot_B.e2;
            matrix_rot[index_B].e3 = rot_B.e3;
        }
    }
}

/**
 * @brief Calculates the JKR contact radius.
 * 
 * @param compression_lenght: The compression lenght of the 
 * @param r0: The equilibrium contact surface radius.
 * @param R: Reduced radius (?)
 */
inline double cpu_getJKRContactRadius(double compression_length, double r0, double R)
{
    double c1_contact_radius = sqrt(r0);
    double c2_contact_radius = 0.25 * 2.0 / 3.0 * pow(r0, 1.5);

    // contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
    // use equilibrium contact radius as starting value
    double k = compression_length * R / 3.0;
    double x_pow3;
    double x_new;
    double x_old = c1_contact_radius;

    // TODO: Check this algorithm.
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

// TODO: Check this algorithm.
/**
 * @brief Calculates the JKR contact radius.
 * 
 * @param compression_lenght: The compression lenght of the 
 * @param r0: (?)
 * @param R: Reduced radius (?)
 */
__device__ double gpu_getJKRContactRadius(double compression_length, double r0, double R)
{
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
 * @brief Calculates the forces and torques updates the contact normal and compression lenght.
 * 
 * @param *pos_new:
 * @param *force_new:
 * @param *torque_old:
 * @param *dMdt_new:
 * @param *matrix_con:
 * @param *matrix_norm:
 * @param *omega:
 * @param *omega_tot:
 * @param *mag:
 * @param *matrix_rot:
 * @param *matrix_comp:
 * @param *matrix_twist:
 * @param *amon:
 * @param *moment:
 * @param *mat:
 * @param *matIDs:
 * @param B_ext:
 * @param N_mon:
 * @param time_step:
 */
inline void cpu_updateParticleInteraction(vec3D* pos_new, vec3D* force_new, vec3D* torque_old, vec3D* torque_new, vec3D* dMdt_new, vec3D* matrix_con,
    vec3D* matrix_norm, vec3D* omega, vec3D* omega_tot, vec3D* mag, quat* matrix_rot, double* matrix_comp, double* matrix_twist, double* amon, double* moment, material* mat, int* matIDs, vec3D B_ext, int Nmon, double time_step)
{
    // Initialize forces, torques and chenge in magnetization to zero.
    memset(force_new, 0, Nmon * sizeof(vec3D));
    memset(torque_new, 0, Nmon * sizeof(vec3D));
    memset(dMdt_new, 0, Nmon * sizeof(vec3D));

    // Iterate over all particle pairs.
    for (int i = 0; i < Nmon; i++)
    {
        for (int j = 0; j < Nmon; j++)
        {
            if (i == j)
                continue;

            //current position
            vec3D pos_A = pos_new[i];
            vec3D pos_B = pos_new[j];

            //material IDs
            int mat_id_A = matIDs[i];
            int mat_id_B = matIDs[j];

            //monomer radius
            double a_mon_A = amon[i];
            double a_mon_B = amon[j];

            // Calculate current contact pointer and monomer distance.
            vec3D n_c = cpu_vec3D_diff(pos_A, pos_B);
            double particle_distance = cpu_vec3D_length(n_c);
            cpu_vec3D_normalize(n_c);

            // TODO: What does this do?? It initializes force_new to an extremely small value proportional to n_c. And why do you need force_temp here?
            vec3D force_tmp;

            force_tmp.x = 1e-12 * n_c.x;
            force_tmp.y = 1e-12 * n_c.y;
            force_tmp.z = 1e-12 * n_c.z;

            force_new[i].x += force_tmp.x;
            force_new[i].y += force_tmp.y;
            force_new[i].z += force_tmp.z;

            // FIXME: Activate magnetization calculation.
            // TODO: Document the magnetization code.
            // ###############    MAGNETIZATION CALCULATIONS FROM HERE    ###############
            /*double chi_A = mat[mat_id_A].chi;

            if (abs(chi_A) > 0)
            {
                double len_B = cpu_vec3D_length(B_ext);
                double Vmon_A = 4. / 3. * PI * a_mon_A * a_mon_A * a_mon_A;
                double Vmon_B = 4. / 3. * PI * a_mon_B * a_mon_B * a_mon_B;

                double Msat_A = mat[mat_id_A].Msat;

                double tau_ss_A = mat[mat_id_A].tss;
                double tau_sl_A = mat[mat_id_A].tsl;

                vec3D mag_A = mag[i];
                vec3D mag_B = mag[j];

                vec3D omega_tot_A = omega_tot[i];

                vec3D Mfin, Mpara, Mperp;

                cpu_vec3D_set(Mfin, 0);
                cpu_vec3D_set(Mpara, 0);
                cpu_vec3D_set(Mperp, 0);

                double len_mag_A = cpu_vec3D_length(mag_A);
                double len_mag_B = cpu_vec3D_length(mag_B);

                if (abs(chi_A) > LIMIT_FER) //ferromagnetic material
                {
                    double len_omega = cpu_vec3D_length(omega_tot_A);

                    if (len_mag_A > 0) //rescale current magnetization
                    {
                        mag_A.x = Msat_A * mag_A.x / len_mag_A;
                        mag_A.y = Msat_A * mag_A.y / len_mag_A;
                        mag_A.z = Msat_A * mag_A.z / len_mag_A;
                    }

                    if (len_B > 0) //future magnetization
                    {
                        Mfin.x = Msat_A * B_ext.x / len_B;
                        Mfin.y = Msat_A * B_ext.y / len_B;
                        Mfin.z = Msat_A * B_ext.z / len_B;
                    }
                    else
                    {
                        if (len_omega > 0)
                        {
                            Mfin.x = Msat_A * omega_tot_A.x / len_omega;
                            Mfin.y = Msat_A * omega_tot_A.y / len_omega;
                            Mfin.z = Msat_A * omega_tot_A.z / len_omega;
                        }
                        else //just put Mfin in x-direction for now if B_ext == 0 and omega_tot == 0
                        {
                            Mfin.x = Msat_A;
                            Mfin.y = 0;
                            Mfin.z = 0;
                        }
                    }

                    double len_Mfin = cpu_vec3D_length(Mfin);
                    double dot_magA_Mfin = cpu_vec3D_dot(mag_A, Mfin);

                    Mpara.x = dot_magA_Mfin * Mfin.x / (len_Mfin * len_Mfin);
                    Mpara.y = dot_magA_Mfin * Mfin.y / (len_Mfin * len_Mfin);
                    Mpara.z = dot_magA_Mfin * Mfin.z / (len_Mfin * len_Mfin);

                    Mperp.x = mag_A.x - Mpara.x;
                    Mperp.y = mag_A.y - Mpara.y;
                    Mperp.z = mag_A.z - Mpara.z;
                }
                else //paramegnetic material
                {
                    double chi_factor_A = chi_A / (chi_A + 1.);

                    //induced moment
                    Mfin.x = chi_factor_A * B_ext.x / mu0;
                    Mfin.y = chi_factor_A * B_ext.y / mu0;
                    Mfin.z = chi_factor_A * B_ext.z / mu0;

                    //Barnett moment
                    Mfin.x += chi_A * omega_tot_A.x / PROD_BARR;
                    Mfin.y += chi_A * omega_tot_A.y / PROD_BARR;
                    Mfin.z += chi_A * omega_tot_A.z / PROD_BARR;

                    double len_Mfin = cpu_vec3D_length(Mfin);

                    if (len_mag_A > Msat_A) //rescale current mag. if saturation is reached
                    {
                        if (len_mag_A > 0)
                        {
                            mag_A.x = Msat_A * mag_A.x / len_mag_A;
                            mag_A.y = Msat_A * mag_A.y / len_mag_A;
                            mag_A.z = Msat_A * mag_A.z / len_mag_A;
                        }
                    }

                    /*if (len_Mfin > Msat_A) //rescale future  mag. if saturation is reached
                    {
                        if (len_Mfin > 0)
                        {
                            Mfin.x = Msat_A * Mfin.x / len_Mfin;
                            Mfin.y = Msat_A * Mfin.y / len_Mfin;
                            Mfin.z = Msat_A * Mfin.z / len_Mfin;
                        }
                    }*/

                    /*double dot_magA_Mfin = cpu_vec3D_dot(mag_A, Mfin);

                    if (len_Mfin > 0)
                    {
                        Mpara.x = dot_magA_Mfin * Mfin.x / (len_Mfin * len_Mfin);
                        Mpara.y = dot_magA_Mfin * Mfin.y / (len_Mfin * len_Mfin);
                        Mpara.z = dot_magA_Mfin * Mfin.z / (len_Mfin * len_Mfin);

                        Mperp.x = mag_A.x - Mpara.x;
                        Mperp.y = mag_A.y - Mpara.y;
                        Mperp.z = mag_A.z - Mpara.z;
                    }
                    else
                    {
                        Mpara.x = 0;
                        Mpara.y = 0;
                        Mpara.z = 0;

                        Mperp.x = 0;
                        Mperp.y = 0;
                        Mperp.z = 0;
                    }
                }

                //solve Bloch equation
                dMdt_new[i].x = 0;
                dMdt_new[i].y = 0;
                dMdt_new[i].z = 0;

                if (len_B > 0)
                {
                    vec3D cross = cpu_vec3D_cross(mag_A, B_ext);

                    dMdt_new[i].x += GAMMA_E * cross.x;
                    dMdt_new[i].y += GAMMA_E * cross.y;
                    dMdt_new[i].z += GAMMA_E * cross.z;
                }

                if (tau_sl_A > 0)
                {
                    dMdt_new[i].x -= (mag_A.x - Mpara.x) / tau_sl_A;
                    dMdt_new[i].y -= (mag_A.y - Mpara.y) / tau_sl_A;
                    dMdt_new[i].z -= (mag_A.z - Mpara.z) / tau_sl_A;
                }

                if (tau_ss_A > 0)
                {
                    dMdt_new[i].x -= Mperp.x / tau_ss_A;
                    dMdt_new[i].y -= Mperp.y / tau_ss_A;
                    dMdt_new[i].z -= Mperp.z / tau_ss_A;
                }

                if (len_mag_A * len_mag_B > 0)
                {
                    vec3D mu_A, mu_B;

                    mu_A.x = mag_A.x * Vmon_A;
                    mu_A.y = mag_A.y * Vmon_A;
                    mu_A.z = mag_A.z * Vmon_A;

                    mu_B.x = mag_B.x * Vmon_B;
                    mu_B.y = mag_B.y * Vmon_B;
                    mu_B.z = mag_B.z * Vmon_B;

                    vec3D force_D, torque_D, torque_EdH;
                    vec3D torque_B = cpu_vec3D_cross(mu_A, B_ext);

                    vec3D vec_d = cpu_vec3D_diff(pos_B, pos_A); //todo: check for correct direction
                    double d = cpu_vec3D_length(vec_d);
                    double d2 = d * d;
                    double d3 = d2 * d;
                    double d4 = d2 * d2;
                    double d5 = d2 * d3;

                    //double fD = (3. * mu0) / (PIx4 * d5);
                    double fD = (3. * mu0) / (PIx4 * d4);
                    double tD = (mu0 / PIx4);
                    double tEdH = (Vmon_A / GAMMA_E);

                    double dot_muA_d = cpu_vec3D_dot(mu_A, vec_d);
                    double dot_muB_d = cpu_vec3D_dot(mu_B, vec_d);
                    double dot_muA_muB = cpu_vec3D_dot(mu_A, mu_B);

                    vec3D cross_d_muB = cpu_vec3D_cross(vec_d, mu_B);
                    vec3D cross_muA_muB = cpu_vec3D_cross(mu_A, mu_B);

                    //force_D.x = fD * (5. * vec_d.x * dot_muA_d * dot_muB_d / d2 - vec_d.x * dot_muA_muB - mu_A.x * dot_muB_d - mu_B.x * dot_muA_d);
                    //force_D.y = fD * (5. * vec_d.y * dot_muA_d * dot_muB_d / d2 - vec_d.y * dot_muA_muB - mu_A.y * dot_muB_d - mu_B.y * dot_muA_d);
                    //force_D.z = fD * (5. * vec_d.z * dot_muA_d * dot_muB_d / d2 - vec_d.z * dot_muA_muB - mu_A.z * dot_muB_d - mu_B.z * dot_muA_d);

                    force_D.x = fD * (-5. * vec_d.x * dot_muA_d * dot_muB_d + vec_d.x * dot_muA_muB + mu_A.x * dot_muB_d + mu_B.x * dot_muA_d);
                    force_D.y = fD * (-5. * vec_d.y * dot_muA_d * dot_muB_d + vec_d.y * dot_muA_muB + mu_A.y * dot_muB_d + mu_B.y * dot_muA_d);
                    force_D.z = fD * (-5. * vec_d.z * dot_muA_d * dot_muB_d + vec_d.z * dot_muA_muB + mu_A.z * dot_muB_d + mu_B.z * dot_muA_d);

                    torque_D.x = tD * 3. * dot_muB_d * cross_d_muB.x / d5 - cross_muA_muB.x / d3;
                    torque_D.y = tD * 3. * dot_muB_d * cross_d_muB.y / d5 - cross_muA_muB.y / d3;
                    torque_D.z = tD * 3. * dot_muB_d * cross_d_muB.z / d5 - cross_muA_muB.z / d3;

                    torque_EdH.x = -tEdH * dMdt_new[i].x;
                    torque_EdH.y = -tEdH * dMdt_new[i].y;
                    torque_EdH.z = -tEdH * dMdt_new[i].z;

                    force_new[i].x += force_D.x;
                    force_new[i].y += force_D.y;
                    force_new[i].z += force_D.z;

                    torque_new[i].x += torque_D.x;
                    torque_new[i].y += torque_D.y;
                    torque_new[i].z += torque_D.z;

                    torque_new[i].x += torque_EdH.x;
                    torque_new[i].y += torque_EdH.y;
                    torque_new[i].z += torque_EdH.z;

                    torque_new[i].x += torque_B.x;
                    torque_new[i].y += torque_B.y;
                    torque_new[i].z += torque_B.z;
                }
            }*/

            // Skip monomer pairs that are not in contact.
            if (matrix_comp[i * Nmon + j] == -1.)
                continue;

            // This tracks if any of the displacements have reached a critical value.
            bool update_contact_pointers = false;

            // Determine monomer indices in the contact matrices.
            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

            // Calculate required quantities for the force calculation.
            vec3D omega_A = omega[i];
            vec3D omega_B = omega[j];

            double R = (a_mon_A * a_mon_B) / (a_mon_A + a_mon_B);

            double moment_A = moment[i];
            double moment_B = moment[j];

            double nu_A = mat[mat_id_A].nu;
            double nu_B = mat[mat_id_B].nu;

            double E_A = mat[mat_id_A].E;
            double E_B = mat[mat_id_B].E;

            double T_vis_A = mat[mat_id_A].tvis;
            double T_vis_B = mat[mat_id_B].tvis;
            double T_vis = 0.5 * (T_vis_A + T_vis_B);

            double Es = (1 - nu_A * nu_A) / E_A + (1 - nu_B * nu_B) / E_B;
            Es = 1. / Es;

            double G_A = 0.5 * E_A / (1. + nu_A);
            double G_B = 0.5 * E_B / (1. + nu_B);

            double Gs = (1. - nu_A * nu_A) / G_A + (1. - nu_B * nu_B) / G_B;
            Gs = 1. / Gs;

            double gamma_A = mat[mat_id_A].gamma;
            double gamma_B = mat[mat_id_B].gamma;
            double gamma = gamma_A + gamma_B - 2. / (1. / gamma_A + 1. / gamma_B);

            double a0 = pow(9 * PI * gamma * R * R / Es, 1. / 3.);
            double delta_c = 0.5 * a0 * a0 / (R * pow(6.0, 1.0 / 3.0));

            // Calculate the current contact pointers.
            vec3D n_A = matrix_con[index_A];
            vec3D n_B = matrix_con[index_B];
            vec3D delta_n = cpu_vec3D_diff(n_A, n_B);

            // ###############      FORCE CALCULATIONS FROM HERE      ###############
            double compression_length = a_mon_A + a_mon_B - particle_distance;
            double contact_radius = cpu_getJKRContactRadius(compression_length, a0, R);

            // Calculate the COMPRESSION FORCE.
            double Fc = 3 * PI * gamma * R;
            double force_elastic = 4.0 * Fc * (pow(contact_radius / a0, 3.0) - pow(contact_radius / a0, 3.0 / 2.0));

            force_new[i].x += force_elastic * n_c.x;
            force_new[i].y += force_elastic * n_c.y;
            force_new[i].z += force_elastic * n_c.z;

            // TODO: Research this.
            // Calculate the DAMPING FORCE.
            double old_compression_length = matrix_comp[i * Nmon + j];
            double vis_damp_const = 2.0 * T_vis / (nu_A * nu_B) * Es;
            double delta_dot = (compression_length - old_compression_length) / time_step;
            double force_damp = vis_damp_const * contact_radius * delta_dot;

            force_new[i].x += force_damp * n_c.x;
            force_new[i].y += force_damp * n_c.y;
            force_new[i].z += force_damp * n_c.z;

            // Determine if sliding and rolling displacements are critical
            double dot = cpu_vec3D_dot(delta_n, n_c);

            // TODO: Is a_mon_A correct here?
            vec3D displacement;
            displacement.x = a_mon_A * (delta_n.x - dot * n_c.x);
            displacement.y = a_mon_A * (delta_n.y - dot * n_c.y);
            displacement.z = a_mon_A * (delta_n.z - dot * n_c.z);

            double displacement_norm = cpu_vec3D_length(displacement);

            double crit_sliding_displacement_modifier = 1.0;
            double crit_sliding_displacement = crit_sliding_displacement_modifier * (2.0 - nu_A) / (16.0 * PI) * a0;

            // When the inelastic sliding regime is reached
            if (displacement_norm > crit_sliding_displacement)
            {
                // Determine the correction to contact pointers.
                vec3D displacement_correction;
                displacement_correction.x = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.x;
                displacement_correction.y = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.y;
                displacement_correction.z = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.z;

                // Calculate the correction factor (see Wada et al. 2007 appendix).
                double inv_norm = 1.0 / cpu_vec3D_length_sq(displacement_correction);

                dot = cpu_vec3D_dot(n_A, displacement_correction);
                double alpha_A = 1.0 / (1.0 - dot * dot * inv_norm);

                dot = cpu_vec3D_dot(n_A, displacement_correction);
                double alpha_B = 1.0 / (1.0 - dot * dot * inv_norm);

                // Apply contact pointer corrections to the current contact pointers.
                double particle_radius_inv = 1.0 / a_mon_A;
                n_A.x -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.x;
                n_A.y -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.y;
                n_A.z -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.z;

                n_B.x += 0.5 * particle_radius_inv * alpha_B * displacement_correction.x;
                n_B.y += 0.5 * particle_radius_inv * alpha_B * displacement_correction.y;
                n_B.z += 0.5 * particle_radius_inv * alpha_B * displacement_correction.z;

                cpu_vec3D_normalize(n_A);
                cpu_vec3D_normalize(n_B);

                // Track that the contact pointers need to be updated.
                update_contact_pointers = true;
            }

            double crit_rolling_displacement = 0.5 * (mat[mat_id_A].xi + mat[mat_id_B].xi);

            displacement.x = R * (n_A.x + n_B.x);
            displacement.y = R * (n_A.y + n_B.y);
            displacement.z = R * (n_A.z + n_B.z);
            displacement_norm = cpu_vec3D_length(displacement);

            // When the inelastic rolling regime is reached.
            if (displacement_norm > crit_rolling_displacement)
            {
                // Determine the correction to the contact pointers.
                vec3D displacement_correction;
                displacement_correction.x = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.x;
                displacement_correction.y = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.y;
                displacement_correction.z = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.z;

                // Calculate the correction factor (see Wada et al. _B007 appendix).
                double inv_norm = 1.0 / cpu_vec3D_length_sq(displacement_correction);

                dot = cpu_vec3D_dot(n_A, displacement_correction);
                double alpha_A = 1.0 / (1.0 - dot * dot * inv_norm);

                dot = cpu_vec3D_dot(n_B, displacement_correction);
                double alpha_B = 1.0 / (1.0 - dot * dot * inv_norm);

                // Apply the correction to the current contact pointers.
                double particle_radius_inv = 1.0 / a_mon_A;
                n_A.x -= particle_radius_inv * alpha_A * displacement_correction.x;
                n_A.y -= particle_radius_inv * alpha_A * displacement_correction.y;
                n_A.z -= particle_radius_inv * alpha_A * displacement_correction.z;

                n_B.x -= particle_radius_inv * alpha_B * displacement_correction.x;
                n_B.y -= particle_radius_inv * alpha_B * displacement_correction.y;
                n_B.z -= particle_radius_inv * alpha_B * displacement_correction.z;

                cpu_vec3D_normalize(n_A);
                cpu_vec3D_normalize(n_B);

                // Track that the contact pointers need to be updated.
                update_contact_pointers = true;
            }

            // Update the contact pointers when neccessary.
            if (update_contact_pointers)
            {
                cpu_updateNormal(n_A, n_B, matrix_con, matrix_rot, i, j, Nmon);
            }

            // Calculate the SLIDING FORCE
            double sliding_modifier = 1.0; // TODO: Change this into a macro or remove it.
            double k_s = sliding_modifier * 8.0 * Gs * a0;

            // Sliding displacement
            vec3D displacement_zeta;
            double tmp_A = a_mon_A * cpu_vec3D_dot(n_A, n_c);
            double tmp_B = a_mon_B * cpu_vec3D_dot(n_B, n_c);

            // FIXME: Check this in Wada07, but: This formula is incorrect it should read r_i * n_i - r_j * n_j + (r_i + r_j) * n_c
            displacement_zeta.x = a_mon_A * n_A.x - a_mon_B * n_B.x - (tmp_A - tmp_B) * n_c.x;
            displacement_zeta.y = a_mon_A * n_A.y - a_mon_B * n_B.y - (tmp_A - tmp_B) * n_c.y;
            displacement_zeta.z = a_mon_A * n_A.z - a_mon_B * n_B.z - (tmp_A - tmp_B) * n_c.z;

            // FIXME: The following equations are using zeta and not zeta_0, see Wada07

            vec3D tmp;
            tmp.x = a_mon_B * n_B.x - a_mon_A * n_A.x;
            tmp.y = a_mon_B * n_B.y - a_mon_A * n_A.y;
            tmp.z = a_mon_B * n_B.z - a_mon_A * n_A.z;

            // FIXME: This can't be correct, right?
            double force_sliding = -k_s * cpu_vec3D_dot(displacement_zeta, tmp) / particle_distance;

            // FIXME: This might resolve itself once the bugs are fixed.
            // Clamp the sliding force to a maximum value
            if (abs(force_sliding) > 1.0e-10)
                force_sliding = 1e-10;

            force_new[i].x += force_sliding * n_c.x;
            force_new[i].y += force_sliding * n_c.y;
            force_new[i].z += force_sliding * n_c.z;

            // Calculate the SLIDING TORQUE
            vec3D torque_sliding;
            tmp = cpu_vec3D_cross(n_A, displacement_zeta); // FIXME: This uses zeta_0 and not zeta. See Wada07.

            torque_sliding.x = -a_mon_A * k_s * tmp.x;
            torque_sliding.y = -a_mon_A * k_s * tmp.y;
            torque_sliding.z = -a_mon_A * k_s * tmp.z;

            torque_new[i].x += torque_sliding.x;
            torque_new[i].y += torque_sliding.y;
            torque_new[i].z += torque_sliding.z;

            // Calculate the ROLLING TORQUE.
            vec3D displacement_xi;
            vec3D torque_rolling;

            double rolling_modifier = 1.0; // TODO: Change this into a macro or remove it.
            double k_r = rolling_modifier * 4.0 * Fc / R; // TODO: Remove /R here, because R is multiplied later anyways ?

            // The rolling displacement.
            displacement_xi.x = R * (n_A.x + n_B.x);
            displacement_xi.y = R * (n_A.y + n_B.y);
            displacement_xi.z = R * (n_A.z + n_B.z);

            tmp = cpu_vec3D_cross(n_A, displacement_xi);
            torque_rolling.x = -k_r * R * tmp.x;
            torque_rolling.y = -k_r * R * tmp.y;
            torque_rolling.z = -k_r * R * tmp.z;

            torque_new[i].x += torque_rolling.x;
            torque_new[i].y += torque_rolling.y;
            torque_new[i].z += torque_rolling.z;

            // Calculating the TWISTING FORCE.
            vec3D delta_omega_old, delta_omega_new, twisting_torque;
            double crit_twisting_displacement = 1.0 / (16.0 * PI);

            double twisting_modifier = 1.0; // TODO: Change this into a macro or remove it.
            double k_t = twisting_modifier * 16.0 / 3.0 * Gs * a0 * a0 * a0; // FIXME: This equation does not use Gs but rather the reduced shear modulus, see Wada07.

            // TODO: Check if this is true, when exactly does matrix_twist get updated?
            double twisting_displacement = matrix_twist[i * Nmon + j];
            double moment_inv_A = 1.0 / moment_A;
            double moment_inv_B = 1.0 / moment_B;

            // Store current contact normal.
            vec3D n_c_old = matrix_norm[i * Nmon + j];

            // Difference in angular momenta, ie change in 
            delta_omega_old.x = omega_A.x - omega_B.x;
            delta_omega_old.y = omega_A.y - omega_B.y;
            delta_omega_old.z = omega_A.z - omega_B.z;

            // update twisting displacement - use second order integration: omega^n+1 = omega^n + 0.5 (<delta_omega^n, n_c^n> + <delta_omega^n+1, n_c^n+1>)
            delta_omega_new.x = delta_omega_old.x + time_step * (moment_inv_A * torque_old[i].x - moment_inv_B * torque_old[j].x);
            delta_omega_new.y = delta_omega_old.y + time_step * (moment_inv_A * torque_old[i].y - moment_inv_B * torque_old[j].y);
            delta_omega_new.z = delta_omega_old.z + time_step * (moment_inv_A * torque_old[i].z - moment_inv_B * torque_old[j].z);

            twisting_displacement += 0.5 * time_step * (cpu_vec3D_dot(delta_omega_old, n_c_old) + cpu_vec3D_dot(delta_omega_new, n_c));

            // Clamps the twisting displacement to the critical displacement, but only in the positive direction...
            if (twisting_displacement > crit_twisting_displacement)
                twisting_displacement = crit_twisting_displacement;

            twisting_torque.x = k_t * twisting_displacement * n_c.x;
            twisting_torque.y = k_t * twisting_displacement * n_c.y;
            twisting_torque.z = k_t * twisting_displacement * n_c.z;

            torque_new[i].x -= twisting_torque.x;
            torque_new[i].y -= twisting_torque.y;
            torque_new[i].z -= twisting_torque.z;

            // Update the contact normal and compression lenghts.
            matrix_norm[i * Nmon + j].x = n_c.x;
            matrix_norm[i * Nmon + j].y = n_c.y;
            matrix_norm[i * Nmon + j].z = n_c.z;
            matrix_comp[i * Nmon + j] = compression_length;
        }
    }
}

/**
 * @brief Calculates the forces and torques updates the contact normal and compression lenght.
 * 
 * @param *pos_new:
 * @param *force_new:
 * @param *torque_old:
 * @param *dMdt_new:
 * @param *matrix_con:
 * @param *matrix_norm:
 * @param *omega:
 * @param *omega_tot:
 * @param *mag:
 * @param *matrix_rot:
 * @param *matrix_comp:
 * @param *matrix_twist:
 * @param *amon:
 * @param *moment:
 * @param *mat:
 * @param *matIDs:
 * @param B_ext:
 * @param N_mon:
 * @param time_step:
 */
__global__ void gpu_updateParticleInteraction(vec3D* pos_new, vec3D* force_new, vec3D* torque_old, vec3D* torque_new, vec3D* dMdt_new, vec3D* matrix_con,
    vec3D* matrix_norm, vec3D* omega, vec3D* omega_tot, vec3D* mag, quat* matrix_rot, double* matrix_comp, double* matrix_twist, double* amon, double* moment, material* mat, int* matIDs, vec3D B_ext, int Nmon, double time_step)
{
    // Determines the thread index which corresponds to a single monomer.
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // Skip threads that dont correspond to a particle.
    if (i < Nmon)
    {
        // Initialize some quantities to zero.
        force_new[i].x = 0;
        force_new[i].y = 0;
        force_new[i].z = 0;

        torque_new[i].x = 0;
        torque_new[i].y = 0;
        torque_new[i].z = 0;

        dMdt_new[i].x = 0;
        dMdt_new[i].y = 0;
        dMdt_new[i].z = 0;

        // Iterate over all monomer pairs.
        for (int j = 0; j < Nmon; j++)
        {
            if (i == j)
                continue;

            //current position
            vec3D pos_A = pos_new[i];
            vec3D pos_B = pos_new[j];

            //material IDs
            int mat_id_A = matIDs[i];
            int mat_id_B = matIDs[j];

            //monomer radius
            double a_mon_A = amon[i];
            double a_mon_B = amon[j];

            // Calculate the current contact pointer and monomer distance.
            vec3D n_c = gpu_vec3D_diff(pos_A, pos_B);
            double particle_distance = gpu_vec3D_length(n_c);
            gpu_vec3D_normalize(n_c);

            // TODO: What does this do?? It initializes force_new to an extremely small value proportional to n_c. And why do you need force_temp here?
            vec3D force_tmp;

            force_tmp.x = 1e-12 * n_c.x;
            force_tmp.y = 1e-12 * n_c.y;
            force_tmp.z = 1e-12 * n_c.z;

            force_new[i].x += force_tmp.x;
            force_new[i].y += force_tmp.y;
            force_new[i].z += force_tmp.z;

            // FIXME: Activate magnetization calculation.
            // TODO: Document the magnetization code.
            // ###############    MAGNETIZATION CALCULATIONS FROM HERE    ###############
            /*double chi_A = mat[mat_id_A].chi;

            if (abs(chi_A) > 0)
            {
                double len_B = gpu_vec3D_length(B_ext);
                double Vmon_A = 4. / 3. * PI * a_mon_A * a_mon_A * a_mon_A;
                double Vmon_B = 4. / 3. * PI * a_mon_B * a_mon_B * a_mon_B;

                double Msat_A = mat[mat_id_A].Msat;

                double tau_ss_A = mat[mat_id_A].tss;
                double tau_sl_A = mat[mat_id_A].tsl;

                vec3D mag_A = mag[i];
                vec3D mag_B = mag[j];

                vec3D omega_tot_A = omega_tot[i];

                vec3D Mfin, Mpara, Mperp;

                gpu_vec3D_set(Mfin, 0);
                gpu_vec3D_set(Mpara, 0);
                gpu_vec3D_set(Mperp, 0);

                double len_mag_A = gpu_vec3D_length(mag_A);
                double len_mag_B = gpu_vec3D_length(mag_B);

                if (abs(chi_A) > LIMIT_FER) //ferromagnetic material
                {
                    double len_omega = gpu_vec3D_length(omega_tot_A);

                    if (len_mag_A > 0) //rescale current magnetization
                    {
                        mag_A.x = Msat_A * mag_A.x / len_mag_A;
                        mag_A.y = Msat_A * mag_A.y / len_mag_A;
                        mag_A.z = Msat_A * mag_A.z / len_mag_A;
                    }

                    if (len_B > 0) //future magnetization
                    {
                        Mfin.x = Msat_A * B_ext.x / len_B;
                        Mfin.y = Msat_A * B_ext.y / len_B;
                        Mfin.z = Msat_A * B_ext.z / len_B;
                    }
                    else
                    {
                        if (len_omega > 0)
                        {
                            Mfin.x = Msat_A * omega_tot_A.x / len_omega;
                            Mfin.y = Msat_A * omega_tot_A.y / len_omega;
                            Mfin.z = Msat_A * omega_tot_A.z / len_omega;
                        }
                        else //just put Mfin in x-direction for now if B_ext == 0 and omega_tot == 0
                        {
                            Mfin.x = Msat_A;
                            Mfin.y = 0;
                            Mfin.z = 0;
                        }
                    }

                    double len_Mfin = gpu_vec3D_length(Mfin);
                    double dot_magA_Mfin = gpu_vec3D_dot(mag_A, Mfin);

                    Mpara.x = dot_magA_Mfin * Mfin.x / (len_Mfin * len_Mfin);
                    Mpara.y = dot_magA_Mfin * Mfin.y / (len_Mfin * len_Mfin);
                    Mpara.z = dot_magA_Mfin * Mfin.z / (len_Mfin * len_Mfin);

                    Mperp.x = mag_A.x - Mpara.x;
                    Mperp.y = mag_A.y - Mpara.y;
                    Mperp.z = mag_A.z - Mpara.z;
                }
                else //paramegnetic material
                {
                    double chi_factor_A = chi_A / (chi_A + 1.);

                    //induced moment
                    Mfin.x = chi_factor_A * B_ext.x / mu0;
                    Mfin.y = chi_factor_A * B_ext.y / mu0;
                    Mfin.z = chi_factor_A * B_ext.z / mu0;

                    //Barnett moment
                    Mfin.x += chi_A * omega_tot_A.x / PROD_BARR;
                    Mfin.y += chi_A * omega_tot_A.y / PROD_BARR;
                    Mfin.z += chi_A * omega_tot_A.z / PROD_BARR;

                    double len_Mfin = gpu_vec3D_length(Mfin);

                    if (len_mag_A > Msat_A) //rescale current mag. if saturation is reached
                    {
                        if (len_mag_A > 0)
                        {
                            mag_A.x = Msat_A * mag_A.x / len_mag_A;
                            mag_A.y = Msat_A * mag_A.y / len_mag_A;
                            mag_A.z = Msat_A * mag_A.z / len_mag_A;
                        }
                    }

                    /*if (len_Mfin > Msat_A) //rescale future  mag. if saturation is reached
                    {
                        if (len_Mfin > 0)
                        {
                            Mfin.x = Msat_A * Mfin.x / len_Mfin;
                            Mfin.y = Msat_A * Mfin.y / len_Mfin;
                            Mfin.z = Msat_A * Mfin.z / len_Mfin;
                        }
                    }*/

                    /*double dot_magA_Mfin = gpu_vec3D_dot(mag_A, Mfin);

                    if (len_Mfin > 0)
                    {
                        Mpara.x = dot_magA_Mfin * Mfin.x / (len_Mfin * len_Mfin);
                        Mpara.y = dot_magA_Mfin * Mfin.y / (len_Mfin * len_Mfin);
                        Mpara.z = dot_magA_Mfin * Mfin.z / (len_Mfin * len_Mfin);

                        Mperp.x = mag_A.x - Mpara.x;
                        Mperp.y = mag_A.y - Mpara.y;
                        Mperp.z = mag_A.z - Mpara.z;
                    }
                    else
                    {
                        Mpara.x = 0;
                        Mpara.y = 0;
                        Mpara.z = 0;

                        Mperp.x = 0;
                        Mperp.y = 0;
                        Mperp.z = 0;
                    }
                }

                //solve Bloch equation
                dMdt_new[i].x = 0;
                dMdt_new[i].y = 0;
                dMdt_new[i].z = 0;

                if (len_B > 0)
                {
                    vec3D cross = gpu_vec3D_cross(mag_A, B_ext);

                    dMdt_new[i].x += GAMMA_E * cross.x;
                    dMdt_new[i].y += GAMMA_E * cross.y;
                    dMdt_new[i].z += GAMMA_E * cross.z;
                }

                if (tau_sl_A > 0)
                {
                    dMdt_new[i].x -= (mag_A.x - Mpara.x) / tau_sl_A;
                    dMdt_new[i].y -= (mag_A.y - Mpara.y) / tau_sl_A;
                    dMdt_new[i].z -= (mag_A.z - Mpara.z) / tau_sl_A;
                }

                if (tau_ss_A > 0)
                {
                    dMdt_new[i].x -= Mperp.x / tau_ss_A;
                    dMdt_new[i].y -= Mperp.y / tau_ss_A;
                    dMdt_new[i].z -= Mperp.z / tau_ss_A;
                }

                if (len_mag_A * len_mag_B > 0)
                {
                    vec3D mu_A, mu_B;

                    mu_A.x = mag_A.x * Vmon_A;
                    mu_A.y = mag_A.y * Vmon_A;
                    mu_A.z = mag_A.z * Vmon_A;

                    mu_B.x = mag_B.x * Vmon_B;
                    mu_B.y = mag_B.y * Vmon_B;
                    mu_B.z = mag_B.z * Vmon_B;

                    vec3D force_D, torque_D, torque_EdH;
                    vec3D torque_B = gpu_vec3D_cross(mu_A, B_ext);

                    vec3D vec_d = gpu_vec3D_diff(pos_B, pos_A); //todo: check for correct direction
                    double d = gpu_vec3D_length(vec_d);
                    double d2 = d * d;
                    double d3 = d2 * d;
                    double d4 = d2 * d2;
                    double d5 = d2 * d3;

                    //double fD = (3. * mu0) / (PIx4 * d5);
                    double fD = (3. * mu0) / (PIx4 * d4);
                    double tD = (mu0 / PIx4);
                    double tEdH = (Vmon_A / GAMMA_E);

                    double dot_muA_d = gpu_vec3D_dot(mu_A, vec_d);
                    double dot_muB_d = gpu_vec3D_dot(mu_B, vec_d);
                    double dot_muA_muB = gpu_vec3D_dot(mu_A, mu_B);

                    vec3D cross_d_muB = gpu_vec3D_cross(vec_d, mu_B);
                    vec3D cross_muA_muB = gpu_vec3D_cross(mu_A, mu_B);

                    //force_D.x = fD * (5. * vec_d.x * dot_muA_d * dot_muB_d / d2 - vec_d.x * dot_muA_muB - mu_A.x * dot_muB_d - mu_B.x * dot_muA_d);
                    //force_D.y = fD * (5. * vec_d.y * dot_muA_d * dot_muB_d / d2 - vec_d.y * dot_muA_muB - mu_A.y * dot_muB_d - mu_B.y * dot_muA_d);
                    //force_D.z = fD * (5. * vec_d.z * dot_muA_d * dot_muB_d / d2 - vec_d.z * dot_muA_muB - mu_A.z * dot_muB_d - mu_B.z * dot_muA_d);

                    force_D.x = fD * (-5. * vec_d.x * dot_muA_d * dot_muB_d + vec_d.x * dot_muA_muB + mu_A.x * dot_muB_d + mu_B.x * dot_muA_d);
                    force_D.y = fD * (-5. * vec_d.y * dot_muA_d * dot_muB_d + vec_d.y * dot_muA_muB + mu_A.y * dot_muB_d + mu_B.y * dot_muA_d);
                    force_D.z = fD * (-5. * vec_d.z * dot_muA_d * dot_muB_d + vec_d.z * dot_muA_muB + mu_A.z * dot_muB_d + mu_B.z * dot_muA_d);

                    torque_D.x = tD * 3. * dot_muB_d * cross_d_muB.x / d5 - cross_muA_muB.x / d3;
                    torque_D.y = tD * 3. * dot_muB_d * cross_d_muB.y / d5 - cross_muA_muB.y / d3;
                    torque_D.z = tD * 3. * dot_muB_d * cross_d_muB.z / d5 - cross_muA_muB.z / d3;

                    torque_EdH.x = -tEdH * dMdt_new[i].x;
                    torque_EdH.y = -tEdH * dMdt_new[i].y;
                    torque_EdH.z = -tEdH * dMdt_new[i].z;

                    force_new[i].x += force_D.x;
                    force_new[i].y += force_D.y;
                    force_new[i].z += force_D.z;

                    torque_new[i].x += torque_D.x;
                    torque_new[i].y += torque_D.y;
                    torque_new[i].z += torque_D.z;

                    torque_new[i].x += torque_EdH.x;
                    torque_new[i].y += torque_EdH.y;
                    torque_new[i].z += torque_EdH.z;

                    torque_new[i].x += torque_B.x;
                    torque_new[i].y += torque_B.y;
                    torque_new[i].z += torque_B.z;
                }
            }*/

            // Skip unconnected monomer pairs.
            if (matrix_comp[i * Nmon + j] == -1.)
                continue;

            // This tracks if any of the displacements have reached a critical value.
            bool update_contact_pointers = false;

            // Determine monomer indices in the contact matrix.
            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

            // Calculate required quantities for force calculations.
            vec3D omega_A = omega[i];
            vec3D omega_B = omega[j];

            double R = (a_mon_A * a_mon_B) / (a_mon_A + a_mon_B);

            double moment_A = moment[i];
            double moment_B = moment[j];

            double nu_A = mat[mat_id_A].nu;
            double nu_B = mat[mat_id_B].nu;

            double E_A = mat[mat_id_A].E;
            double E_B = mat[mat_id_B].E;

            double T_vis_A = mat[mat_id_A].tvis;
            double T_vis_B = mat[mat_id_B].tvis;
            double T_vis = 0.5 * (T_vis_A + T_vis_B);

            double Es = (1 - nu_A * nu_A) / E_A + (1 - nu_B * nu_B) / E_B;
            Es = 1. / Es;

            double G_A = 0.5 * E_A / (1. + nu_A);
            double G_B = 0.5 * E_B / (1. + nu_B);

            double Gs = (1. - nu_A * nu_A) / G_A + (1. - nu_B * nu_B) / G_B;
            Gs = 1. / Gs;

            double gamma_A = mat[mat_id_A].gamma;
            double gamma_B = mat[mat_id_B].gamma;
            double gamma = gamma_A + gamma_B - 2. / (1. / gamma_A + 1. / gamma_B);

            double a0 = pow(9 * PI * gamma * R * R / Es, 1. / 3.);
            double delta_c = 0.5 * a0 * a0 / (R * pow(6.0, 1.0 / 3.0));

            // Calculate current contact pointers.
            vec3D n_A = matrix_con[index_A];
            vec3D n_B = matrix_con[index_B];
            vec3D delta_n = gpu_vec3D_diff(n_A, n_B);

            // ###############      FORCE CALCULATIONS FROM HERE      ###############
            double compression_length = a_mon_A + a_mon_B - particle_distance;
            double contact_radius = gpu_getJKRContactRadius(compression_length, a0, R);
            
            // Calculate the NORMAL FORCE.
            double Fc = 3 * PI * gamma * R;
            double force_elastic = 4.0 * Fc * (pow(contact_radius / a0, 3.0) - pow(contact_radius / a0, 3.0 / 2.0));

            force_new[i].x += force_elastic * n_c.x;
            force_new[i].y += force_elastic * n_c.y;
            force_new[i].z += force_elastic * n_c.z;

            // TODO: Research this.
            // Calculate the DAMPING FORCE.
            double old_compression_length = matrix_comp[i * Nmon + j];
            double vis_damp_const = 2.0 * T_vis / (nu_A * nu_B) * Es;
            double delta_dot = (compression_length - old_compression_length) / time_step;
            double force_damp = vis_damp_const * contact_radius * delta_dot;

            force_new[i].x += force_damp * n_c.x;
            force_new[i].y += force_damp * n_c.y;
            force_new[i].z += force_damp * n_c.z;

            // Determine sliding and rolling displacement.
            double dot = gpu_vec3D_dot(delta_n, n_c);

            // TODO: Is a_mon_A correct here?
            vec3D displacement;
            displacement.x = a_mon_A * (delta_n.x - dot * n_c.x);
            displacement.y = a_mon_A * (delta_n.y - dot * n_c.y);
            displacement.z = a_mon_A * (delta_n.z - dot * n_c.z);

            double displacement_norm = gpu_vec3D_length(displacement);

            double crit_sliding_displacement_modifier = 1.0;
            double crit_sliding_displacement = crit_sliding_displacement_modifier * (2.0 - nu_A) / (16.0 * PI) * a0;

            // When the inelastic sliding regime is regime is reached.
            if (displacement_norm > crit_sliding_displacement)
            {
                // Determine correction of contact pointers
                vec3D displacement_correction;
                displacement_correction.x = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.x;
                displacement_correction.y = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.y;
                displacement_correction.z = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.z;

                // Calculate correction factor (see Wada et al. 2007 appendix for details)
                double inv_norm = 1.0 / gpu_vec3D_length_sq(displacement_correction);

                dot = gpu_vec3D_dot(n_A, displacement_correction);
                double alpha_A = 1.0 / (1.0 - dot * dot * inv_norm);

                dot = gpu_vec3D_dot(n_A, displacement_correction);
                double alpha_B = 1.0 / (1.0 - dot * dot * inv_norm);

                // Apply contact pointer corrections to the current contact pointers.
                double particle_radius_inv = 1.0 / a_mon_A;
                n_A.x -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.x;
                n_A.y -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.y;
                n_A.z -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.z;

                n_B.x += 0.5 * particle_radius_inv * alpha_B * displacement_correction.x;
                n_B.y += 0.5 * particle_radius_inv * alpha_B * displacement_correction.y;
                n_B.z += 0.5 * particle_radius_inv * alpha_B * displacement_correction.z;

                gpu_vec3D_normalize(n_A);
                gpu_vec3D_normalize(n_B);

                // Track that the contact pointers need to be updates.
                update_contact_pointers = true;
            }

            double crit_rolling_displacement = 0.5 * (mat[mat_id_A].xi + mat[mat_id_B].xi);

            displacement.x = R * (n_A.x + n_B.x);
            displacement.y = R * (n_A.y + n_B.y);
            displacement.z = R * (n_A.z + n_B.z);
            displacement_norm = gpu_vec3D_length(displacement);

            // When the inelastic rolling regime is reached.
            if (displacement_norm > crit_rolling_displacement)
            {
                // Determine correction of contact pointers
                vec3D displacement_correction;
                displacement_correction.x = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.x;
                displacement_correction.y = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.y;
                displacement_correction.z = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.z;

                // Calculate correction factor (see Wada et al. _B007 appendix for details)
                double inv_norm = 1.0 / gpu_vec3D_length_sq(displacement_correction);

                dot = gpu_vec3D_dot(n_A, displacement_correction);
                double alpha_A = 1.0 / (1.0 - dot * dot * inv_norm);

                dot = gpu_vec3D_dot(n_B, displacement_correction);
                double alpha_B = 1.0 / (1.0 - dot * dot * inv_norm);

                // Apply the correction to the current contact pointers.
                double particle_radius_inv = 1.0 / a_mon_A;
                n_A.x -= particle_radius_inv * alpha_A * displacement_correction.x;
                n_A.y -= particle_radius_inv * alpha_A * displacement_correction.y;
                n_A.z -= particle_radius_inv * alpha_A * displacement_correction.z;

                n_B.x -= particle_radius_inv * alpha_B * displacement_correction.x;
                n_B.y -= particle_radius_inv * alpha_B * displacement_correction.y;
                n_B.z -= particle_radius_inv * alpha_B * displacement_correction.z;

                gpu_vec3D_normalize(n_A);
                gpu_vec3D_normalize(n_B);

                // Track that the contact pointers need to be updated.
                update_contact_pointers = true;
            }

            // Update the contact pointers when neccessary.
            if (update_contact_pointers)
            {
                gpu_updateNormal(n_A, n_B, matrix_con, matrix_rot, i, j, Nmon);
            }

            // Calculate the SLIDING FORCE
            double sliding_modifier = 1.0; // TODO: Make this a macro or maybe just remove it.
            double k_s = sliding_modifier * 8.0 * Gs * a0;

            // Sliding displacement
            vec3D displacement_zeta;
            double tmp_A = a_mon_A * gpu_vec3D_dot(n_A, n_c);
            double tmp_B = a_mon_B * gpu_vec3D_dot(n_B, n_c);

            // FIXME: Check this in Wada07, but: This formula is incorrect it should read r_i * n_i - r_j * n_j + (r_i + r_j) * n_c
            displacement_zeta.x = a_mon_A * n_A.x - a_mon_B * n_B.x - (tmp_A - tmp_B) * n_c.x;
            displacement_zeta.y = a_mon_A * n_A.y - a_mon_B * n_B.y - (tmp_A - tmp_B) * n_c.y;
            displacement_zeta.z = a_mon_A * n_A.z - a_mon_B * n_B.z - (tmp_A - tmp_B) * n_c.z;

            // FIXME: The following equations are using zeta and not zeta_0, see Wada07
            
            vec3D tmp;
            tmp.x = a_mon_B * n_B.x - a_mon_A * n_A.x;
            tmp.y = a_mon_B * n_B.y - a_mon_A * n_A.y;
            tmp.z = a_mon_B * n_B.z - a_mon_A * n_A.z;

            // FIXME: I dont understand what this does, but it cant be right?
            double force_sliding = -k_s * gpu_vec3D_dot(displacement_zeta, tmp) / particle_distance;

            // FIXME: This is a hack, might resolve itself once the formulaes are correct.
            // Clamp the sliding force to a maximum value
            if (abs(force_sliding) > 1.0e-10)
                force_sliding = 1e-10;

            force_new[i].x += force_sliding * n_c.x;
            force_new[i].y += force_sliding * n_c.y;
            force_new[i].z += force_sliding * n_c.z;

            // Calculate the SLIDING FORCE
            vec3D torque_sliding;
            tmp = gpu_vec3D_cross(n_A, displacement_zeta); // FIXME: This uses zeta_0 and not zeta. See Wada07

            torque_sliding.x = -a_mon_A * k_s * tmp.x;
            torque_sliding.y = -a_mon_A * k_s * tmp.y;
            torque_sliding.z = -a_mon_A * k_s * tmp.z;

            torque_new[i].x += torque_sliding.x;
            torque_new[i].y += torque_sliding.y;
            torque_new[i].z += torque_sliding.z;

            // Calculate the ROLLING TORQUE
            vec3D displacement_xi;
            vec3D torque_rolling;

            double rolling_modifier = 1.0; // TODO: Change this into a macro or just remove it.
            double k_r = rolling_modifier * 4.0 * Fc / R; // Remove /R from here because it gets devided out later anyway.

            // The rolling displacement.
            displacement_xi.x = R * (n_A.x + n_B.x);
            displacement_xi.y = R * (n_A.y + n_B.y);
            displacement_xi.z = R * (n_A.z + n_B.z);

            tmp = gpu_vec3D_cross(n_A, displacement_xi);
            torque_rolling.x = -k_r * R * tmp.x;
            torque_rolling.y = -k_r * R * tmp.y;
            torque_rolling.z = -k_r * R * tmp.z;

            torque_new[i].x += torque_rolling.x;
            torque_new[i].y += torque_rolling.y;
            torque_new[i].z += torque_rolling.z;

            // Calculate the TWISTING FORCE
            vec3D delta_omega_old, delta_omega_new, twisting_torque;
            double crit_twisting_displacement = 1.0 / (16.0 * PI);

            double twisting_modifier = 1.0; // TODO: Change this into a macro or just remove it.
            double k_t = twisting_modifier * 16.0 / 3.0 * Gs * a0 * a0 * a0;
            
            // TODO: Check if this is true, when exactly does matrix_twist get modified?
            double twisting_displacement = matrix_twist[i * Nmon + j];
            double moment_inv_A = 1.0 / moment_A;
            double moment_inv_B = 1.0 / moment_B;

            // Store current contact normal.
            vec3D n_c_old = matrix_norm[i * Nmon + j];

            // Difference in angular momenta, ie change in twisting displacement.
            delta_omega_old.x = omega_A.x - omega_B.x;
            delta_omega_old.y = omega_A.y - omega_B.y;
            delta_omega_old.z = omega_A.z - omega_B.z;

            // update twisting displacement - use second order integration: omega^n+1 = omega^n + 0.5 (<delta_omega^n, n_c^n> + <delta_omega^n+1, n_c^n+1>)
            delta_omega_new.x = delta_omega_old.x + time_step * (moment_inv_A * torque_old[i].x - moment_inv_B * torque_old[j].x);
            delta_omega_new.y = delta_omega_old.y + time_step * (moment_inv_A * torque_old[i].y - moment_inv_B * torque_old[j].y);
            delta_omega_new.z = delta_omega_old.z + time_step * (moment_inv_A * torque_old[i].z - moment_inv_B * torque_old[j].z);

            twisting_displacement += 0.5 * time_step * (gpu_vec3D_dot(delta_omega_old, n_c_old) + gpu_vec3D_dot(delta_omega_new, n_c));

            // Clamps the twisting displacement // TODO: the clamping is only in the positive direction...
            if (twisting_displacement > crit_twisting_displacement)
                twisting_displacement = crit_twisting_displacement;

            twisting_torque.x = k_t * twisting_displacement * n_c.x;
            twisting_torque.y = k_t * twisting_displacement * n_c.y;
            twisting_torque.z = k_t * twisting_displacement * n_c.z;

            torque_new[i].x -= twisting_torque.x;
            torque_new[i].y -= twisting_torque.y;
            torque_new[i].z -= twisting_torque.z;/**/

            // Update the contact normal and compression lenghts.
            matrix_norm[i * Nmon + j].x = n_c.x;
            matrix_norm[i * Nmon + j].y = n_c.y;
            matrix_norm[i * Nmon + j].z = n_c.z;
            matrix_comp[i * Nmon + j] = compression_length;
        }
    }
}

