#pragma once
#include "typedefs.h"
#include <cmath>

#ifndef CPHYSICS
#define CPHYSICS

// Function to perform Depth-First Search (DFS) using recursion
inline void dfs(int node, int n, const double * matrix, int* cluster, int currentCluster) 
{
    // Assign the current cluster ID to the node
    cluster[node] = currentCluster;

    // Explore all neighbors of the node
    for (int i = 0; i < n; ++i)
    {
        // Check if there's an edge and the node has not been assigned to a cluster
        if (matrix[node * n + i] != -1 && cluster[i] == -1)
            dfs(i, n, matrix, cluster, currentCluster);
    }
}

inline void findConnectedComponents(int Nmon, const double* matrix, int* cluster)
{
    int currentCluster = 0;  // Cluster ID counter

    // Iterate through all nodes
    for (int i = 0; i < Nmon; i++)
    {
        if (cluster[i] == -1)
        {  // If node hasn't been assigned to a cluster
            dfs(i, Nmon, matrix, cluster, currentCluster);  // Perform DFS
            currentCluster++;  // Move to the next cluster ID
        }
    }
}

//update sticking
inline void updateNeighbourhoodRelations(vec3D* pos, vec3D* matrix_con, vec3D* matrix_norm, quat* matrix_rot, double* matrix_comp, double* matrix_twist, double* amon, material* mat, int* matIDs, int Nmon)
{
    for (int i = 0; i < Nmon; i++)
    {
        for (int j = 1; j < Nmon; j++)
        {
            if (i == j)
                continue;

            vec3D pos_A = pos[i];
            vec3D pos_B = pos[j];

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
            double distance = vec3D_distance(pos_A, pos_B);

            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

            /*if (index_A > 2 * Nmon * Nmon)
                int tt = 0;

            if (index_B > 2 * Nmon * Nmon)
                int tt = 0;*/

            if (distance < contact_distance)
            {
                //new connection is established
                if (matrix_comp[i * Nmon + j] == -1.)
                {
                    vec3D n = vec3D_get_normal(pos_A, pos_B);

                    // init. contact pointer
                    matrix_con[index_A].x = -n.x;
                    matrix_con[index_A].y = -n.y;
                    matrix_con[index_A].z = -n.z;

                    matrix_con[index_B].x = n.x;
                    matrix_con[index_B].y = n.y;
                    matrix_con[index_B].z = n.z;

                    matrix_rot[index_A].e0 = 1;
                    matrix_rot[index_A].e1 = 0;
                    matrix_rot[index_A].e2 = 0;
                    matrix_rot[index_A].e3 = 0;

                    matrix_rot[index_B].e0 = 1;
                    matrix_rot[index_B].e1 = 0;
                    matrix_rot[index_B].e2 = 0;
                    matrix_rot[index_B].e3 = 0;

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
                matrix_con[index_A].x = 0;
                matrix_con[index_A].y = 0;
                matrix_con[index_A].z = 0;

                matrix_con[index_B].x = 0;
                matrix_con[index_B].y = 0;
                matrix_con[index_B].z = 0;

                matrix_comp[i * Nmon + j] = -1.; //mark as disconnected
                matrix_twist[i * Nmon + j] = 0;
            }
        }
    }
}

inline void updateNormal(vec3D& n_A, vec3D& n_B, vec3D* matrix_con, quat* matrix_rot, int i, int j, int Nmon)
{
    int index_A = 0 * Nmon * Nmon + i * Nmon + j;
    int index_B = 1 * Nmon * Nmon + i * Nmon + j;

    quat rot_A = matrix_rot[index_A];
    quat rot_B = matrix_rot[index_B];

    vec3D init_n_A, init_n_B;

    init_n_A.x = 2.0 * ((0.5 - rot_A.e2 * rot_A.e2 - rot_A.e3 * rot_A.e3) * n_A.x + (rot_A.e1 * rot_A.e2 + rot_A.e3 * rot_A.e0) * n_A.y + (rot_A.e1 * rot_A.e3 - rot_A.e2 * rot_A.e0) * n_A.z);
    init_n_A.y = 2.0 * ((rot_A.e1 * rot_A.e2 - rot_A.e3 * rot_A.e0) * n_A.x + (0.5 - rot_A.e1 * rot_A.e1 - rot_A.e3 * rot_A.e3) * n_A.y + (rot_A.e2 * rot_A.e3 + rot_A.e1 * rot_A.e0) * n_A.z);
    init_n_A.z = 2.0 * ((rot_A.e1 * rot_A.e3 + rot_A.e2 * rot_A.e0) * n_A.x + (rot_A.e2 * rot_A.e3 - rot_A.e1 * rot_A.e0) * n_A.y + (0.5 - rot_A.e1 * rot_A.e1 - rot_A.e2 * rot_A.e2) * n_A.z);

    init_n_B.x = 2.0 * ((0.5 - rot_B.e2 * rot_B.e2 - rot_B.e3 * rot_B.e3) * n_B.x + (rot_B.e1 * rot_B.e2 + rot_B.e3 * rot_B.e0) * n_B.y + (rot_B.e1 * rot_B.e3 - rot_B.e2 * rot_B.e0) * n_B.z);
    init_n_B.y = 2.0 * ((rot_B.e1 * rot_B.e2 - rot_B.e3 * rot_B.e0) * n_B.x + (0.5 - rot_B.e1 * rot_B.e1 - rot_B.e3 * rot_B.e3) * n_B.y + (rot_B.e2 * rot_B.e3 + rot_B.e1 * rot_B.e0) * n_B.z);
    init_n_B.z = 2.0 * ((rot_B.e1 * rot_B.e3 + rot_B.e2 * rot_B.e0) * n_B.x + (rot_B.e2 * rot_B.e3 - rot_B.e1 * rot_B.e0) * n_B.y + (0.5 - rot_B.e1 * rot_B.e1 - rot_B.e2 * rot_B.e2) * n_B.z);

    matrix_con[index_A].x = init_n_A.x;
    matrix_con[index_A].y = init_n_A.y;
    matrix_con[index_A].z = init_n_A.z;

    matrix_con[index_B].x = init_n_B.x;
    matrix_con[index_B].y = init_n_B.y;
    matrix_con[index_B].z = init_n_B.z;
}

inline void predictor(vec3D* pos_old, vec3D* pos_new, vec3D* force_old, vec3D* vel, double* mass, double time_step, int Nmon)
{
    for (int i = 0; i < Nmon; i++)
    {
        double mass_inv = 1. / mass[i];

        pos_new[i].x = pos_old[i].x + time_step * vel[i].x + 0.5 * mass_inv * time_step * time_step * force_old[i].x;
        pos_new[i].y = pos_old[i].y + time_step * vel[i].y + 0.5 * mass_inv * time_step * time_step * force_old[i].y;
        pos_new[i].z = pos_old[i].z + time_step * vel[i].z + 0.5 * mass_inv * time_step * time_step * force_old[i].z;
    }
}



inline void corrector(vec3D* force_old, vec3D* force_new, vec3D* torque_old, vec3D* torque_new, vec3D* dMdt_old, vec3D* dMdt_new, vec3D* vel, vec3D* omega, vec3D* omega_tot, vec3D* mag, double* mass, double* moment, material* mat, int* matIDs, double time_step, int Nmon)
{
    for (int i = 0; i < Nmon; i++)
    {
        double mass_inv = 1. / mass[i];
        double moment_of_inertia_inv = 1. / moment[i];
        
        vec3D acc;//acceleration

        acc.x = 0.5 * mass_inv * (force_new[i].x + force_old[i].x);
        acc.y = 0.5 * mass_inv * (force_new[i].y + force_old[i].y);
        acc.z = 0.5 * mass_inv * (force_new[i].z + force_old[i].z);

        //vel[i].x += 0.5 * mass_inv * time_step * (force_new[i].x + force_old[i].x);
        //vel[i].y += 0.5 * mass_inv * time_step * (force_new[i].y + force_old[i].y);
        //vel[i].z += 0.5 * mass_inv * time_step * (force_new[i].z + force_old[i].z);

        vel[i].x += acc.x * time_step;
        vel[i].y += acc.y * time_step;
        vel[i].z += acc.z * time_step;

        omega[i].x += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].x + torque_old[i].x);
        omega[i].y += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].y + torque_old[i].y);
        omega[i].z += 0.5 * moment_of_inertia_inv * time_step * (torque_new[i].z + torque_old[i].z);

        double v_sq = vec3D_length_sq(vel[i]);
        vec3D cross = vec3D_cross(vel[i], acc);

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

        int mat_id = matIDs[i];
        double chi = mat[mat_id].chi;
        double Msat = mat[mat_id].Msat;

        if (abs(chi) > 0) //magnetic material
        {
            mag[i].x += 0.5 * time_step * (dMdt_new[i].x + dMdt_old[i].x);
            mag[i].y += 0.5 * time_step * (dMdt_new[i].y + dMdt_old[i].y);
            mag[i].z += 0.5 * time_step * (dMdt_new[i].z + dMdt_old[i].z);

            double len_mag = vec3D_length(mag[i]);

            if (len_mag > 0)
            {
                if (abs(chi) > LIMIT_FER) //re-scale ferromagnetic material to Msat
                {
                    mag[i].x = Msat * mag[i].x / len_mag;
                    mag[i].y = Msat * mag[i].y / len_mag;
                    mag[i].z = Msat * mag[i].z / len_mag;
                }
                else
                {
                    if (len_mag > Msat) //rescale current mag. if saturation is reached
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

inline void updateContacts(vec3D* omega, vec3D* omega_tot, vec3D* torque, vec3D* mag, quat* matrix_rot, double* matrix_comp, double* moment, int Nmon, double time_step)
{
    for (int i = 0; i < Nmon; i++)
    {
        vec3D omega_A;
        vec3D omega_A_dot;
        double moment_A_inv = 1.0 / moment[i];

        quat e_dot;
        quat e_ddot;
        double temp = 0;

        omega_A_dot.x = moment_A_inv* torque[i].x;
        omega_A_dot.y = moment_A_inv* torque[i].y;
        omega_A_dot.z = moment_A_inv* torque[i].z;

        omega_A.x = omega_tot[i].x;
        omega_A.y = omega_tot[i].y;
        omega_A.z = omega_tot[i].z;

        double len_mag = vec3D_length(mag[i]);
        double len_omega_A = vec3D_length(omega_A);
        //double len_omega_A_dot = vec3D_length(omega_A_dot);//len_omega_A_dot

        if(len_mag * len_omega_A  > 0)
        {
            quat q_mag;
            vec3D tmp_mag;

            tmp_mag.x = mag[i].x;
            tmp_mag.y = mag[i].y;
            tmp_mag.z = mag[i].z;

            vec3D_normalize(tmp_mag);

            q_mag.e0 = 0;
            q_mag.e1 = tmp_mag.x;
            q_mag.e2 = tmp_mag.y;
            q_mag.e3 = tmp_mag.z;

            e_dot.e0 = -0.5 * (q_mag.e1 * omega_A.x + q_mag.e2 * omega_A.y + q_mag.e3 * omega_A.z);
            e_dot.e1 = 0.5 * (q_mag.e0 * omega_A.x - q_mag.e2 * omega_A.z + q_mag.e3 * omega_A.y);
            e_dot.e2 = 0.5 * (q_mag.e0 * omega_A.y - q_mag.e3 * omega_A.x + q_mag.e1 * omega_A.z);
            e_dot.e3 = 0.5 * (q_mag.e0 * omega_A.z - q_mag.e1 * omega_A.y + q_mag.e2 * omega_A.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (q_mag.e0 * vec3D_length_sq(omega_A) + 2.0 * (q_mag.e1 * omega_A_dot.x + q_mag.e2 * omega_A_dot.y + q_mag.e3 * omega_A_dot.z));
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

            double len_tmp = vec3D_length(tmp_mag);

            /*if (len_tmp > 0)
            {
                mag[i].x = len_mag * tmp_mag.x / len_tmp;
                mag[i].y = len_mag * tmp_mag.y / len_tmp;
                mag[i].z = len_mag * tmp_mag.z / len_tmp;
            }*/
        }

        for (int j = 1; j < Nmon; j++)
        {
            if (i == j)
                continue;

            if (matrix_comp[i * Nmon + j] == -1.)
                continue;


            //--- > Rotate contact pointer
       
            vec3D omega_B;
            vec3D omega_B_dot;
            

            omega_A_dot.x = moment_A_inv* torque[i].x;
            omega_A_dot.y = moment_A_inv* torque[i].y;
            omega_A_dot.z = moment_A_inv* torque[i].z;

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

            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

            quat rot_A = matrix_rot[index_A];
            quat rot_B = matrix_rot[index_B];

            // first particle
            e_dot.e0 = -0.5 * (rot_A.e1 * omega_A.x + rot_A.e2 * omega_A.y + rot_A.e3 * omega_A.z);
            e_dot.e1 = 0.5 * (rot_A.e0 * omega_A.x - rot_A.e2 * omega_A.z + rot_A.e3 * omega_A.y);
            e_dot.e2 = 0.5 * (rot_A.e0 * omega_A.y - rot_A.e3 * omega_A.x + rot_A.e1 * omega_A.z);
            e_dot.e3 = 0.5 * (rot_A.e0 * omega_A.z - rot_A.e1 * omega_A.y + rot_A.e2 * omega_A.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (rot_A.e0 * vec3D_length_sq(omega_A) + 2.0 * (rot_A.e1 * omega_A_dot.x + rot_A.e2 * omega_A_dot.y + rot_A.e3 * omega_A_dot.z));
            e_ddot.e1 = temp * omega_A.x + 0.5 * (rot_A.e0 * omega_A_dot.x - rot_A.e2 * omega_A_dot.z + rot_A.e3 * omega_A_dot.y);
            e_ddot.e2 = temp * omega_A.y + 0.5 * (rot_A.e0 * omega_A_dot.y - rot_A.e3 * omega_A_dot.x + rot_A.e1 * omega_A_dot.z);
            e_ddot.e3 = temp * omega_A.z + 0.5 * (rot_A.e0 * omega_A_dot.z - rot_A.e1 * omega_A_dot.y + rot_A.e2 * omega_A_dot.x);

            rot_A.e0 += time_step * e_dot.e0 + 0.5 * time_step * time_step * e_ddot.e0;
            rot_A.e1 += time_step * e_dot.e1 + 0.5 * time_step * time_step * e_ddot.e1;
            rot_A.e2 += time_step * e_dot.e2 + 0.5 * time_step * time_step * e_ddot.e2;
            rot_A.e3 += time_step * e_dot.e3 + 0.5 * time_step * time_step * e_ddot.e3;

            // second particle
            e_dot.e0 = -0.5 * (rot_B.e1 * omega_B.x + rot_B.e2 * omega_B.y + rot_B.e3 * omega_B.z);
            e_dot.e1 = 0.5 * (rot_B.e0 * omega_B.x - rot_B.e2 * omega_B.z + rot_B.e3 * omega_B.y);
            e_dot.e2 = 0.5 * (rot_B.e0 * omega_B.y - rot_B.e3 * omega_B.x + rot_B.e1 * omega_B.z);
            e_dot.e3 = 0.5 * (rot_B.e0 * omega_B.z - rot_B.e1 * omega_B.y + rot_B.e2 * omega_B.x);

            temp = 0.5 * e_dot.e0;

            e_ddot.e0 = -0.25 * (rot_B.e0 * vec3D_length_sq(omega_B) + 2.0 * (rot_B.e1 * omega_B_dot.x + rot_B.e2 * omega_B_dot.y + rot_B.e3 * omega_B_dot.z));
            e_ddot.e1 = temp * omega_B.x + 0.5 * (rot_B.e0 * omega_B_dot.x - rot_B.e2 * omega_B_dot.z + rot_B.e3 * omega_B_dot.y);
            e_ddot.e2 = temp * omega_B.y + 0.5 * (rot_B.e0 * omega_B_dot.y - rot_B.e3 * omega_B_dot.x + rot_B.e1 * omega_B_dot.z);
            e_ddot.e3 = temp * omega_B.z + 0.5 * (rot_B.e0 * omega_B_dot.z - rot_B.e1 * omega_B_dot.y + rot_B.e2 * omega_B_dot.x);

            rot_B.e0 += time_step * e_dot.e0 + 0.5 * time_step * time_step * e_ddot.e0;
            rot_B.e1 += time_step * e_dot.e1 + 0.5 * time_step * time_step * e_ddot.e1;
            rot_B.e2 += time_step * e_dot.e2 + 0.5 * time_step * time_step * e_ddot.e2;
            rot_B.e3 += time_step * e_dot.e3 + 0.5 * time_step * time_step * e_ddot.e3;

            quat_normalize(rot_A);
            quat_normalize(rot_B);

            matrix_rot[index_A].e0 = rot_A.e0;
            matrix_rot[index_A].e1 = rot_A.e1;
            matrix_rot[index_A].e2 = rot_A.e2;
            matrix_rot[index_A].e3 = rot_A.e3;

            matrix_rot[index_B].e0 = rot_B.e0;
            matrix_rot[index_B].e1 = rot_B.e1;
            matrix_rot[index_B].e2 = rot_B.e2;
            matrix_rot[index_B].e3 = rot_B.e3;


            //--- > Rotate magnetization
        }
    }
}

inline double getJKRContactRadius(double compression_length, double r0, double R)
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

    /*    int i = 100;
        do
        {
            x_pow3 = x_old * x_old * x_old;
            x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius);

            if (std::abs(x_new - x_old) / particle_radius < 0.0001)
                break;

            x_old = x_new;
        } while (--i > 0);*/


    return x_new * x_new;
}




inline void updateParticleInteraction(vec3D* pos_new, vec3D* force_new, vec3D* torque_old, vec3D* torque_new, vec3D* dMdt_new, vec3D* matrix_con, 
    vec3D* matrix_norm, vec3D* omega, vec3D* omega_tot, vec3D* mag, quat* matrix_rot, double* matrix_comp, double* matrix_twist, double* amon, double* moment, material* mat,  int* matIDs, vec3D B_ext, int Nmon, double time_step)
{
    // init forces & torques with 0
    memset(force_new, 0, Nmon * sizeof(vec3D));
    memset(torque_new, 0, Nmon * sizeof(vec3D));
    memset(dMdt_new, 0, Nmon * sizeof(vec3D));

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

            vec3D n_c = vec3D_diff(pos_A, pos_B);
            double particle_distance = vec3D_length(n_c);
            vec3D_normalize(n_c);

            vec3D force_tmp;

            force_tmp.x = 1e-12 * n_c.x;
            force_tmp.y = 1e-12 * n_c.y;
            force_tmp.z = 1e-12 * n_c.z;

            force_new[i].x += force_tmp.x;
            force_new[i].y += force_tmp.y;
            force_new[i].z += force_tmp.z;

            //calculate magnetization
            /*double chi_A = mat[mat_id_A].chi;

            if (abs(chi_A) > 0)
            {
                double len_B = vec3D_length(B_ext);
                double Vmon_A = 4. / 3. * PI * a_mon_A * a_mon_A * a_mon_A;
                double Vmon_B = 4. / 3. * PI * a_mon_B * a_mon_B * a_mon_B;

                double Msat_A = mat[mat_id_A].Msat;

                double tau_ss_A = mat[mat_id_A].tss;
                double tau_sl_A = mat[mat_id_A].tsl;

                vec3D mag_A = mag[i];
                vec3D mag_B = mag[j];

                vec3D omega_tot_A = omega_tot[i];

                vec3D Mfin, Mpara, Mperp;

                vec3D_set(Mfin, 0);
                vec3D_set(Mpara, 0);
                vec3D_set(Mperp, 0);

                double len_mag_A = vec3D_length(mag_A);
                double len_mag_B = vec3D_length(mag_B);

                if (abs(chi_A) > LIMIT_FER) //ferromagnetic material
                {
                    double len_omega = vec3D_length(omega_tot_A);

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

                    double len_Mfin = vec3D_length(Mfin);
                    double dot_magA_Mfin = vec3D_dot(mag_A, Mfin);

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

                    double len_Mfin = vec3D_length(Mfin);

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

                    /*double dot_magA_Mfin = vec3D_dot(mag_A, Mfin);

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
                    vec3D cross = vec3D_cross(mag_A, B_ext);

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
                    vec3D torque_B = vec3D_cross(mu_A, B_ext);

                    vec3D vec_d = vec3D_diff(pos_B, pos_A); //todo: check for correct direction
                    double d = vec3D_length(vec_d);
                    double d2 = d * d;
                    double d3 = d2 * d;
                    double d4 = d2 * d2;
                    double d5 = d2 * d3;

                    //double fD = (3. * mu0) / (PIx4 * d5);
                    double fD = (3. * mu0) / (PIx4 * d4);
                    double tD = (mu0 / PIx4);
                    double tEdH = (Vmon_A / GAMMA_E);

                    double dot_muA_d = vec3D_dot(mu_A, vec_d);
                    double dot_muB_d = vec3D_dot(mu_B, vec_d);
                    double dot_muA_muB = vec3D_dot(mu_A, mu_B);

                    vec3D cross_d_muB = vec3D_cross(vec_d, mu_B);
                    vec3D cross_muA_muB = vec3D_cross(mu_A, mu_B);

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
            // end of magnetization
            

            //calc. surface forces only for con. monomers
            if (matrix_comp[i * Nmon + j] == -1.)
                continue;

            bool update_contact_pointers = false;

            int index_A = 0 * Nmon * Nmon + i * Nmon + j;
            int index_B = 1 * Nmon * Nmon + i * Nmon + j;

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

            // calculate current contact pointers
            vec3D n_A = matrix_con[index_A];
            vec3D n_B = matrix_con[index_B];
            vec3D delta_n = vec3D_diff(n_A, n_B);

            // determine distance between particles & contact normal


            // -> elastic force
            double compression_length = a_mon_A + a_mon_B - particle_distance;
            double contact_radius = getJKRContactRadius(compression_length, a0, R);
            double Fc = 3 * PI * gamma * R;
            double force_elastic = 4.0 * Fc * (pow(contact_radius / a0, 3.0) - pow(contact_radius / a0, 3.0 / 2.0));

            force_new[i].x += force_elastic * n_c.x;
            force_new[i].y += force_elastic * n_c.y;
            force_new[i].z += force_elastic * n_c.z;

            // -> damping force
            double old_compression_length = matrix_comp[i * Nmon + j];
            double vis_damp_const = 2.0 * T_vis / (nu_A * nu_B) * Es;
            double delta_dot = (compression_length - old_compression_length) / time_step;
            double force_damp = vis_damp_const * contact_radius * delta_dot;

            force_new[i].x += force_damp * n_c.x;
            force_new[i].y += force_damp * n_c.y;
            force_new[i].z += force_damp * n_c.z;

            // -> sliding
            double dot = vec3D_dot(delta_n, n_c);

            //is a_mon_A correct here?
            vec3D displacement;
            displacement.x = a_mon_A * (delta_n.x - dot * n_c.x);
            displacement.y = a_mon_A * (delta_n.y - dot * n_c.y);
            displacement.z = a_mon_A * (delta_n.z - dot * n_c.z);

            double displacement_norm = vec3D_length(displacement);

            double crit_sliding_displacement_modifier = 1.0;
            double crit_sliding_displacement = crit_sliding_displacement_modifier * (2.0 - nu_A) / (16.0 * PI) * a0;

            // check if we are in the inelastic regime
            if (displacement_norm > crit_sliding_displacement)
            {
                // determine correction of contact pointers
                vec3D displacement_correction;
                displacement_correction.x = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.x;
                displacement_correction.y = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.y;
                displacement_correction.z = (1.0 - crit_sliding_displacement / displacement_norm) * displacement.z;

                // calculate correction factor (see Wada et al. 2007 appendix for details)
                double inv_norm = 1.0 / vec3D_length_sq(displacement_correction);

                dot = vec3D_dot(n_A, displacement_correction);
                double alpha_A = 1.0 / (1.0 - dot * dot * inv_norm);

                dot = vec3D_dot(n_A, displacement_correction);
                double alpha_B = 1.0 / (1.0 - dot * dot * inv_norm);

                // correct contact pointers
                double particle_radius_inv = 1.0 / a_mon_A;
                n_A.x -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.x;
                n_A.y -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.y;
                n_A.z -= 0.5 * particle_radius_inv * alpha_A * displacement_correction.z;

                n_B.x += 0.5 * particle_radius_inv * alpha_B * displacement_correction.x;
                n_B.y += 0.5 * particle_radius_inv * alpha_B * displacement_correction.y;
                n_B.z += 0.5 * particle_radius_inv * alpha_B * displacement_correction.z;

                vec3D_normalize(n_A);
                vec3D_normalize(n_B);

                update_contact_pointers = true;
            }

            // -> rolling
            double crit_rolling_displacement = 0.5 * (mat[mat_id_A].xi + mat[mat_id_B].xi);

            displacement.x = R * (n_A.x + n_B.x);
            displacement.y = R * (n_A.y + n_B.y);
            displacement.z = R * (n_A.z + n_B.z);
            displacement_norm = vec3D_length(displacement);

            if (displacement_norm > crit_rolling_displacement)
            {
                // determine correction of contact pointers
                vec3D displacement_correction;
                displacement_correction.x = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.x;
                displacement_correction.y = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.y;
                displacement_correction.z = (1.0 - crit_rolling_displacement / displacement_norm) * displacement.z;

                // calculate correction factor (see Wada et al. _B007 appendix for details)
                double inv_norm = 1.0 / vec3D_length_sq(displacement_correction);

                dot = vec3D_dot(n_A, displacement_correction);
                double alpha_A = 1.0 / (1.0 - dot * dot * inv_norm);

                dot = vec3D_dot(n_B, displacement_correction);
                double alpha_B = 1.0 / (1.0 - dot * dot * inv_norm);

                // correct contact pointers
                double particle_radius_inv = 1.0 / a_mon_A;
                n_A.x -= particle_radius_inv * alpha_A * displacement_correction.x;
                n_A.y -= particle_radius_inv * alpha_A * displacement_correction.y;
                n_A.z -= particle_radius_inv * alpha_A * displacement_correction.z;

                n_B.x -= particle_radius_inv * alpha_B * displacement_correction.x;
                n_B.y -= particle_radius_inv * alpha_B * displacement_correction.y;
                n_B.z -= particle_radius_inv * alpha_B * displacement_correction.z;

                vec3D_normalize(n_A);
                vec3D_normalize(n_B);

                update_contact_pointers = true;
            }

            if (update_contact_pointers)
            {
                updateNormal(n_A, n_B, matrix_con, matrix_rot, i, j, Nmon);
                //delta_n = vec3D_diff(n_A, n_B);
            }

            //sliding force
            double sliding_modifier = 1.0;
            double k_s = sliding_modifier * 8.0 * Gs * a0;
            vec3D displacement_zeta;
            double tmp_A = a_mon_A * vec3D_dot(n_A, n_c);
            double tmp_B = a_mon_B * vec3D_dot(n_B, n_c);

            displacement_zeta.x = a_mon_A * n_A.x - a_mon_B * n_B.x - (tmp_A - tmp_B) * n_c.x;
            displacement_zeta.y = a_mon_A * n_A.y - a_mon_B * n_B.y - (tmp_A - tmp_B) * n_c.y;
            displacement_zeta.z = a_mon_A * n_A.z - a_mon_B * n_B.z - (tmp_A - tmp_B) * n_c.z;

            //-> sliding force sliding torque
            vec3D tmp; //helper variable

            tmp.x = a_mon_B * n_B.x - a_mon_A * n_A.x;
            tmp.y = a_mon_B * n_B.y - a_mon_A * n_A.y;
            tmp.z = a_mon_B * n_B.z - a_mon_A * n_A.z;

            double force_sliding = -k_s * vec3D_dot(displacement_zeta, tmp) / particle_distance;

            if (abs(force_sliding) > 1.0e-10)
                force_sliding = 1e-10;

            //cout << force_sliding << endl << flush;
            force_new[i].x += force_sliding * n_c.x;
            force_new[i].y += force_sliding * n_c.y;
            force_new[i].z += force_sliding * n_c.z;

            vec3D torque_sliding;
            tmp = vec3D_cross(n_A, displacement_zeta);

            torque_sliding.x = -a_mon_A * k_s * tmp.x;
            torque_sliding.y = -a_mon_A * k_s * tmp.y;
            torque_sliding.z = -a_mon_A * k_s * tmp.z;

            torque_new[i].x += torque_sliding.x;
            torque_new[i].y += torque_sliding.y;
            torque_new[i].z += torque_sliding.z;

            //-> rolling torque
            vec3D displacement_xi;
            vec3D torque_rolling;
            double rolling_modifier = 1.0;
            double k_r = rolling_modifier * 4.0 * Fc / R;

            displacement_xi.x = R * (n_A.x + n_B.x);
            displacement_xi.y = R * (n_A.y + n_B.y);
            displacement_xi.z = R * (n_A.z + n_B.z);

            tmp = vec3D_cross(n_A, displacement_xi);

            torque_rolling.x = -k_r * R * tmp.x;
            torque_rolling.y = -k_r * R * tmp.y;
            torque_rolling.z = -k_r * R * tmp.z;

            torque_new[i].x += torque_rolling.x;
            torque_new[i].y += torque_rolling.y;
            torque_new[i].z += torque_rolling.z;

            //-> twisting torque
            vec3D delta_omega_old, delta_omega_new, twisting_torque;
            double crit_twisting_displacement = 1.0 / (16.0 * PI);
            double twisting_modifier = 1.0;
            double k_t = twisting_modifier * 16.0 / 3.0 * Gs * a0 * a0 * a0;
            double twisting_displacement = matrix_twist[i * Nmon + j];
            double moment_inv_A = 1.0 / moment_A;
            double moment_inv_B = 1.0 / moment_B;

            //store old normal vector
            vec3D n_c_old = matrix_norm[i * Nmon + j];

            delta_omega_old.x = omega_A.x - omega_B.x;
            delta_omega_old.y = omega_A.y - omega_B.y;
            delta_omega_old.z = omega_A.z - omega_B.z;

            // update twisting displacement - use second order integration: omega^n+1 = omega^n + 0.5 (<delta_omega^n, n_c^n> + <delta_omega^n+1, n_c^n+1>)
            delta_omega_new.x = delta_omega_old.x + time_step * (moment_inv_A * torque_old[i].x - moment_inv_B * torque_old[j].x);
            delta_omega_new.y = delta_omega_old.y + time_step * (moment_inv_A * torque_old[i].y - moment_inv_B * torque_old[j].y);
            delta_omega_new.z = delta_omega_old.z + time_step * (moment_inv_A * torque_old[i].z - moment_inv_B * torque_old[j].z);

            twisting_displacement += 0.5 * time_step * (vec3D_dot(delta_omega_old, n_c_old) + vec3D_dot(delta_omega_new, n_c));

            if (twisting_displacement > crit_twisting_displacement)
                twisting_displacement = crit_twisting_displacement;

            twisting_torque.x = k_t * twisting_displacement * n_c.x;
            twisting_torque.y = k_t * twisting_displacement * n_c.y;
            twisting_torque.z = k_t * twisting_displacement * n_c.z;

            torque_new[i].x -= twisting_torque.x;
            torque_new[i].y -= twisting_torque.y;
            torque_new[i].z -= twisting_torque.z;/**/

            // store current state for next update step
            matrix_norm[i * Nmon + j].x = n_c.x;
            matrix_norm[i * Nmon + j].y = n_c.y;
            matrix_norm[i * Nmon + j].z = n_c.z;
            matrix_comp[i * Nmon + j] = compression_length;
        }
    }
}

#endif