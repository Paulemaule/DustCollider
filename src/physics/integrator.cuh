#pragma once

#include "utils/vector.cuh"
#include "utils/typedefs.cuh"
#include "utils/constant.cuh"

/**
 * The predictor step of the PECE scheme.
 */
__global__ void gpu_predictor(
    const double3* pos_curr,
    double3* pos_next,
    const double3* force_curr,
    const double3* vel_curr,
    const double* mass,
    const double timestep,
    const int Nmon
) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    if (threadID < Nmon) return;

    double inv_mass = 1.0 / mass[threadID];

    pos_next[threadID].x = pos_curr[threadID].x + timestep * vel_curr[threadID].x + 0.5 * inv_mass * timestep * timestep * force_curr[threadID].x;
    pos_next[threadID].y = pos_curr[threadID].y + timestep * vel_curr[threadID].y + 0.5 * inv_mass * timestep * timestep * force_curr[threadID].y;
    pos_next[threadID].z = pos_curr[threadID].z + timestep * vel_curr[threadID].z + 0.5 * inv_mass * timestep * timestep * force_curr[threadID].z;
}

/**
 * The corrector step of the PECE scheme.
 */
__global__ void gpu_corrector(
    double3* force_curr,
    double3* force_next,
    double3* torque_curr,
    double3* torque_next,
    double3* dMdt_curr,
    double3* dMdt_next,
    double3* vel_curr,
    double3* omega_curr,
    double3* omega_tot_curr,
    double3* mag_curr,
    double* mass,
    double* moment,
    material* mat,
    int* matIDs,
    double timestep,
    int Nmon
) {
    // TODO: I should attempt to reduce global memory accesses as much as possible...
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    if (threadID < Nmon) return;

    double inv_mass = 1.0 / mass[threadID];
    double inv_moment = 1.0 / moment[threadID];

    // Calculate acceleration
    double3 acc;
    acc.x = 0.5 * inv_mass * (force_curr[threadID].x + force_next[threadID].x);
    acc.y = 0.5 * inv_mass * (force_curr[threadID].y + force_next[threadID].y);
    acc.z = 0.5 * inv_mass * (force_curr[threadID].z + force_next[threadID].z);

    // Update velocity
    vel_curr[threadID].x += timestep * acc.x;
    vel_curr[threadID].y += timestep * acc.y;
    vel_curr[threadID].z += timestep * acc.z;

    // Update angular momentum
    omega_curr[threadID].x += 0.5 * inv_moment * timestep * (torque_curr[threadID].x + torque_next[threadID].x);
    omega_curr[threadID].y += 0.5 * inv_moment * timestep * (torque_curr[threadID].y + torque_next[threadID].y);
    omega_curr[threadID].z += 0.5 * inv_moment * timestep * (torque_curr[threadID].z + torque_next[threadID].z);

    // Calculate trajectories curvature
    double vel_sq = vec_lenght_sq(vel_curr[threadID]);
    double3 cross = vec_cross(vel_curr[threadID], acc);

    if (vel_sq > 0) {
        omega_tot_curr[threadID].x = omega_curr[threadID].x + cross.x / vel_sq;
        omega_tot_curr[threadID].y = omega_curr[threadID].y + cross.y / vel_sq;
        omega_tot_curr[threadID].z = omega_curr[threadID].z + cross.z / vel_sq;
    } else {
        omega_tot_curr[threadID].x = omega_curr[threadID].x;
        omega_tot_curr[threadID].y = omega_curr[threadID].y;
        omega_tot_curr[threadID].z = omega_curr[threadID].z;
    }

    // Update magnetization
    int mat_id = matIDs[threadID];
    double chi = mat[mat_id].chi;
    double Msat = mat[mat_id].Msat;

    // TODO: Check if this leads to significant warp divergence
    if (abs(chi) > 0) {
        mag_curr[threadID].x += 0.5 * timestep * (dMdt_curr[threadID].x + dMdt_next[threadID].x);
        mag_curr[threadID].y += 0.5 * timestep * (dMdt_curr[threadID].y + dMdt_next[threadID].y);
        mag_curr[threadID].z += 0.5 * timestep * (dMdt_curr[threadID].z + dMdt_next[threadID].z);
    
        double len_mag = vec_lenght(mag_curr[threadID]);

        if (len_mag > 0) {
            if (abs(chi) > LIMIT_FER) {
                mag_curr[threadID].x = Msat * mag_curr[threadID].x / len_mag;
                mag_curr[threadID].y = Msat * mag_curr[threadID].y / len_mag;
                mag_curr[threadID].z = Msat * mag_curr[threadID].z / len_mag;
            } else {
                if (len_mag > Msat) {
                    mag_curr[threadID].x = Msat * mag_curr[threadID].x / len_mag;
                    mag_curr[threadID].y = Msat * mag_curr[threadID].y / len_mag;
                    mag_curr[threadID].z = Msat * mag_curr[threadID].z / len_mag;
                }
            }
        }  
    }
}