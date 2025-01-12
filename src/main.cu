#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <string>
//#include <studio.h>
#include <iostream>
#include<fstream>
#include <algorithm>

using namespace std;
using namespace std::chrono;

#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "typedefs.cuh"
#include "vector.cuh"
#include "pipeline.cuh"
#include "physics.cuh"

int main(const int argc, const char** argv)
{
    CPipeline pipeline;

    auto start = high_resolution_clock::now();

    // Check command line input for command file.
    if (!pipeline.init(argc, argv))
        return -1;

    // Read the command file.
    if (!pipeline.parse())
        return -1;

    // Check run parametes and create output directories.
    if (!pipeline.checkParameters())
        return -1;

    // FIXME: Change nullpointers to 'nullptr' and check for dangling pointers.
    material* mat = 0;          // An array containing the material parameters.
    material* buff_mat = 0;     // An array containing the material parameters in GPU memory.

    int Nmat = 0;               // The number of materials.

    // Reading and preparing the material parameters.
    pipeline.prepareMaterial(mat, Nmat);

    cudaMalloc(&buff_mat, Nmat * sizeof(material));

    // Definition of variables.
    vec3D* vel = 0;                 // An array of the monomer velocities.
    // FIXME: This is unphysical ?
    vec3D* omega = 0;               // An array of the monomer angular velocities due to self rotation.
    vec3D* omega_tot = 0;           // An array of the monomer angular velocities due to self rotation and curved trajectories.
    vec3D* mag = 0;                 // An array of the monomer magnetizations.

    vec3D* pos_old = 0;
    vec3D* pos_new = 0;
    vec3D* force_old = 0;
    vec3D* force_new = 0;
    vec3D* torque_old = 0;
    vec3D* torque_new = 0;

    vec3D* dMdt_old = 0;
    vec3D* dMdt_new = 0;

    vec3D* buff_vel = 0;
    vec3D* buff_omega = 0;
    vec3D* buff_omega_tot = 0;
    vec3D* buff_mag = 0;

    vec3D* buff_pos_old = 0;
    vec3D* buff_pos_new = 0;
    vec3D* buff_force_old = 0;
    vec3D* buff_force_new = 0;
    vec3D* buff_torque_old = 0;
    vec3D* buff_torque_new = 0;

    vec3D* buff_dMdt_old = 0;
    vec3D* buff_dMdt_new = 0;

    vec3D* storage_pos = 0;         // Array of monomer position that are to be stored after simulation.
    vec3D* storage_vel = 0;         // Array of monomer velocities that are to be stored after simulation.
    vec3D* storage_force = 0;       // Array of monomer forces that are to be stored after simulation.
    vec3D* storage_torque = 0;      // Array of monomer torques that are to be stored after simulation.
    vec3D* storage_omega = 0;       // Array of monomer angular velocities that are to be stored after simulation.
    vec3D* storage_mag = 0;         // Array of monomer magnetizations that are to be stored after simulation.
    int* storage_cluster = 0;       // Array of monomer (// TODO: ?) that are to be stored after simulation.

    vec3D* matrix_con = 0;          // contact pointer between monomers
    vec3D* matrix_norm = 0;         // normal vectors between monomers
    quat* matrix_rot = 0;           // contact pointer rotation direction
    double* matrix_comp = 0;        // old compression lengths, also used to track connection
    double* matrix_twist = 0;       // old twisting displacement

    vec3D* buff_matrix_con = 0;     // GPU contact pointer
    vec3D* buff_matrix_norm = 0;    // GPU normal vectors
    quat* buff_matrix_rot = 0;      // GPU contact pointer rotation
    double* buff_matrix_comp = 0;   // GPU old compression lenghts
    double* buff_matrix_twist = 0;  // GPU old twisting displacement

    int* matIDs = 0;                // material IDs
    double* amon = 0;               // Monomer radii
    double* moment = 0;             // Monomer moments of inertia
    double* mass = 0;               // Monomer masses
    int* clusterIDs = 0;            // Monomer cluster membership

    int* buff_matIDs = 0;           // GPU material IDs
    double* buff_amon = 0;          // GPU monomer radii
    double* buff_moment = 0;        // GPU monomer moments of inertia
    double* buff_mass = 0;          // GPU monomer masses
    
    int Nmon = 0;                   // Number of monomers

    // Read the aggregate files and initialize the state.
    pipeline.prepareData(pos_old, vel, omega_tot, mag, amon, mass, moment, matIDs, Nmon);

    // Print run summary.
    pipeline.printParameters();

    // Initialize the vectors.
    omega = new vec3D[Nmon];
    pos_new = new vec3D[Nmon];
    force_old = new vec3D[Nmon];
    force_new = new vec3D[Nmon];
    torque_old = new vec3D[Nmon];
    torque_new = new vec3D[Nmon];
    dMdt_old = new vec3D[Nmon];
    dMdt_new = new vec3D[Nmon];

    // Sets the values of the vectors in CPU memory to 0
    memset(omega, 0, Nmon * sizeof(vec3D));
    memset(pos_new, 0, Nmon * sizeof(vec3D));
    memset(force_old, 0, Nmon * sizeof(vec3D));
    memset(force_new, 0, Nmon * sizeof(vec3D));
    memset(torque_old, 0, Nmon * sizeof(vec3D));
    memset(torque_new, 0, Nmon * sizeof(vec3D));
    memset(dMdt_old, 0, Nmon * sizeof(vec3D));
    memset(dMdt_new, 0, Nmon * sizeof(vec3D));

    // Allocate memory in the GPU
    cudaMalloc(&buff_vel, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_omega, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_omega_tot, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_mag, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_pos_old, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_pos_new, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_force_old, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_force_new, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_torque_old, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_torque_new, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_dMdt_old, Nmon*sizeof(vec3D));
    cudaMalloc(&buff_dMdt_new, Nmon*sizeof(vec3D));

    int Ncon = Nmon * Nmon; // Square of the number of monomers Nmon^2.

    // Initialize the matrices.
    matrix_con = new vec3D[2 * Ncon];
    matrix_rot = new quat[2 * Ncon];
    matrix_norm = new vec3D[Ncon];
    matrix_comp = new double[Ncon];
    matrix_twist = new double[Ncon];

    memset(matrix_con, 0, 2 * Ncon * sizeof(vec3D));
    memset(matrix_rot, 0, 2 * Ncon * sizeof(quat));
    memset(matrix_norm, 0, Ncon * sizeof(vec3D));
    fill(matrix_comp, matrix_comp + Ncon, -1.0);
    memset(matrix_twist, 0, Ncon * sizeof(double));

    cudaMalloc(&buff_matrix_con, Ncon * sizeof(vec3D));
    cudaMalloc(&buff_matrix_rot, Ncon * sizeof(quat));
    cudaMalloc(&buff_matrix_norm, Ncon * sizeof(vec3D));
    cudaMalloc(&buff_matrix_comp, Ncon * sizeof(double));
    cudaMalloc(&buff_matrix_twist, Ncon * sizeof(double));
    cudaMalloc(&buff_matIDs, Nmon * sizeof(int));
    cudaMalloc(&buff_amon, Nmon * sizeof(double));
    cudaMalloc(&buff_moment, Nmon * sizeof(double));
    cudaMalloc(&buff_mass, Nmon * sizeof(double));

    // Get some run parameters.
    ullong N_iter = pipeline.getNIter();
    ullong N_save = pipeline.getNSave();
    double time_step = pipeline.getTimeStep();
    bool save_ovito = pipeline.saveOvito();
    vec3D B_ext = pipeline.getBext();

    int N_store = 0; // The number of individual monomer states that need to be stored.

    if (N_save > 0)
    {
        N_store = Nmon * int((double(N_iter) / double(N_save) + 0.5));

        clusterIDs = new int[Nmon];
        fill(clusterIDs, clusterIDs + Nmon, -1);

        if (pipeline.savePos())
        {
            storage_pos = new vec3D[N_store];
            memset(storage_pos, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveVel())
        {
            storage_vel = new vec3D[N_store];
            memset(storage_vel, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveForce())
        {
            storage_force = new vec3D[N_store];
            memset(storage_force, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveTorque())
        {
            storage_torque = new vec3D[N_store];
            memset(storage_torque, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveOmega())
        {
            storage_omega = new vec3D[N_store];
            memset(storage_omega, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveMag())
        {
            storage_mag = new vec3D[N_store];
            memset(storage_mag, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveCluster())
        {
            storage_cluster = new int[N_store];
            memset(storage_cluster, 0, N_store * sizeof(int));
        }
    }

    // Calculate the number of blocks
    int nBlocks = (Nmon + BLOCK_SIZE + 1) / BLOCK_SIZE; // FIXME: This should be Nmon + BLOCK_SIZE - 1

    // Copy memory into GPU
    cudaMemcpy(buff_mat, mat, Nmat * sizeof(material), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_pos_old, pos_old, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_pos_new, pos_new, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_force_old, force_old, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_force_new, force_new, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_torque_old, torque_old, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_torque_new, torque_new, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_dMdt_old, dMdt_old, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_dMdt_new, dMdt_new, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_vel, vel, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_omega, omega, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_omega_tot, omega_tot, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_mag, mag, Nmon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_matIDs, matIDs, Nmon * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_amon, amon, Nmon * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_moment, moment, Nmon * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_mass, mass, Nmon * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_matrix_con, matrix_con, Ncon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_matrix_rot, matrix_rot, Ncon * sizeof(quat), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_matrix_norm, matrix_norm, Ncon * sizeof(vec3D), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_matrix_comp, matrix_comp, Ncon * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(buff_matrix_twist, matrix_twist, Ncon * sizeof(double), cudaMemcpyHostToDevice);

    ullong counter_save = 0;

    // THE MAIN SIMULATION LOOP
    for (ullong iter = 0; iter < N_iter; iter++) // The simulation iteration count
    {
        #ifdef RUN_ON_GPU

        gpu_predictor << < nBlocks, BLOCK_SIZE >> > (buff_pos_old, buff_pos_new, buff_force_old, buff_vel, buff_mass, time_step, Nmon);
        cudaDeviceSynchronize();

        gpu_updateNeighbourhoodRelations << < nBlocks, BLOCK_SIZE >> >  (buff_pos_new, buff_matrix_con, buff_matrix_norm, buff_matrix_rot,
                                buff_matrix_comp, buff_matrix_twist, buff_amon, buff_mat, buff_matIDs, Nmon);
        cudaDeviceSynchronize();

        gpu_updateContacts << < nBlocks, BLOCK_SIZE >> >  (buff_omega, buff_omega_tot, buff_torque_old, buff_mag, buff_matrix_rot, buff_matrix_comp, buff_moment, Nmon, time_step);
        cudaDeviceSynchronize();

        gpu_updateParticleInteraction << < nBlocks, BLOCK_SIZE >> > (buff_pos_new, buff_force_new, buff_torque_old, buff_torque_new, buff_dMdt_new, 
                                buff_matrix_con, buff_matrix_norm, buff_omega, buff_omega_tot, buff_mag, buff_matrix_rot, buff_matrix_comp, buff_matrix_twist, 
                                buff_amon, buff_moment, buff_mat, buff_matIDs, B_ext, Nmon, time_step);
        cudaDeviceSynchronize();

        gpu_corrector << < nBlocks, BLOCK_SIZE >> > (buff_force_old, buff_force_new, buff_torque_old, buff_torque_new, buff_dMdt_old, 
                                buff_dMdt_new, buff_vel, buff_omega, buff_omega_tot, buff_mag, buff_mass, buff_moment, buff_mat, 
                                buff_matIDs, time_step, Nmon);
        cudaDeviceSynchronize();

        switch_pointer(buff_pos_old, buff_pos_new, buff_force_old, buff_force_new, buff_torque_old, buff_torque_new, buff_dMdt_old, buff_dMdt_new);
        
        #else
        cpu_predictor(pos_old, pos_new, force_old, vel, mass, time_step, Nmon);
        cpu_updateNeighbourhoodRelations(pos_new, matrix_con, matrix_norm, matrix_rot, matrix_comp, matrix_twist, amon, mat, matIDs, Nmon);
        cpu_updateContacts(omega, omega_tot, torque_old, mag, matrix_rot, matrix_comp, moment, Nmon, time_step);

        cpu_updateParticleInteraction(pos_new, force_new, torque_old, torque_new, dMdt_new, matrix_con, matrix_norm, omega, omega_tot, mag, matrix_rot, matrix_comp, matrix_twist, amon, moment, mat, matIDs, B_ext, Nmon, time_step);
        cpu_corrector(force_old, force_new, torque_old, torque_new, dMdt_old, dMdt_new, vel, omega, omega_tot, mag, mass, moment, mat, matIDs, time_step, Nmon);

        switch_pointer(pos_old, pos_new, force_old, force_new, torque_old, torque_new, dMdt_old, dMdt_new);
        #endif

        // 
        if (N_save > 0)
        {
            if (iter % N_save == 0)
            {
                ullong start_index = counter_save * Nmon;

                if (long(N_store) - long(start_index) < Nmon)
                {
                    cout << "ERROR: 'Storage overrun!  \n";
                    cout << long(N_store) - long(start_index) << "\t" << Nmon << endl;
                    return -1;
                }

                if (start_index < N_store) //just a failsave if there was a rounding error in N_store
                {
                    if (storage_pos != 0)
                    {
                        #ifdef RUN_ON_GPU
                        cudaMemcpy(pos_old, buff_pos_old, Nmon * sizeof(vec3D), cudaMemcpyDeviceToHost);
                        #endif
                        copy(pos_old, pos_old + Nmon, storage_pos + start_index);
                    }

                    if (storage_vel != 0)
                    {
                        #ifdef RUN_ON_GPU
                        cudaMemcpy(vel, buff_vel, Nmon * sizeof(vec3D), cudaMemcpyDeviceToHost);
                        #endif
                        copy(vel, vel + Nmon, storage_vel + start_index);
                    }

                    if (storage_force != 0)
                    {
                        #ifdef RUN_ON_GPU
                        cudaMemcpy(force_old, buff_force_old, Nmon * sizeof(vec3D), cudaMemcpyDeviceToHost);
                        #endif
                        copy(force_old, force_old + Nmon, storage_force + start_index);
                    }

                    if (storage_torque != 0)
                    {
                        #ifdef RUN_ON_GPU
                        cudaMemcpy(torque_old, buff_torque_old, Nmon * sizeof(vec3D), cudaMemcpyDeviceToHost);
                        #endif
                        copy(torque_old, torque_old + Nmon, storage_torque + start_index);
                    }

                    if (storage_omega != 0)
                    {
                        #ifdef RUN_ON_GPU
                        cudaMemcpy(omega_tot, buff_omega_tot, Nmon * sizeof(vec3D), cudaMemcpyDeviceToHost);
                        #endif
                        copy(omega_tot, omega_tot + Nmon, storage_omega + start_index);
                    }

                    if (storage_mag != 0)
                    {
                        #ifdef RUN_ON_GPU
                        cudaMemcpy(mag, buff_mag, Nmon * sizeof(vec3D), cudaMemcpyDeviceToHost);
                        #endif
                        copy(mag, mag + Nmon, storage_mag + start_index);
                    }

                    if (storage_cluster != 0)
                    {
                        fill(clusterIDs, clusterIDs + Nmon, -1);
                        cpu_findConnectedComponents(Nmon, matrix_comp, clusterIDs);
                        copy(clusterIDs, clusterIDs + Nmon, storage_cluster + start_index);
                    }
                }
                counter_save++;
            }
        }

        if (iter % 2000 == 0)
        {
            printf("-> Simulation progress: %.4f %%      \r", 100.0 * float(iter) / N_iter);
        }
    }

    cout << CLR_LINE;

    if (N_save > 0)
    {
        cout << SEP_LINE;

        cout << "Writing simulation data:    \n" << flush;

        // Write ouput files
        if (save_ovito && storage_pos != 0)
        {
            if (!pipeline.writeAllOVITO(storage_pos, storage_vel, storage_force, storage_torque, storage_omega, storage_mag, storage_cluster, amon, matIDs, Nmon, N_store))
                return -1;
        }

        if (!pipeline.writeHeader())
            return -1;

        if (!pipeline.writeBinaryDouble("agg_a_mon.bin", amon, Nmon))
            return -1;

        if (!pipeline.writeBinaryDouble("agg_mass_mon.bin", mass, Nmon))
            return -1;

        for (int i = 0; i < Nmon; i++)
            matIDs[i] = matIDs[i] + 1;

        if (!pipeline.writeBinaryInt("agg_matid_mon.bin", matIDs, Nmon))
            return -1;

        if (storage_pos != 0)
        {
            if (!pipeline.writeBinaryVec("sim_pos.bin", storage_pos, N_store))
                return -1;
        }

        if (storage_vel != 0)
            if (!pipeline.writeBinaryVec("sim_vel.bin", storage_vel, N_store))
                return -1;

        if (storage_force != 0)
            if (!pipeline.writeBinaryVec("sim_force.bin", storage_force, N_store))
                return -1;

        if (storage_torque != 0)
            if (!pipeline.writeBinaryVec("sim_torque.bin", storage_torque, N_store))
                return -1;

        if (storage_cluster != 0)
            if (!pipeline.writeBinaryInt("sim_cluster.bin", storage_cluster, N_store))
                return -1;
    }

    cout << "-> Final cleanup ...              \r" << flush;

    // Free allocated memory
    if (omega != 0)
        delete[] omega;

    if (omega_tot != 0)
        delete[] omega_tot;

    if (vel != 0)
        delete[] vel;

    if (mag != 0)
        delete[] mag;

    if (pos_old != 0)
        delete[] pos_old;

    if (pos_new != 0)
        delete[] pos_new;

    if (matIDs != 0)
        delete[] matIDs;

    if (amon != 0)
        delete[] amon;

    if (moment != 0)
        delete[] moment;

    if (matrix_con != 0)
        delete[] matrix_con;

    if (matrix_norm != 0)
        delete[] matrix_norm;

    if (matrix_rot != 0)
        delete[] matrix_rot;

    if (matrix_comp != 0)
        delete[] matrix_comp;

    if (matrix_twist != 0)
        delete[] matrix_twist;

    if (mass != 0)
        delete[] mass;

    if (force_old != 0)
        delete[] force_old;

    if (force_new != 0)
        delete[] force_new;

    if (torque_old != 0)
        delete[] torque_old;

    if (torque_new != 0)
        delete[] torque_new;

    if (dMdt_old != 0)
        delete[] dMdt_old;

    if (dMdt_new != 0)
        delete[] dMdt_new;

    if (mat != 0)
        delete[] mat;

    if (storage_pos != 0)
        delete[] storage_pos;

    if (storage_vel != 0)
        delete[] storage_vel;

    if (storage_force != 0)
        delete[] storage_force;

    if (storage_torque != 0)
        delete[] storage_torque;

    if (storage_omega != 0)
        delete[] storage_omega;

    if (storage_mag != 0)
        delete[] storage_mag;

    if (storage_cluster != 0)
        delete[] storage_cluster;

    cudaFree(buff_mat);

    cudaFree(buff_vel);
    cudaFree(buff_omega);
    cudaFree(buff_omega_tot);
    cudaFree(buff_mag);

    cudaFree(buff_pos_old);
    cudaFree(buff_pos_new);
    cudaFree(buff_force_old);
    cudaFree(buff_force_new);
    cudaFree(buff_torque_old);
    cudaFree(buff_torque_new);

    cudaFree(buff_dMdt_old);
    cudaFree(buff_dMdt_new);

    cudaFree(buff_matrix_con);
    cudaFree(buff_matrix_rot);
    cudaFree(buff_matrix_norm);
    cudaFree(buff_matrix_comp);
    cudaFree(buff_matrix_twist);

    cudaFree(buff_matIDs);
    cudaFree(buff_amon);
    cudaFree(buff_moment);
    cudaFree(buff_mass);

    cout << SEP_LINE;
    cout << "  - Final clanup: done              \n" << flush;
    cout << SEP_LINE;

    auto end = high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("Run time for %llu iterations : %.3f seconds.\n", N_iter, elapsed.count() * 1e-9);
    cout << SEP_LINE;

    return 0;
}