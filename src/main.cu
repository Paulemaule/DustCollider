// TODO: Add std:: , ... back for clarity
// TODO: Es macht wahrscheinlich Sinn die material eigenschaften nicht in einem seperaten array *mat sondern in arrays the form *pos zu speichern weil die Pattern chi_i = mat[matIDs[i]].chi wahrscheinlich zu race conditions führt.
// TODO: Einheitliche Variablennamen für monomer radius (r_i) und contact surface radius (a_ij)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <format>

#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "utils/config.cuh"
#include "utils/constant.cuh"
#include "utils/errors.cuh"
#include "utils/printing.cuh"
#include "utils/typedefs.cuh"

#include "pipeline.cuh"

#include "physics/state.cuh"
#include "physics/materialProperties.cuh"
#include "physics/integrator.cuh"

int main(const int argc, const char** argv)
{
    /* ############################################################
            INITIALIZATION
       ############################################################ */


    CPipeline pipeline;

    auto start = std::chrono::high_resolution_clock::now();

    // Check command line input for command file.
    if (!pipeline.init(argc, argv))
        return -1;

    PRINT_CLR_LINE();
    PRINT_TITLE("INITIALIZATION")
    PRINT_CLR_LINE();

    // Read the command file.
    if (!pipeline.parse())
        return -1;

    // Check run parametes and create output directories.
    if (!pipeline.checkParameters())
        return -1;

    // Reading material properties from the command file.
    material* mat = nullptr;          // An array containing the material parameters.
    int Nmat = 0;                     // The number of different materials.

    pipeline.prepareMaterial(mat, Nmat);

    // Declaring variables for the snapshots.
    double3* storage_pos = nullptr;         // Array of monomer position that are to be stored after simulation.
    double3* storage_vel = nullptr;         // Array of monomer velocities that are to be stored after simulation.
    double3* storage_force = nullptr;       // Array of monomer forces that are to be stored after simulation.
    double3* storage_torque = nullptr;      // Array of monomer torques that are to be stored after simulation.
    double3* storage_omega = nullptr;       // Array of monomer angular velocities that are to be stored after simulation.
    double3* storage_mag = nullptr;         // Array of monomer magnetizations that are to be stored after simulation.
    double* storage_potential_N = nullptr;  // Array of potential energies in normal direction.
    double* storage_potential_S = nullptr;  // Array of potential energies in sliding direction.
    double* storage_potential_R = nullptr;  // Array of potential energies in rolling direction.
    double* storage_potential_T = nullptr;  // Array of potential energies in twisting direction.
    double* storage_inelastic_N = nullptr;  // Array of dissipated energies in normal direction.
    double* storage_inelastic_S = nullptr;  // Array of dissipated energies in sliding direction.
    double* storage_inelastic_R = nullptr;  // Array of dissipated energies in rolling direction.
    double* storage_inelastic_T = nullptr;  // Array of dissipated energies in twisting direction.
    int* storage_cluster = nullptr;         // TODO
    int* clusterIDs = nullptr;              // TODO

    // Read the initial state of the system from the monomer files
    // Prepare arrays for the monomer data
    double3* initial_position = nullptr;
    double3* initial_velocity = nullptr;
    double3* initial_omega_tot = nullptr;
    double3* initial_magnetization = nullptr;
    double* monomer_radius = nullptr;
    double* monomer_mass = nullptr;
    double* monomer_moment = nullptr;
    int* monomer_matID = nullptr;

    int Nmon = 0; // The number of monomers.

    // Read aggregate files and prepare initial state.
    pipeline.prepareData(initial_position, initial_velocity, initial_omega_tot, initial_magnetization, 
                            monomer_radius, monomer_mass, monomer_moment, monomer_matID, Nmon);

    // Initialize the pointer containers
    // Allocating memory on the host
    hostState host_state_curr; // A container for pointers to the state of the monomers in the current timestep, in host memory.
    hostState host_state_next; // A container for pointers to the state of the monomers in the next timestep, in host memory.
    state_allocateHostMemory(host_state_curr, Nmon);
    state_allocateHostMemory(host_state_next, Nmon);

    hostMatProperties host_matProperties; // A container for pointers to the material properties of the monomers, in host memory.
    matProperties_allocateHostMemory(host_matProperties, Nmon);

    // Allocating memory on the device
    deviceState device_state_curr; // A container for pointers to the state of the monomers in the current timestep, in device memory.
    deviceState device_state_next; // A container for pointers to the state of the monomers in the next timestep, in device memory.
    state_allocateDeviceMemory(device_state_curr, Nmon);
    state_allocateDeviceMemory(device_state_next, Nmon);

    deviceMatProperties device_matProperties; // A container for pointers to the material properties of the monomers, in device memory.
    matProperties_allocateDeviceMemory(device_matProperties, Nmon);

    // Zeroing the arrays
    state_clearHost(host_state_curr, Nmon);
    state_clearHost(host_state_next, Nmon);
    matProperties_clearHost(host_matProperties, Nmon);

    // Copy initial values into state container
    CHECK_CUDA(cudaMemcpy(host_state_curr.position, initial_position, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost));
    CHECK_CUDA(cudaMemcpy(host_state_curr.velocity, initial_velocity, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost));
    CHECK_CUDA(cudaMemcpy(host_state_curr.omega, initial_omega_tot, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost));
    CHECK_CUDA(cudaMemcpy(host_state_curr.magnetization, initial_magnetization, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost));
    // TODO: This is no longer necessary, I think.
    // Initially all monomer pairs are marked as unconnected.
    for (int i = 0; i < Nmon * Nmon; i++) {
        host_state_curr.contact_compression[i] = -1.0;
    }

    // Push the initialized states to the device
    state_pushToDevice(host_state_curr, device_state_curr, Nmon);
    state_pushToDevice(host_state_next, device_state_next, Nmon);

    // Copy material values into properties container
    CHECK_CUDA(cudaMemcpy(host_matProperties.radius, monomer_radius, Nmon * sizeof(double), cudaMemcpyKind::cudaMemcpyHostToHost));
    CHECK_CUDA(cudaMemcpy(host_matProperties.mass, monomer_mass, Nmon * sizeof(double), cudaMemcpyKind::cudaMemcpyHostToHost));
    CHECK_CUDA(cudaMemcpy(host_matProperties.moment, monomer_moment, Nmon * sizeof(double), cudaMemcpyKind::cudaMemcpyHostToHost));
    CHECK_CUDA(cudaMemcpy(host_matProperties.matID, monomer_matID, Nmon * sizeof(int), cudaMemcpyKind::cudaMemcpyHostToHost));

    for (int i = 0; i < Nmon; i++) {
        host_matProperties.density[i] = mat[monomer_matID[i]].rho;
    }
    for (int i = 0; i < Nmon; i++) {
        host_matProperties.surface_energy[i] = mat[monomer_matID[i]].gamma;
    }
    for (int i = 0; i < Nmon; i++) {
        host_matProperties.youngs_modulus[i] = mat[monomer_matID[i]].E;
    }
    for (int i = 0; i < Nmon; i++) {
        host_matProperties.poisson_number[i] = mat[monomer_matID[i]].nu;
    }
    for (int i = 0; i < Nmon; i++) {
        host_matProperties.damping_timescale[i] = mat[monomer_matID[i]].tvis;
    }
    for (int i = 0; i < Nmon; i++) {
        host_matProperties.crit_rolling_disp[i] = mat[monomer_matID[i]].xi;
    }

    // Push the material values to the device
    matProperties_pushToDevice(host_matProperties, device_matProperties, Nmon);

    // Free memory of initial values.
    free(initial_position);
    free(initial_velocity);
    free(initial_omega_tot);
    free(initial_magnetization);

    // Initialize potential energy tracker
    double4* device_potential_energy = nullptr;
    CHECK_CUDA(cudaMalloc(& device_potential_energy, sizeof(double4)));
    CHECK_CUDA(cudaMemset(device_potential_energy, 0, sizeof(double4)));

    // Initialize inelastic motion counters
    double4* device_inelastic_counter = nullptr;
    CHECK_CUDA(cudaMalloc(& device_inelastic_counter, sizeof(double4)));
    CHECK_CUDA(cudaMemset(device_inelastic_counter, 0, sizeof(double4)));

    // Print run summary.
    pipeline.printParameters();
    
    // Get some run parameters.
    ullong N_iter = pipeline.getNIter(); // Number of total simulation steps.
    ullong N_save = pipeline.getNSave(); // Frequency of snapshots.
    double time_step = pipeline.getTimeStep(); // Timestep of the simulation in [s].
    bool save_ovito = pipeline.saveOvito();
    double3 B_ext = pipeline.getBext();

    int N_store = 0; // The total number of simulation steps where the state is stored.
    int N_store_mon = 0; // The total number of individual monomer states that need to be stored.

    if (N_save > 0)
    {
        N_store = int((double(N_iter) / double(N_save) + 0.5));
        N_store_mon = Nmon * N_store;

        clusterIDs = new int[Nmon];
        fill(clusterIDs, clusterIDs + Nmon, -1);

        if (pipeline.savePos())
        {
            storage_pos = new double3[N_store_mon];
            memset(storage_pos, 0, N_store_mon * sizeof(double3));
        }

        if (pipeline.saveVel())
        {
            storage_vel = new double3[N_store_mon];
            memset(storage_vel, 0, N_store_mon * sizeof(double3));
        }

        if (pipeline.saveForce())
        {
            storage_force = new double3[N_store_mon];
            memset(storage_force, 0, N_store_mon * sizeof(double3));
        }

        if (pipeline.saveTorque())
        {
            storage_torque = new double3[N_store_mon];
            memset(storage_torque, 0, N_store_mon * sizeof(double3));
        }

        if (pipeline.saveOmega())
        {
            storage_omega = new double3[N_store_mon];
            memset(storage_omega, 0, N_store_mon * sizeof(double3));
        }

        if (pipeline.saveMag())
        {
            storage_mag = new double3[N_store_mon];
            memset(storage_mag, 0, N_store_mon * sizeof(double3));
        }

        if (true) {
            storage_potential_N = new double[N_store];
            memset(storage_potential_N, 0, N_store * sizeof(double));
            storage_potential_S = new double[N_store];
            memset(storage_potential_S, 0, N_store * sizeof(double));
            storage_potential_R = new double[N_store];
            memset(storage_potential_R, 0, N_store * sizeof(double));
            storage_potential_T = new double[N_store];
            memset(storage_potential_T, 0, N_store * sizeof(double));
        }

        if (true) {
            storage_inelastic_N = new double[N_store];
            memset(storage_inelastic_N, 0, N_store * sizeof(double));
            storage_inelastic_S = new double[N_store];
            memset(storage_inelastic_S, 0, N_store * sizeof(double));
            storage_inelastic_R = new double[N_store];
            memset(storage_inelastic_R, 0, N_store * sizeof(double));
            storage_inelastic_T = new double[N_store];
            memset(storage_inelastic_T, 0, N_store * sizeof(double));
        }

        if (pipeline.saveCluster())
        {
            storage_cluster = new int[N_store_mon];
            memset(storage_cluster, 0, N_store_mon * sizeof(int));
        }
    }

    // Print an overview of the current compute device.
    PRINT_TITLE("OVERVIEW OF COMPUTE DEVICE")
    PRINT_CLR_LINE();

    int device_count;
    cudaGetDeviceCount(&device_count);
    int active_device_ID;
    cudaGetDevice(&active_device_ID);

    printf("Active device: %d of %d\n", active_device_ID + 1, device_count);

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, active_device_ID);

    printf("   - Name:                  %s\n", prop.name);
    printf("   - Compute capability:    %d.%d\n", prop.major, prop.minor);
    printf("   - Total global mem:      %lu bytes\n", prop.totalGlobalMem);
    printf("   - Warp size:             %d threads\n", prop.warpSize);
    printf("   - Max threads / block:   %d threads\n", prop.maxThreadsPerBlock);
    
    // Calculate the number of blocks the kernels need.
    int nBlocks_single = (Nmon + BLOCK_SIZE - 1) / BLOCK_SIZE; // The number of blocks needed to process all monomers.
    int nBlocks_pair = (Nmon * Nmon + BLOCK_SIZE - 1) / BLOCK_SIZE; // The number of blocks needed to process all monomer pairs.

    // Some helper variables needed during simulation.
    ullong counter_save = 0;        // The number of snapshots that have been saved.
    ullong ns_per_iteration = 0;    // The estimated time each simulation step takes in [ns].

    PRINT_CLR_LINE();
    PRINT_TITLE("SIMULATING");
    PRINT_CLR_LINE();

#ifdef DEBUG
    // This change the following line to enable/disable code that pulls the current and next state during the simulation. This has a significant resource drain and should only be used for debugging.
    #define DEBUG_PULL_STATES
#endif

#ifdef DEBUG_PULL_STATES
    hostState curr, next;
    state_allocateHostMemory(curr, Nmon);
    state_allocateHostMemory(next, Nmon);
#endif

    /* ############################################################
            SIMULATING
       ############################################################ */

    for (ullong iter = 0; iter < N_iter; iter++) // The simulations iteration count.
    {
        auto iteration_start = std::chrono::high_resolution_clock::now();

#ifdef DEBUG
        printf("Starting iteration %llu \n", iter);
        std::cout << std::flush;
#endif

#ifdef DEBUG_PULL_STATES
        state_pullFromDevice(device_state_curr, curr, Nmon);
        state_pullFromDevice(device_state_next, next, Nmon);
#endif

        predictor <<<nBlocks_single, BLOCK_SIZE>>> (
            device_state_curr.position, device_state_curr.velocity, device_state_curr.force,
            device_state_next.position,
            device_matProperties.mass, time_step, Nmon
        );
        
        predictor_pointer <<<nBlocks_pair, BLOCK_SIZE>>> (
            device_state_curr.contact_rotation, device_state_curr.contact_twist, device_state_curr.omega, device_state_curr.torque,
            device_state_next.contact_rotation, device_state_next.contact_twist, 
            device_matProperties.moment, time_step, Nmon
        );

        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

#ifdef DEBUG_PULL_STATES
        state_pullFromDevice(device_state_curr, curr, Nmon);
        state_pullFromDevice(device_state_next, next, Nmon);
#endif

        evaluate <<<nBlocks_pair, BLOCK_SIZE>>> (
            device_state_next.position, device_state_curr.contact_pointer, device_state_next.contact_rotation, device_state_next.contact_twist, device_state_curr.contact_compression,
            device_state_next.force, device_state_next.torque, device_potential_energy,
            device_matProperties.mass, device_matProperties.radius, device_matProperties.youngs_modulus, device_matProperties.poisson_number, device_matProperties.surface_energy, device_matProperties.crit_rolling_disp, device_matProperties.damping_timescale,
            time_step, Nmon
        );
        
        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

#ifdef DEBUG_PULL_STATES
        state_pullFromDevice(device_state_curr, curr, Nmon);
        state_pullFromDevice(device_state_next, next, Nmon);
#endif

        corrector <<<nBlocks_single, BLOCK_SIZE>>> (
            device_state_curr.velocity, device_state_curr.omega,
            device_state_curr.force, device_state_next.force,
            device_state_curr.torque, device_state_next.torque,
            device_state_next.velocity, device_state_next.omega,
            device_matProperties.mass, device_matProperties.moment,
            time_step, Nmon
        );

        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

#ifdef DEBUG_PULL_STATES
        state_pullFromDevice(device_state_curr, curr, Nmon);
        state_pullFromDevice(device_state_next, next, Nmon);
#endif

        updatePointers <<<nBlocks_pair, BLOCK_SIZE>>> (
            device_state_next.position, device_state_curr.contact_pointer, device_state_curr.contact_rotation, device_state_curr.contact_compression,
            device_state_next.contact_pointer, device_state_next.contact_rotation, device_state_next.contact_twist, device_state_next.contact_compression, device_inelastic_counter,
            device_matProperties.radius, device_matProperties.youngs_modulus, device_matProperties.poisson_number, device_matProperties.surface_energy, device_matProperties.crit_rolling_disp,
            Nmon
        );

        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

#ifdef DEBUG_PULL_STATES
        state_pullFromDevice(device_state_curr, curr, Nmon);
        state_pullFromDevice(device_state_next, next, Nmon);
        if (iter >= 0) {
            printf("pos:\n");
            print_double3(next.position, 0, min(Nmon, 5));
            printf("F:\n");
            print_double3(next.force, 0, min(Nmon, 5));
            printf("vel:\n");
            print_double3(next.velocity, 0, min(Nmon, 5));
            printf("M:\n");
            print_double3(next.torque, 0, min(Nmon, 5));
            printf("omega:\n");
            print_double3(next.omega, 0, min(Nmon, 5));
            printf("rotation:\n");
            print_double4(next.contact_rotation, 0, min(Nmon * Nmon, 5));
            printf("pointer:\n");
            print_double3(next.contact_pointer, 0, min(Nmon * Nmon, 5));
        }
        cout << endl;
#endif
        
        // Switch the current and next state.
        deviceState tmp = device_state_curr;
        device_state_curr = device_state_next;
        device_state_next = tmp;

        cudaMemset(device_state_next.force, 0, Nmon * sizeof(double3));
        cudaMemset(device_state_next.torque, 0, Nmon * sizeof(double3));

        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

#ifdef DEBUG_PULL_STATES
        state_pullFromDevice(device_state_curr, curr, Nmon);
        state_pullFromDevice(device_state_next, next, Nmon);
#endif

        // Store state snapshots.
        if (N_save > 0)
        {
            if (iter % N_save == 0)
            {
                // Pull the current state information form the device.
                state_pullFromDevice(device_state_curr, host_state_curr, Nmon);

                ullong start_index = counter_save * Nmon;

                if (long(N_store_mon) - long(start_index) < Nmon)
                {
                    PRINT_ERROR("Storage overrun!")
                    return -1;
                }

                if (start_index < N_store_mon) //just a failsave if there was a rounding error in N_store_mon
                {
                    if (storage_pos != 0)
                    {
                        CHECK_CUDA(cudaMemcpy(storage_pos + start_index, host_state_curr.position, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost))
                    }

                    if (storage_vel != 0)
                    {
                        CHECK_CUDA(cudaMemcpy(storage_vel + start_index, host_state_curr.velocity, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost))
                    }

                    if (storage_force != 0)
                    {
                        CHECK_CUDA(cudaMemcpy(storage_force + start_index, host_state_curr.force, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost))
                    }

                    if (storage_torque != 0)
                    {
                        CHECK_CUDA(cudaMemcpy(storage_torque + start_index, host_state_curr.torque, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost))
                    }

                    if (storage_omega != 0)
                    {
                        CHECK_CUDA(cudaMemcpy(storage_omega + start_index, host_state_curr.omega, Nmon * sizeof(double3), cudaMemcpyKind::cudaMemcpyHostToHost))
                    }

                    if (device_potential_energy != nullptr) {
                        CHECK_CUDA(cudaMemcpy(storage_potential_N + counter_save, &device_potential_energy->w, sizeof(double), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                        CHECK_CUDA(cudaMemcpy(storage_potential_S + counter_save, &device_potential_energy->x, sizeof(double), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                        CHECK_CUDA(cudaMemcpy(storage_potential_R + counter_save, &device_potential_energy->y, sizeof(double), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                        CHECK_CUDA(cudaMemcpy(storage_potential_T + counter_save, &device_potential_energy->z, sizeof(double), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                        // Reset the potential energy counter.
                        CHECK_CUDA(cudaMemset(device_potential_energy, 0, sizeof(double4)));
                    }

                    if (device_inelastic_counter != nullptr) {
                        CHECK_CUDA(cudaMemcpy(storage_inelastic_N + counter_save, &device_inelastic_counter->w, sizeof(double), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                        CHECK_CUDA(cudaMemcpy(storage_inelastic_S + counter_save, &device_inelastic_counter->x, sizeof(double), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                        CHECK_CUDA(cudaMemcpy(storage_inelastic_R + counter_save, &device_inelastic_counter->y, sizeof(double), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                        CHECK_CUDA(cudaMemcpy(storage_inelastic_T + counter_save, &device_inelastic_counter->z, sizeof(double), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                        // Reset the inelastic motion counter.
                        CHECK_CUDA(cudaMemset(device_inelastic_counter, 0, sizeof(double4)));
                    }

                    if (storage_cluster != 0)
                    {
                        fill(clusterIDs, clusterIDs + Nmon, -1);
                        cpu_findConnectedComponents(Nmon, host_state_curr.contact_compression, clusterIDs);
                        copy(clusterIDs, clusterIDs + Nmon, storage_cluster + start_index);
                    }
                }

                // Increment the number of states stored in the storage arrays.
                counter_save++;
            }
        }

        // Calculate the time ellapsed during this iteration.
        auto iteration_current = std::chrono::high_resolution_clock::now();
        auto iteration_ellapsed = std::chrono::duration_cast <std::chrono::nanoseconds> (iteration_current - iteration_start);

        // Use rolling average to estimate the time per iteration.
        float percentage = 100 * float(iter) / N_iter;
        float weight = ROLLING_AVERAGE_WEIGHT;
        
        // The first iteration has volatile timing and should be skipped, this can be achieved with a rolling average weight of 1.0 for the next iteration.
        if (iter < 2) {
            weight = 1.;
        }

        ns_per_iteration = (ulong) ((1.0 - weight) * ((double)ns_per_iteration) + weight * ((double)iteration_ellapsed.count()));

        // Print the simulation progress to console.
        if (((iter - PROGRESS_LOG_OFFSET) % (N_iter / PROGRESS_LOG_AMMOUNT)) == 0) {
            ullong remaining_ns = ns_per_iteration * (N_iter - iter);

            char buffer[14];
            ns_to_time_string(remaining_ns, buffer, 14);

            printf("Simulation progress  :%5.1f %%\n      Remaining time ~ %s\n", percentage, buffer);
        }
    }

    PRINT_CLR_LINE();

    /* ############################################################
            STORING THE RESULTS
       ############################################################ */

    if (N_save > 0)
    {
        PRINT_TITLE("WRITING SIMULATION DATA");
        PRINT_CLR_LINE();

        // Write ouput files
        if (save_ovito && storage_pos != 0)
        {
            if (!pipeline.writeAllOVITO(storage_pos, storage_vel, storage_force, storage_torque, storage_omega, storage_mag, storage_cluster, monomer_radius, monomer_matID, Nmon, N_store_mon))
                return -1;
        }

        if (!pipeline.writeHeader())
            return -1;

        if (!pipeline.writeBinaryDouble("agg_a_mon.bin", monomer_radius, Nmon))
            return -1;

        if (!pipeline.writeBinaryDouble("agg_mass_mon.bin", monomer_mass, Nmon))
            return -1;

        for (int i = 0; i < Nmon; i++)
            monomer_matID[i] = monomer_matID[i] + 1;

        if (!pipeline.writeBinaryInt("agg_matid_mon.bin", monomer_matID, Nmon))
            return -1;

        if (storage_pos != 0)
        {
            if (!pipeline.writeBinaryVec("sim_pos.bin", storage_pos, N_store_mon))
                return -1;
        }

        if (storage_vel != 0)
            if (!pipeline.writeBinaryVec("sim_vel.bin", storage_vel, N_store_mon))
                return -1;

        if (storage_omega != 0)
            if (!pipeline.writeBinaryVec("sim_omega.bin", storage_omega, N_store_mon))
                return -1;
                
        if (storage_force != 0)
            if (!pipeline.writeBinaryVec("sim_force.bin", storage_force, N_store_mon))
                return -1;

        if (storage_torque != 0)
            if (!pipeline.writeBinaryVec("sim_torque.bin", storage_torque, N_store_mon))
                return -1;

        if (storage_cluster != 0)
            if (!pipeline.writeBinaryInt("sim_cluster.bin", storage_cluster, N_store_mon))
                return -1;

        if (storage_potential_N != 0) {
            for (int i = 0; i < N_store; i++) {
                storage_potential_N[i] /= N_save;
                storage_potential_R[i] /= N_save;
                storage_potential_S[i] /= N_save;
                storage_potential_T[i] /= N_save;
            }

            if (!pipeline.writeBinaryDouble("sim_normal_pot.bin", storage_potential_N, N_store)) {
                return -1;
            }
        }

        if (storage_potential_S != 0) {
            if (!pipeline.writeBinaryDouble("sim_sliding_pot.bin", storage_potential_S, N_store)) {
                return -1;
            }
        }

        if (storage_potential_R != 0) {
            if (!pipeline.writeBinaryDouble("sim_rolling_pot.bin", storage_potential_R, N_store)) {
                return -1;
            }
        }

        if (storage_potential_T != 0) {
            if (!pipeline.writeBinaryDouble("sim_twisting_pot.bin", storage_potential_T, N_store)) {
                return -1;
            }
        }

        if (storage_inelastic_N != 0) {
            if (!pipeline.writeBinaryDouble("sim_normal_diss.bin", storage_inelastic_N, N_store)) {
                return -1;
            }
        }
        
        if (storage_inelastic_S != 0) {
            if (!pipeline.writeBinaryDouble("sim_sliding_diss.bin", storage_inelastic_S, N_store)) {
                return -1;
            }
        }
        
        if (storage_inelastic_R != 0) {
            if (!pipeline.writeBinaryDouble("sim_rolling_diss.bin", storage_inelastic_R, N_store)) {
                return -1;
            }
        }
        
        if (storage_inelastic_T != 0) {
            if (!pipeline.writeBinaryDouble("sim_twisting_diss.bin", storage_inelastic_T, N_store)) {
                return -1;
            }
        }
        

        PRINT_CLR_LINE();
    }

    /* ############################################################
            FINAL CLEANUP
       ############################################################ */

    PRINT_TITLE("FINAL CLEANUP");
    PRINT_CLR_LINE();

    // Freeing host and device memory.
    state_freeHost(host_state_curr);
    state_freeHost(host_state_next);
    state_freeDevice(device_state_curr);
    state_freeDevice(device_state_next);
    matProperties_freeHostMemory(host_matProperties);
    matProperties_freeDeviceMemory(device_matProperties);

    free(monomer_radius);
    free(monomer_mass);
    free(monomer_moment);
    free(monomer_matID); 

    // Check if there were CUDA errors during memory deallocation.
    CUDA_LAST_ERROR_CHECK();

    // Calculate and print the final runtime.
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    char buffer[14];
    ns_to_time_string(elapsed.count(), buffer, 14);
    printf("Total runtime : %s .\n", buffer);
    
    PRINT_CLR_LINE();
    PRINT_TITLE("DONE");

    return 0;
}