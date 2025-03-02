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

using namespace std;
using namespace std::chrono;

#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "utils/config.cuh"
#include "utils/constant.cuh"
#include "utils/errors.cuh"
#include "utils/printing.cuh"
#include "utils/typedefs.cuh"

#include "physics/state.cuh"
#include "physics/materialProperties.cuh"
#include "physics/integrator.cuh"

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

    PRINT_CLR_LINE();
    PRINT_TITLE("INITIALIZATION")
    PRINT_CLR_LINE();

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

    /*
    // Definition of variables.
    // The old variables:
    vec3D* vel = nullptr;                 // An array of the monomer velocities.
    vec3D* omega = nullptr;               // An array of the monomer angular velocities due to self rotation.
    vec3D* omega_tot = nullptr;           // An array of the monomer angular velocities due to self rotation and curved trajectories.
    vec3D* mag = nullptr;                 // An array of the monomer magnetizations.

    vec3D* pos_old = nullptr;
    vec3D* pos_new = nullptr;
    vec3D* force_old = nullptr;
    vec3D* force_new = nullptr;
    vec3D* torque_old = nullptr;
    vec3D* torque_new = nullptr;

    vec3D* dMdt_old = nullptr;
    vec3D* dMdt_new = nullptr;

    vec3D* buff_vel = nullptr;
    vec3D* buff_omega = nullptr;
    vec3D* buff_omega_tot = nullptr;
    vec3D* buff_mag = nullptr;

    vec3D* buff_pos_old = nullptr;
    vec3D* buff_pos_new = nullptr;
    vec3D* buff_force_old = nullptr;
    vec3D* buff_force_new = nullptr;
    vec3D* buff_torque_old = nullptr;
    vec3D* buff_torque_new = nullptr;

    vec3D* buff_dMdt_old = nullptr;
    vec3D* buff_dMdt_new = nullptr;
    */

    double3* storage_pos = nullptr;         // Array of monomer position that are to be stored after simulation.
    double3* storage_vel = nullptr;         // Array of monomer velocities that are to be stored after simulation.
    double3* storage_force = nullptr;       // Array of monomer forces that are to be stored after simulation.
    double3* storage_torque = nullptr;      // Array of monomer torques that are to be stored after simulation.
    double3* storage_omega = nullptr;       // Array of monomer angular velocities that are to be stored after simulation.
    double3* storage_mag = nullptr;         // Array of monomer magnetizations that are to be stored after simulation.
    double* storage_inelastic_N = nullptr;
    double* storage_inelastic_S = nullptr;
    double* storage_inelastic_R = nullptr;
    double* storage_inelastic_T = nullptr;
    int* storage_cluster = nullptr;         // Array of monomer cluster membership that are to be stored after simulation.
    int* clusterIDs = nullptr;

    /*
    vec3D* matrix_con = nullptr;          // contact pointer between monomers
    vec3D* matrix_norm = nullptr;         // normal vectors between monomers
    quat* matrix_rot = nullptr;           // contact pointer rotation direction
    double* matrix_comp = nullptr;        // old compression lengths, also used to track connection
    double* matrix_twist = nullptr;       // old twisting displacement

    vec3D* buff_matrix_con = nullptr;     // GPU contact pointer
    vec3D* buff_matrix_norm = nullptr;    // GPU normal vectors
    quat* buff_matrix_rot = nullptr;      // GPU contact pointer rotation
    double* buff_matrix_comp = nullptr;   // GPU old compression lenghts
    double* buff_matrix_twist = nullptr;  // GPU old twisting displacement

    int* matIDs = nullptr;                // material IDs
    double* amon = nullptr;               // Monomer radii
    double* moment = nullptr;             // Monomer moments of inertia
    double* mass = nullptr;               // Monomer masses
    int* clusterIDs = nullptr;            // Monomer cluster membership

    int* buff_matIDs = nullptr;           // GPU material IDs
    double* buff_amon = nullptr;          // GPU monomer radii
    double* buff_moment = nullptr;        // GPU monomer moments of inertia
    double* buff_mass = nullptr;          // GPU monomer masses
    */

    // Read the aggregate files and initialize the state.
    //pipeline.prepareData(pos_old, vel, omega_tot, mag, amon, mass, moment, matIDs, Nmon);

    /*
    PRINT_LOG("Allocating memory on host.", 2);
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

    PRINT_LOG("Allocating memory on device.", 2)
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

    // Check for CUDA errors during memory allocation.
    CUDA_LAST_ERROR_CHECK();

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
    */

    // Check for CUDA errors during memory allocation
    CUDA_LAST_ERROR_CHECK();

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
    int Nmon = 0;                   // Number of monomers

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

    // Initialize inelastic motion counters
    double4* device_inelastic_counter = nullptr;
    CHECK_CUDA(cudaMalloc(& device_inelastic_counter, sizeof(double4)));
    CHECK_CUDA(cudaMemset(device_inelastic_counter, 0, sizeof(double4)));

    // Initialize the energy counter
    double* device_total_energy = nullptr;
    CHECK_CUDA(cudaMalloc(& device_total_energy, sizeof(double)));
    CHECK_CUDA(cudaMemset(device_total_energy, 0, sizeof(double)));

    // Print run summary.
    pipeline.printParameters();
    
    // Get some run parameters.
    ullong N_iter = pipeline.getNIter();
    ullong N_save = pipeline.getNSave();
    double time_step = pipeline.getTimeStep();
    bool save_ovito = pipeline.saveOvito();
    double3 B_ext = pipeline.getBext();

    int N_store = 0; // The number of simulation steps where the state is stored.
    int N_store_mon = 0; // The number of individual monomer states that need to be stored.

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

    // Calculate the number of blocks
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
    
    int nBlocks_single = (Nmon + BLOCK_SIZE - 1) / BLOCK_SIZE; // The number of blocks needed to process all monomers.
    int nBlocks_pair = (Nmon * Nmon + BLOCK_SIZE - 1) / BLOCK_SIZE; // The number of blocks needed to process all monomer pairs.

    /*
    PRINT_LOG("Push initial state to device.", 2);
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

    // Check for CUDA errors during the push of initial state to device.
    CUDA_LAST_ERROR_CHECK();
    */

    ullong counter_save = 0;

    PRINT_CLR_LINE();
    PRINT_TITLE("SIMULATING");
    PRINT_CLR_LINE();

    ullong ns_per_iteration = 0;

    // THE MAIN SIMULATION LOOP
    for (ullong iter = 0; iter < N_iter; iter++) // The simulation iteration count.
    {
        auto iteration_start = std::chrono::high_resolution_clock::now();

#ifdef DEBUG
        #define DEBUG_PULL_STATES
#endif

#ifdef DEBUG_PULL_STATES
        printf("Starting iteration %llu \n", iter);
        hostState curr, next;
        state_allocateHostMemory(curr, Nmon);
        state_allocateHostMemory(next, Nmon);
        state_pullFromDevice(device_state_curr, curr, Nmon);
        state_pullFromDevice(device_state_next, next, Nmon);
        //printf("    Predictor:\n");
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
        //printf("    Evaluator:\n");
#endif

        evaluate <<<nBlocks_pair, BLOCK_SIZE>>> (
            device_state_next.position, device_state_curr.contact_pointer, device_state_next.contact_rotation, device_state_next.contact_twist, device_state_curr.contact_compression,
            device_state_next.force, device_state_next.torque,
            device_matProperties.mass, device_matProperties.radius, device_matProperties.youngs_modulus, device_matProperties.poisson_number, device_matProperties.surface_energy, device_matProperties.crit_rolling_disp, device_matProperties.damping_timescale,
            time_step, Nmon
        );
        
        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

#ifdef DEBUG_PULL_STATES
        state_pullFromDevice(device_state_curr, curr, Nmon);
        state_pullFromDevice(device_state_next, next, Nmon);
        //printf("    Corrector:\n");
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
        //printf("    Pointer Update:\n");
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
        //printf("    Switch:\n");
#endif
        
        // Switch pointers
        //PRINT_LOG("Pointer switch", 1);
        deviceState tmp = device_state_curr;
        device_state_curr = device_state_next;
        device_state_next = tmp;

        // TODO: Check if other fields need to be zeroed
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

                    if (device_inelastic_counter!= nullptr) {
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

    PRINT_TITLE("FINAL CLEANUP");
    PRINT_CLR_LINE();

    /*
    PRINT_LOG("Freeing host memory.", 3);
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

    PRINT_LOG("Freeing device memory.\n", 3);

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
    */

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
    auto end = high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    char buffer[14];
    ns_to_time_string(elapsed.count(), buffer, 14);
    printf("Total runtime : %s .\n", buffer);
    
    PRINT_CLR_LINE();
    PRINT_TITLE("DONE");

    return 0;
}