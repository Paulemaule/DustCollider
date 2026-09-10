#pragma once

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <cuda_runtime.h>

#include "utils/logging.cuh"
#include "utils/config.cuh"
#include "utils/errors.cuh"
#include "utils/printing.cuh"

#include "physics/state.cuh"
#include "physics/materials.cuh"
#include "physics/integrator.cuh"

#include "simulationSetup/simulationConfig.cuh"

/**
 * @brief Owns all runtime memory and drives the PEC time loop.
 *
 * Constructed from a fully-initialised SimulationConfig (moved in).
 * Call run() to execute the simulation, then write_output() to flush snapshots.
 */
class Simulator {
    SimulationConfig      config;
    size_t                Nmon;

    // Single host buffer
    HostState             host_state;

    // RAII buffers for device memory
    DeviceState           device_curr;      // The current system state in device memory
    DeviceState           device_next;      // The next system state in device memory

    HostMaterials         host_materials;   // The material properties of the monomers in host memory.
    DeviceMaterials       device_materials; // The material properties of the monomers in device memory.

    // RAII buffers for the energy accumulators
    DeviceBuffer<double4> device_potential; // Tracker for the potential energies in device memory.
    DeviceBuffer<double4> device_inelastic; // Accumulator for the energy loss in device memory.

    // Storage for the system state snapshots
    size_t N_store_ = 0;                    // TODO

    // Kinematic snapshots — allocated only when the matching output flag is set
    std::optional<std::vector<double3>> snap_pos_;
    std::optional<std::vector<double3>> snap_vel_;
    std::optional<std::vector<double3>> snap_force_;
    std::optional<std::vector<double3>> snap_torque_;
    std::optional<std::vector<double3>> snap_omega_;

    // Energy / dissipation — always stored when N_store_ > 0
    std::vector<double> snap_pot_N_, snap_pot_S_, snap_pot_R_, snap_pot_T_;
    std::vector<double> snap_dis_N_, snap_dis_S_, snap_dis_R_, snap_dis_T_;

    // Cluster IDs — always stored when N_store_ > 0
    std::vector<int> snap_cluster_;

public:
    /**
     * @brief Constructs and fully initialises the simulator.
     *
     * Allocates all host and device memory, copies the initial state and material
     * properties to the device, and pre-allocates snapshot storage.
     *
     * @param cfg  Completed SimulationConfig (taken by value and moved in).
     */
    explicit Simulator(SimulationConfig cfg)
        : config(std::move(cfg))
        , Nmon(config.initial_state.positions.size())
        , host_state(Nmon)
        , device_curr(Nmon)
        , device_next(Nmon)
        , host_materials(Nmon)
        , device_materials(Nmon)
        , device_potential(1)
        , device_inelastic(1)
    {
        init_state();
        init_mat();
        allocate_snapshots();
        log_device_info();
    }

    void run();
    void write_output() const;

private:
    void log_device_info() const;
    void init_state();
    void init_mat();
    void allocate_snapshots();
    void save_snapshot(size_t snap_idx);
    void write_ovito() const;
};

/**
 * @brief Logs the GPU architectures this binary was compiled for and checks for compatibility with the active compute device.
 * 
 * The compiler builds a fat binary holding machine code for all specified architectures (see arch.mk).
 * Additionally forward support for newer architectures through JIT compilation is included.
 * This function logs the supported architectures and checks if the active compute device is compatible.
 *
 * @param prop The properties of the active compute device.
 */
inline void log_compiled_architectures(const cudaDeviceProp& prop) {
#ifdef __CUDA_ARCH_LIST__
    // nvcc defines __CUDA_ARCH_LIST__ as the compute capabilities present in the binary
    constexpr int compiled[] = { __CUDA_ARCH_LIST__ };

    // Convert the list of compute capabilities from the compile flag into a string
    std::string arch_list;
    for ( const int arch : compiled ) {
        if ( !arch_list.empty() ) arch_list += ", ";
        arch_list += std::to_string(arch / 100) + "." + std::to_string((arch / 10) % 10);
    }

    // Adjust the compute capability of the current active device to the same format
    const int device_arch = prop.major * 100 + prop.minor * 10;

    int  newest             = 0;      // The newest architecture in the binary, its PTX is the fallback.
    bool exact_match        = false;  // The binary holds machine code for this exact device.
    bool compatible_match   = false;  // The binary holds machine code of the same GPU generation.

    for ( const int arch : compiled ) {
        if ( arch == device_arch ) exact_match = true;
        // Machine code is compatible upwards within a generation: code for 8.0 also runs on 8.6.
        if ( arch / 100 == prop.major && arch < device_arch ) compatible_match = true;
        if ( arch > newest ) newest = arch;
    }

    // Log the results
    Logger::print("   compute capabilities:          {}+", arch_list);

    if ( exact_match ) return;

    if ( compatible_match ) {
        Logger::warn("The binary contains no machine code for compute capability {}.{}. The device runs code "
                     "of an older minor version, which is not tuned for it. Add {} to CUDA_ARCHS in arch.mk.",
                     prop.major, prop.minor, device_arch / 10);
    } else if ( device_arch > newest ) {
        Logger::warn("This GPU is newer than every architecture the binary was compiled for. The kernels are "
                     "JIT-compiled from PTX on startup, which costs time. Add {} to CUDA_ARCHS in arch.mk.",
                     device_arch / 10);
    } else {
        Logger::error("This binary can not run on compute capability {}.{}. Rebuild with {} added to "
                      "CUDA_ARCHS in arch.mk.", prop.major, prop.minor, device_arch / 10);
    }
#else
    // Toolkits older than CUDA 11.5 do not expose the architecture list to host code.
    Logger::warn("The CUDA version used to compile this code is too old. CUDA 13.0+ is recommended!")
    (void)prop;
#endif
}

/**
 * @brief Prints an overview of the active CUDA compute devices.
 */
inline void Simulator::log_device_info() const {
    Logger::lineBreak();
    Logger::header("OVERVIEW OF COMPUTE DEVICE");
    Logger::lineBreak();

    // Retrieve device information
    int device_count = 0;
    CHECK_CUDA(cudaGetDeviceCount(&device_count));
    int active_device_id = 0;
    CHECK_CUDA(cudaGetDevice(&active_device_id));

    cudaDeviceProp prop;
    CHECK_CUDA(cudaGetDeviceProperties(&prop, active_device_id));

    // Log device information
    Logger::print("Active device: {} of {}", active_device_id + 1, device_count);
    Logger::print("   Name:                  {}", prop.name);
    Logger::print("   Compute capability:    {}.{}", prop.major, prop.minor);
    Logger::print("   Total global mem:      {} bytes", prop.totalGlobalMem);
    Logger::print("   Warp size:             {} threads", prop.warpSize);
    Logger::print("   Max threads / block:   {} threads", prop.maxThreadsPerBlock);

    // Print a log comparing the active devices compute capability against the 
    // available compute capabilites in the compiled fat binary.
    Logger::lineBreak();
    Logger::print("Code compiled for: ");
    log_compiled_architectures(prop);

    Logger::lineBreak();
}

/**
 * @brief Prepares the initial system state in device memory.
 *
 * This function copies the initial state from the simulation config into a host state,
 * initializes additional values properly and then pushes it into the device memory.
 * Initializes the energy trackers.
 */
inline void Simulator::init_state() {
    Logger::log("Preparing initial system state in device memory.");
    
    const InitialState& is = config.initial_state;
    HostStateView hv = host_state.view();

    // Copy the initial state from the config into pinned host memory.
    std::copy(is.positions.begin(),  is.positions.end(),  hv.position);
    std::copy(is.velocities.begin(), is.velocities.end(), hv.velocity);
    std::copy(is.omegas.begin(),     is.omegas.end(),     hv.omega);

    // Initialize the host state with the proper values.
    std::fill(hv.force,               hv.force   + Nmon,              double3{});
    std::fill(hv.torque,              hv.torque  + Nmon,              double3{});
    std::fill(hv.contact_compression, hv.contact_compression + Nmon * Nmon, -1.0);
    std::fill(hv.contact_twist,       hv.contact_twist       + Nmon * Nmon,  0.0);
    std::fill(hv.contact_pointer,     hv.contact_pointer     + Nmon * Nmon, double3{});
    std::fill(hv.contact_normal,      hv.contact_normal      + Nmon * Nmon, double3{});
    std::fill(hv.contact_rotation,    hv.contact_rotation    + Nmon * Nmon, double4{});

    // Push the initial system state to device memory.
    host_state.push_to(device_curr, Nmon);
    host_state.push_to(device_next, Nmon);

    // Initialize the energy trackes
    CHECK_CUDA(cudaMemset(device_potential.data(), 0, sizeof(double4)));
    CHECK_CUDA(cudaMemset(device_inelastic.data(), 0, sizeof(double4)));
}

/**
 * @brief Prepares the material parameters in device memory.
 * 
 * This function takes the per material properties from the simulation config 
 * and prepares a per monomer properties struct in pinned host memory.
 * The material properties are then pushed to device.
 */
inline void Simulator::init_mat() {
    Logger::log("Preparing full material monomer properties in device memory.");

    const InitialState& is = config.initial_state;
    HostMaterialsView mv = host_materials.view();

    // Copy the monomer properties from the initial state into pinned host memory.
    std::copy(is.radii.begin(),        is.radii.end(),        mv.radius);
    std::copy(is.masses.begin(),       is.masses.end(),       mv.mass);
    std::copy(is.moments.begin(),      is.moments.end(),      mv.moment);
    std::copy(is.material_ids.begin(), is.material_ids.end(), mv.matID);

    // Lookup the material properties and place them into pinned host memory.
    for (size_t i = 0; i < Nmon; i++) {
        const MaterialEntry& mat = config.materials[is.material_ids[i]];
        mv.density[i]           = mat.rho;
        mv.surface_energy[i]    = mat.gamma;
        mv.youngs_modulus[i]    = mat.E;
        mv.poisson_number[i]    = mat.nu;
        mv.damping_timescale[i] = mat.tvis;
        mv.crit_rolling_disp[i] = mat.xi;
    }

    // Push the material properties into device memory.
    host_materials.push_to(device_materials, Nmon);
}

/**
 * @brief Allocates host memory for the snapshots and other diagnostic information.
 */
inline void Simulator::allocate_snapshots() {
    Logger::log("Allocating memory for snapshot and diagnostic storage.");

    if (config.output.N_save <= 0) {
        Logger::warn("N_save = {} => no snapshot storage will be allocated.", config.output.N_save);
        return;
    }

    // Calculate the final size of the snapshot vectors.
    N_store_ = (size_t(config.N_iter - 1) / size_t(config.output.N_save)) + 1;
    const size_t N_store_mon = Nmon * N_store_;

    // Preallocate memory for the snapshot vectors
    if (config.output.position) snap_pos_.emplace(N_store_mon,    double3{});
    if (config.output.velocity) snap_vel_.emplace(N_store_mon,    double3{});
    if (config.output.force)    snap_force_.emplace(N_store_mon,  double3{});
    if (config.output.torque)   snap_torque_.emplace(N_store_mon, double3{});
    if (config.output.angular)  snap_omega_.emplace(N_store_mon,  double3{});

    // Preallocate and initialize memory for the energy trackers.
    snap_pot_N_.assign(N_store_, 0.0); snap_pot_S_.assign(N_store_, 0.0);
    snap_pot_R_.assign(N_store_, 0.0); snap_pot_T_.assign(N_store_, 0.0);
    snap_dis_N_.assign(N_store_, 0.0); snap_dis_S_.assign(N_store_, 0.0);
    snap_dis_R_.assign(N_store_, 0.0); snap_dis_T_.assign(N_store_, 0.0);

    // Preallocate and initialize memory for the cluster membership.
    snap_cluster_.assign(N_store_mon, -1);

    // Sum up the allocated memory for logging.
    size_t n_vec3 = 0;
    if (snap_pos_)    n_vec3++;
    if (snap_vel_)    n_vec3++;
    if (snap_force_)  n_vec3++;
    if (snap_torque_) n_vec3++;
    if (snap_omega_)  n_vec3++;

    const size_t bytes = n_vec3 * N_store_mon * sizeof(double3)  // kinematic snapshots
                       + 8 * N_store_ * sizeof(double)           // energy trackers
                       + N_store_mon * sizeof(int);              // cluster membership

    Logger::log("Allocated {} snapshot slots ({} monomers each) using {} of host memory.",
                N_store_, Nmon, bytes_to_string(bytes));
}

/**
 * @brief Save the current system state as a snapshot into memory.
 * 
 * Pulls the system state from device memory and writes them into the corresponding
 * snap_* buffers starting at the position snap_idx * Nmon.
 * Also records the energy accumulators 
 * and calculates and stores the monomer aggregate memberships.
 * 
 * @param snap_idx 0 based index of the snapshot to write.
 */
inline void Simulator::save_snapshot(size_t snap_idx) {
    // Pull the system state from device memory
    host_state.pull_from(device_curr, Nmon);
    HostStateView hv = host_state.view();

    // Calculate the index of the snapshot in the storage structs
    const size_t base = snap_idx * Nmon;

    // Copy the state data into the storage structs
    if (snap_pos_)    std::copy(hv.position, hv.position + Nmon, snap_pos_->data()    + base);
    if (snap_vel_)    std::copy(hv.velocity, hv.velocity + Nmon, snap_vel_->data()    + base);
    if (snap_force_)  std::copy(hv.force,    hv.force    + Nmon, snap_force_->data()  + base);
    if (snap_torque_) std::copy(hv.torque,   hv.torque   + Nmon, snap_torque_->data() + base);
    if (snap_omega_)  std::copy(hv.omega,    hv.omega    + Nmon, snap_omega_->data()  + base);

    // Copy the energy diagnostics from device memory into the storage structs
    double4 pot{};
    CHECK_CUDA(cudaMemcpy(&pot, device_potential.data(), sizeof(double4), cudaMemcpyDeviceToHost));
    snap_pot_N_[snap_idx] = pot.w;
    snap_pot_S_[snap_idx] = pot.x;
    snap_pot_R_[snap_idx] = pot.y;
    snap_pot_T_[snap_idx] = pot.z;
    CHECK_CUDA(cudaMemset(device_potential.data(), 0, sizeof(double4)));

    double4 dis{};
    CHECK_CUDA(cudaMemcpy(&dis, device_inelastic.data(), sizeof(double4), cudaMemcpyDeviceToHost));
    snap_dis_N_[snap_idx] = dis.w;
    snap_dis_S_[snap_idx] = dis.x;
    snap_dis_R_[snap_idx] = dis.y;
    snap_dis_T_[snap_idx] = dis.z;
    CHECK_CUDA(cudaMemset(device_inelastic.data(), 0, sizeof(double4)));

    // Calculate monomer cluster membership and store it into a storage struct
    findMonomerClusters((int)Nmon, hv.contact_pointer, snap_cluster_.data() + base);
}

/**
 * @brief Run the PEC simulation loop.
 */
inline void Simulator::run() {
    // The number of thread blocks necessary for kernels with Nmon threads. Conceptually these are kernels acting on single monomers.
    const int nBlocks_single    = int((Nmon + BLOCK_SIZE - 1) / BLOCK_SIZE);
    // The number of thread blcoks necessary for kernels with Nmon * Nmon threads. Conceptually these are kernels acting on monomer pairs.
    const int nBlocks_pair      = int((Nmon * Nmon + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // The number of snapshots stored in the storage structs.
    size_t              counter_save = 0;
    // Timing variable that containes a moving average of the computation time each iteration takes.
    unsigned long long  ns_per_iter  = 0;

    Logger::header("SIMULATING");
    Logger::lineBreak();
    
    // The main simulation loop
    for (int iter = 0; iter < config.N_iter; iter++) {
        auto iter_start = std::chrono::high_resolution_clock::now();

        // Views of raw pointers over the device memory.
        DeviceStateView     curr = device_curr.view();      // View over the current (t) system state in device memory.
        DeviceStateView     next = device_next.view();      // View over the next (t+dt) system state in device memory.
        DeviceMaterialsView mat = device_materials.view();  // View over the material properties of the monomers in device memory.

        // The PREDICTION steps
        predictor<<<nBlocks_single, BLOCK_SIZE>>>(
            curr.position, curr.velocity, curr.force,
            next.position,
            mat.mass, config.timestep, (int)Nmon
        );

        predictor_pointer<<<nBlocks_pair, BLOCK_SIZE>>>(
            curr.contact_rotation, curr.contact_twist,
            next.position, curr.omega, curr.torque,
            next.contact_rotation, next.contact_twist,
            mat.moment, config.timestep, (int)Nmon
        );

        cudaDeviceSynchronize(); // Likely has no effect as synchronize is implied.
        CUDA_LAST_ERROR_CHECK();

        // The EVALUATION step
        evaluate<<<nBlocks_pair, BLOCK_SIZE>>>(
            next.position, curr.contact_pointer,
            next.contact_rotation, next.contact_twist, curr.contact_compression,
            next.force, next.torque,
            device_potential.data(), device_inelastic.data(),
            mat.mass, mat.radius, mat.youngs_modulus, mat.poisson_number,
            mat.surface_energy, mat.crit_rolling_disp, mat.damping_timescale,
            config.timestep, (int)Nmon
        );

        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

        // The CORRECTOR step
        corrector<<<nBlocks_single, BLOCK_SIZE>>>(
            curr.velocity, curr.omega,
            curr.force, next.force,
            curr.torque, next.torque,
            next.velocity, next.omega,
            mat.mass, mat.moment,
            config.timestep, (int)Nmon
        );

        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

        // Update the contact information based on the new positions.
        // FIXME: Incorrect order — updatePointers should run before evaluate (see main.cu)
        updatePointers<<<nBlocks_pair, BLOCK_SIZE>>>(
            next.position, curr.contact_pointer,
            curr.contact_rotation, curr.contact_compression,
            next.contact_pointer, next.contact_rotation,
            next.contact_twist, next.contact_compression,
            device_inelastic.data(),
            mat.radius, mat.youngs_modulus, mat.poisson_number,
            mat.surface_energy, mat.crit_rolling_disp,
            (int)Nmon
        );

        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

        // Swap buffers and zero forces and torques.
        std::swap(device_curr, device_next);
        {
            DeviceStateView nv = device_next.view();
            CHECK_CUDA(cudaMemset(nv.force,  0, Nmon * sizeof(double3)));
            CHECK_CUDA(cudaMemset(nv.torque, 0, Nmon * sizeof(double3)));
        }

        cudaDeviceSynchronize();
        CUDA_LAST_ERROR_CHECK();

        // Store a snapshot if scheduled for this iteration
        if (N_store_ > 0 && size_t(iter) % size_t(config.output.N_save) == 0) {
            save_snapshot(counter_save);
            counter_save++;
        }

        // Progress reporting
        auto iter_end = std::chrono::high_resolution_clock::now();
        unsigned long long iter_ns = (unsigned long long)
            std::chrono::duration_cast<std::chrono::nanoseconds>(iter_end - iter_start).count();

        // The rolling average weight is set to 1 for the first two iterations as they are very volatile and would throw the average off.
        const double weight = (iter < 2) ? 1.0 : ROLLING_AVERAGE_WEIGHT;
        ns_per_iter = (unsigned long long)((1.0 - weight) * (double)ns_per_iter + weight * (double)iter_ns);

        // Print the progress if scheduled.
        if (config.N_iter >= PROGRESS_LOG_AMOUNT) {
            const int step = config.N_iter / PROGRESS_LOG_AMOUNT;
            if (step > 0 && ((iter - PROGRESS_LOG_OFFSET) % step) == 0) {
                float pct = 100.f * float(iter) / float(config.N_iter);
                unsigned long long remaining_ns = ns_per_iter * (unsigned long long)(config.N_iter - iter);
                char buf[14];
                ns_to_time_string(remaining_ns, buf, 14);
                printf("Simulation progress  :%5.1f %%\n      Remaining time ~ %s\n", pct, buf);
                std::cout << std::flush;
            }
        }
    }

    Logger::lineBreak();
}

/**
 * @brief Writes Ovito readable file for visual analysis.
 * 
 * The state variables are all scaled for better visualisation.
 * Velocity, torque and angular velocity are stored as unit vectors.
 * The force is scaled by |F|^(1/8).
 */
inline void Simulator::write_ovito() const {
    Logger::log("Writing Ovito files.");

    const auto ovito_path = std::filesystem::path(config.output.path) / "ovito";
    std::filesystem::create_directories(ovito_path);

    const size_t N_store_mon = Nmon * N_store_;
    const std::vector<double3>& pos = *snap_pos_;

    // Bounding box: largest (distance from origin + radius) across all snapshots [m],
    // then converted to nm with a 15 % margin.
    double b_size = 0.0;
    for (size_t i = 0; i < N_store_mon; i++) {
        const double r   = config.initial_state.radii[i % Nmon];
        const double3& p = pos[i];
        const double   d = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z) + r;
        if (d > b_size) b_size = d;
    }
    b_size *= 1.15e9;
    if (b_size > 10000.0) b_size = 10000.0;

    char box_line[64];
    std::snprintf(box_line, sizeof(box_line), "%.4f %.4f\n", -b_size, b_size);

    // Iterate over all available snapshots
    for (size_t s = 0; s < N_store_; s++) {
        char fname[32];
        std::snprintf(fname, sizeof(fname), "t_%05d.dump", (int)s);
        const std::string filepath = (ovito_path / fname).string();

        std::ofstream writer(filepath);
        if (!writer) {
            Logger::error("Failed to open OVITO file target: '{}'", filepath);
            continue;
        }

        writer << "ITEM: TIMESTEP\n" << s << "\n";
        writer << "ITEM: NUMBER OF ATOMS\n" << Nmon << "\n";
        writer << "ITEM: BOX BOUNDS pp pp pp\n";
        writer << box_line << box_line << box_line;
        writer << "ITEM: ATOMS id mol type x y z vx vy vz fx fy fz tqx tqy tqz omegax omegay omegaz radius\n";

        const size_t base = s * Nmon;

        for (size_t j = 0; j < Nmon; j++) {
            const double x = 1.0e9 * pos[base + j].x;
            const double y = 1.0e9 * pos[base + j].y;
            const double z = 1.0e9 * pos[base + j].z;

            const int cl_id  = snap_cluster_.empty() ? 0 : snap_cluster_[base + j];
            const int mat_id = config.initial_state.material_ids[j] + 1;
            const double r   = 1.0e9 * config.initial_state.radii[j];

            // Store velocity as a unit vector
            double vx = 0, vy = 0, vz = 0;
            if (snap_vel_) {
                double3 v = (*snap_vel_)[base + j];
                const double len = vec_length(v);
                if (len > 0.0) { vx = v.x/len; vy = v.y/len; vz = v.z/len; }
            }

            // Store force scaled by |F|^(1/8) for visual compression
            double fx = 0, fy = 0, fz = 0;
            if (snap_force_) {
                const double3& fv  = (*snap_force_)[base + j];
                const double   len = vec_length(fv);
                if (len > 0.0) {
                    const double s8 = std::pow(len, 1.0/8.0);
                    fx = fv.x/len*s8; fy = fv.y/len*s8; fz = fv.z/len*s8;
                }
            }

            // Store torque as unit vector
            double tx = 0, ty = 0, tz = 0;
            if (snap_torque_) {
                double3 tv = (*snap_torque_)[base + j];
                vec_normalize(tv);
                tx = tv.x; ty = tv.y; tz = tv.z;
            }

            // Angular velocity as unit vector
            double ox = 0, oy = 0, oz = 0;
            if (snap_omega_) {
                double3 ov = (*snap_omega_)[base + j];
                vec_normalize(ov);
                ox = ov.x; oy = ov.y; oz = ov.z;
            }

            // Write data to file
            char line[256];
            std::snprintf(line, sizeof(line),
                "%d %d %d %.5f %.5f %.5f %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5f\n",
                (int)j, cl_id, mat_id, x, y, z,
                vx, vy, vz, fx, fy, fz, tx, ty, tz, ox, oy, oz, r);
            writer << line;
        }
    }
}

/**
 * @brief Write binary files for the system snapshots and energy diagnostics.
 */
inline void Simulator::write_output() const {
    Logger::header("WRITING SIMULATION DATA");
    Logger::lineBreak();

    if (N_store_ == 0) return;

    // Determine file location
    const auto bin_path = std::filesystem::path(config.output.path) / "binary";
    std::filesystem::create_directories(bin_path);
    const std::string bin = bin_path.string();

    // Helper function that write a specific type of data
    auto write_vec3 = [&](const std::string& name, const std::vector<double3>& v) {
        std::ofstream f(bin + "/" + name, std::ios::binary);
        if (!f) { Logger::error("Failed to open binary target file: '{}'", name); }
        f.write(reinterpret_cast<const char*>(v.data()), (std::streamsize)(v.size() * sizeof(double3)));
    };

    auto write_double = [&](const std::string& name, const std::vector<double>& v) {
        std::ofstream f(bin + "/" + name, std::ios::binary);
        if (!f) { Logger::error("Failed to open binary target file: '{}'", name); }
        f.write(reinterpret_cast<const char*>(v.data()), (std::streamsize)(v.size() * sizeof(double)));
    };

    auto write_int = [&](const std::string& name, const std::vector<int>& v) {
        std::ofstream f(bin + "/" + name, std::ios::binary);
        if (!f) { Logger::error("Failed to open binary target file: '{}'", name); }
        f.write(reinterpret_cast<const char*>(v.data()), (std::streamsize)(v.size() * sizeof(int)));
    };

    // Write the monomer data
    // TODO: Write header.txt (N_iter, Nmon, N_save, timestep, N_mat)
    Logger::log("Writing snapshots to disk.");

    write_double("agg_a_mon.bin",    config.initial_state.radii);
    write_double("agg_mass_mon.bin", config.initial_state.masses);
    {
        std::vector<int> matids_1indexed(Nmon);
        for (size_t i = 0; i < Nmon; i++)
            matids_1indexed[i] = config.initial_state.material_ids[i] + 1;
        write_int("agg_matid_mon.bin", matids_1indexed);
    }

    // Write the systam state snapshots
    if (snap_pos_)    write_vec3("sim_pos.bin",    *snap_pos_);
    if (snap_vel_)    write_vec3("sim_vel.bin",    *snap_vel_);
    if (snap_force_)  write_vec3("sim_force.bin",  *snap_force_);
    if (snap_torque_) write_vec3("sim_torque.bin", *snap_torque_);
    if (snap_omega_)  write_vec3("sim_omega.bin",  *snap_omega_);

    if (!snap_cluster_.empty())
        write_int("sim_cluster.bin", snap_cluster_);

    // Write Ovito files
    if ( config.output.ovito ) {
        if ( snap_pos_ )
            write_ovito();
        else
            Logger::warn("Cannot write Ovito output: No position data available.");
    }

    // Write the energy diagnostics
    Logger::log("Writing energy diagnostics to disk.");
    
    if (!snap_pot_N_.empty()) {
        // The potential energies are accumulated over all iterations and will need to be averaged
        const double inv = 1.0 / config.output.N_save;
        std::vector<double> avg_N(N_store_), avg_S(N_store_), avg_R(N_store_), avg_T(N_store_);
        for (size_t i = 0; i < N_store_; i++) {
            avg_N[i] = snap_pot_N_[i] * inv;
            avg_S[i] = snap_pot_S_[i] * inv;
            avg_R[i] = snap_pot_R_[i] * inv;
            avg_T[i] = snap_pot_T_[i] * inv;
        }

        write_double("sim_normal_pot.bin",   avg_N);
        write_double("sim_sliding_pot.bin",  avg_S);
        write_double("sim_rolling_pot.bin",  avg_R);
        write_double("sim_twisting_pot.bin", avg_T);

        write_double("sim_normal_diss.bin",   snap_dis_N_);
        write_double("sim_sliding_diss.bin",  snap_dis_S_);
        write_double("sim_rolling_diss.bin",  snap_dis_R_);
        write_double("sim_twisting_diss.bin", snap_dis_T_);
    }

    Logger::lineBreak();
}
