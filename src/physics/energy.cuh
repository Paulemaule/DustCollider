/**
 * @file energy.cuh
 * @brief Containers for the energy trackers: the potential and dissipated energy split by dofs and dissipation channels.
 * 
 * The energies are tracked by monomer. A pair kernel thread (i,j) adds half of the energy of the pair into the slot of monomer i, 
 * its partner thread (j,i) adds the other half into the slot of j. The host copies the slots into the snapshot storage when a snapshot is saved and then resets them.
 *
 * The trackers are:
 *  - *_pot:       The potential energies stored in the contacts, evaluated on the stored state at every snapshot.
 *  - *_damp:      Energy dissipated through the viscous damping of the normal motion.
 *  - *_inelastic: When a contact exceeds its critical displacement in tangential motion, the contact moves inelastically, dissipating energy.
 *  - *_break:     Contact breaking, all energy stored in the contact (in all four dofs) is lost.
 *  - *_form:      When a contact forms the potential jumps from 0 to U_N(delta) < 0. This is booked here.
 */

#pragma once

#include "../utils/errors.cuh"
#include "../utils/buffer.cuh"

/**
 * @brief The energy trackers of one monomer, one record per monomer and snapshot.
 */
struct EnergyRecord {
    double normal_pot;
    double sliding_pot;
    double rolling_pot;
    double twisting_pot;
    double normal_damp;
    double sliding_inelastic;
    double rolling_inelastic;
    double twisting_inelastic;
    double normal_break;
    double sliding_break;
    double rolling_break;
    double twisting_break;
    double normal_form;
};

/**
 * @brief Non-owning view of the device energy tracker arrays. Passed directly to kernels.
 */
struct DeviceEnergyView {
    double* normal_pot;                     // Array - Nmon * sizeof(<T>)
    double* sliding_pot;                    // Array - Nmon * sizeof(<T>)
    double* rolling_pot;                    // Array - Nmon * sizeof(<T>)
    double* twisting_pot;                   // Array - Nmon * sizeof(<T>)
    double* normal_damp;                    // Array - Nmon * sizeof(<T>)
    double* sliding_inelastic;              // Array - Nmon * sizeof(<T>)
    double* rolling_inelastic;              // Array - Nmon * sizeof(<T>)
    double* twisting_inelastic;             // Array - Nmon * sizeof(<T>)
    double* normal_break;                   // Array - Nmon * sizeof(<T>)
    double* sliding_break;                  // Array - Nmon * sizeof(<T>)
    double* rolling_break;                  // Array - Nmon * sizeof(<T>)
    double* twisting_break;                 // Array - Nmon * sizeof(<T>)
    double* normal_form;                    // Array - Nmon * sizeof(<T>)
};

/**
 * @brief RAII container for the device memory of the energy trackers.
 *
 * Automatically allocates all arrays on construction and frees them on destruction.
 * To pass the data to the kernel use ::view() to obtain a struct of raw pointers.
 */
class DeviceEnergy {
    DeviceBuffer<double> normal_pot;
    DeviceBuffer<double> sliding_pot;
    DeviceBuffer<double> rolling_pot;
    DeviceBuffer<double> twisting_pot;
    DeviceBuffer<double> normal_damp;
    DeviceBuffer<double> sliding_inelastic;
    DeviceBuffer<double> rolling_inelastic;
    DeviceBuffer<double> twisting_inelastic;
    DeviceBuffer<double> normal_break;
    DeviceBuffer<double> sliding_break;
    DeviceBuffer<double> rolling_break;
    DeviceBuffer<double> twisting_break;
    DeviceBuffer<double> normal_form;

public:
    /**
     * @brief Allocates device memory for the energy trackers of a system with Nmon monomers.
     * @param Nmon The number of monomers in the system.
     */
    explicit DeviceEnergy(size_t Nmon)
        : normal_pot(Nmon),        sliding_pot(Nmon),       rolling_pot(Nmon),   twisting_pot(Nmon)
        , normal_damp(Nmon)
        , sliding_inelastic(Nmon), rolling_inelastic(Nmon), twisting_inelastic(Nmon)
        , normal_break(Nmon),      sliding_break(Nmon),     rolling_break(Nmon), twisting_break(Nmon)
        , normal_form(Nmon)
    {}

    /** @brief Generates a view over the device energy trackers for passing into kernels. */
    DeviceEnergyView view() {
        return {
            normal_pot.data(),        sliding_pot.data(),       rolling_pot.data(),   twisting_pot.data(),
            normal_damp.data(),
            sliding_inelastic.data(), rolling_inelastic.data(), twisting_inelastic.data(),
            normal_break.data(),      sliding_break.data(),     rolling_break.data(), twisting_break.data(),
            normal_form.data()
        };
    }

    /** @brief Sets all energy trackers to zero. */
    void zero() {
        CHECK_CUDA(cudaMemset(normal_pot.data(),         0, normal_pot.size_bytes()));
        CHECK_CUDA(cudaMemset(sliding_pot.data(),        0, sliding_pot.size_bytes()));
        CHECK_CUDA(cudaMemset(rolling_pot.data(),        0, rolling_pot.size_bytes()));
        CHECK_CUDA(cudaMemset(twisting_pot.data(),       0, twisting_pot.size_bytes()));
        CHECK_CUDA(cudaMemset(normal_damp.data(),        0, normal_damp.size_bytes()));
        CHECK_CUDA(cudaMemset(sliding_inelastic.data(),  0, sliding_inelastic.size_bytes()));
        CHECK_CUDA(cudaMemset(rolling_inelastic.data(),  0, rolling_inelastic.size_bytes()));
        CHECK_CUDA(cudaMemset(twisting_inelastic.data(), 0, twisting_inelastic.size_bytes()));
        CHECK_CUDA(cudaMemset(normal_break.data(),       0, normal_break.size_bytes()));
        CHECK_CUDA(cudaMemset(sliding_break.data(),      0, sliding_break.size_bytes()));
        CHECK_CUDA(cudaMemset(rolling_break.data(),      0, rolling_break.size_bytes()));
        CHECK_CUDA(cudaMemset(twisting_break.data(),     0, twisting_break.size_bytes()));
        CHECK_CUDA(cudaMemset(normal_form.data(),        0, normal_form.size_bytes()));
    }
};

/**
 * @brief Non-owning view of the pinned host energy tracker arrays.
 */
struct HostEnergyView {
    double* normal_pot;                     // Array - Nmon * sizeof(<T>)
    double* sliding_pot;                    // Array - Nmon * sizeof(<T>)
    double* rolling_pot;                    // Array - Nmon * sizeof(<T>)
    double* twisting_pot;                   // Array - Nmon * sizeof(<T>)
    double* normal_damp;                    // Array - Nmon * sizeof(<T>)
    double* sliding_inelastic;              // Array - Nmon * sizeof(<T>)
    double* rolling_inelastic;              // Array - Nmon * sizeof(<T>)
    double* twisting_inelastic;             // Array - Nmon * sizeof(<T>)
    double* normal_break;                   // Array - Nmon * sizeof(<T>)
    double* sliding_break;                  // Array - Nmon * sizeof(<T>)
    double* rolling_break;                  // Array - Nmon * sizeof(<T>)
    double* twisting_break;                 // Array - Nmon * sizeof(<T>)
    double* normal_form;                    // Array - Nmon * sizeof(<T>)
};

/**
 * @brief RAII owner of the pinned host memory of the energy trackers.
 *
 * Allocates all arrays on construction and frees them on destruction.
 * Call view() to obtain a HostEnergyView of raw pointers.
 */
class HostEnergy {
    HostBuffer<double> normal_pot;
    HostBuffer<double> sliding_pot;
    HostBuffer<double> rolling_pot;
    HostBuffer<double> twisting_pot;
    HostBuffer<double> normal_damp;
    HostBuffer<double> sliding_inelastic;
    HostBuffer<double> rolling_inelastic;
    HostBuffer<double> twisting_inelastic;
    HostBuffer<double> normal_break;
    HostBuffer<double> sliding_break;
    HostBuffer<double> rolling_break;
    HostBuffer<double> twisting_break;
    HostBuffer<double> normal_form;

public:
    explicit HostEnergy(size_t Nmon)
        : normal_pot(Nmon),        sliding_pot(Nmon),       rolling_pot(Nmon),   twisting_pot(Nmon)
        , normal_damp(Nmon)
        , sliding_inelastic(Nmon), rolling_inelastic(Nmon), twisting_inelastic(Nmon)
        , normal_break(Nmon),      sliding_break(Nmon),     rolling_break(Nmon), twisting_break(Nmon)
        , normal_form(Nmon)
    {}

    /** @brief Generates a view over the pinned host energy trackers for direct data access. */
    HostEnergyView view() {
        return {
            normal_pot.data(),        sliding_pot.data(),       rolling_pot.data(),   twisting_pot.data(),
            normal_damp.data(),
            sliding_inelastic.data(), rolling_inelastic.data(), twisting_inelastic.data(),
            normal_break.data(),      sliding_break.data(),     rolling_break.data(), twisting_break.data(),
            normal_form.data()
        };
    }

    /**
     * @brief Copies the energy trackers from device memory into this host memory.
     *
     * @param d The device energy trackers to pull from.
     * @param Nmon The number of monomers in the system.
     */
    void pull_from(DeviceEnergy& d, size_t Nmon) {
        DeviceEnergyView dv = d.view();
        CHECK_CUDA(cudaMemcpy(normal_pot.data(),         dv.normal_pot,         Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(sliding_pot.data(),        dv.sliding_pot,        Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(rolling_pot.data(),        dv.rolling_pot,        Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(twisting_pot.data(),       dv.twisting_pot,       Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(normal_damp.data(),        dv.normal_damp,        Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(sliding_inelastic.data(),  dv.sliding_inelastic,  Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(rolling_inelastic.data(),  dv.rolling_inelastic,  Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(twisting_inelastic.data(), dv.twisting_inelastic, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(normal_break.data(),       dv.normal_break,       Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(sliding_break.data(),      dv.sliding_break,      Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(rolling_break.data(),      dv.rolling_break,      Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(twisting_break.data(),     dv.twisting_break,     Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(normal_form.data(),        dv.normal_form,        Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    }

    /**
     * @brief Copies the energy trackers into one record per monomer.
     *
     * @param Nmon The number of monomers in the system.
     * @param out  Storage for Nmon records, record i receives the energy trackers of monomer i.
     */
    void records(size_t Nmon, EnergyRecord* out) const {
        for (size_t i = 0; i < Nmon; i++) {
            EnergyRecord& r = out[i];
            r.normal_pot         = normal_pot.data()[i];
            r.sliding_pot        = sliding_pot.data()[i];
            r.rolling_pot        = rolling_pot.data()[i];
            r.twisting_pot       = twisting_pot.data()[i];
            r.normal_damp        = normal_damp.data()[i];
            r.sliding_inelastic  = sliding_inelastic.data()[i];
            r.rolling_inelastic  = rolling_inelastic.data()[i];
            r.twisting_inelastic = twisting_inelastic.data()[i];
            r.normal_break       = normal_break.data()[i];
            r.sliding_break      = sliding_break.data()[i];
            r.rolling_break      = rolling_break.data()[i];
            r.twisting_break     = twisting_break.data()[i];
            r.normal_form        = normal_form.data()[i];
        }
    }
};
