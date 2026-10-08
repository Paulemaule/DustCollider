/**
 * @file energy.cuh
 * @brief Containers for the energy trackers: the potential and dissipated energy split by dofs and dissipation channels.
 * 
 * The energies are tracked by monomer. A pair kernel thread (i,j) adds half of the energy of the pair into the slot of monomer i, 
 * its partner thread (j,i) adds the other half into the slot of j. The host sums the slots over the monomers when a snapshot is saved and then resets them.
 *
 * The trackers are:
 *  - *_pot:   The potential energies stored in the contacts, evaluated on the stored state at every snapshot.
 *  - *_damp:  Energy dissipated through the viscous damping of the normal motion.
 *  - *_slip:  When a contact exceeds its critical displacement in tangential motion, the contact moves inelastically, dissipating energy.
 *  - *_break: Contact breaking, all energy stored in the contact (in all four dofs) is lost.
 *  - *_form:  When a contact forms the potential jumps from 0 to U_N(delta) < 0. This is booked here.
 */

#pragma once

#include "../utils/errors.cuh"
#include "../utils/buffer.cuh"

/**
 * @brief The energy trackers summed over all monomers, one record per snapshot.
 */
struct EnergyRecord {
    double normal_pot;
    double sliding_pot;
    double rolling_pot;
    double twisting_pot;
    double normal_damp;
    double sliding_slip;
    double rolling_slip;
    double twisting_slip;
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
    double* sliding_slip;                   // Array - Nmon * sizeof(<T>)
    double* rolling_slip;                   // Array - Nmon * sizeof(<T>)
    double* twisting_slip;                  // Array - Nmon * sizeof(<T>)
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
    DeviceBuffer<double> sliding_slip;
    DeviceBuffer<double> rolling_slip;
    DeviceBuffer<double> twisting_slip;
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
        : normal_pot(Nmon),   sliding_pot(Nmon),   rolling_pot(Nmon),   twisting_pot(Nmon)
        , normal_damp(Nmon)
        , sliding_slip(Nmon), rolling_slip(Nmon),  twisting_slip(Nmon)
        , normal_break(Nmon), sliding_break(Nmon), rolling_break(Nmon), twisting_break(Nmon)
        , normal_form(Nmon)
    {}

    /** @brief Generates a view over the device energy trackers for passing into kernels. */
    DeviceEnergyView view() {
        return {
            normal_pot.data(),   sliding_pot.data(),   rolling_pot.data(),   twisting_pot.data(),
            normal_damp.data(),
            sliding_slip.data(), rolling_slip.data(),  twisting_slip.data(),
            normal_break.data(), sliding_break.data(), rolling_break.data(), twisting_break.data(),
            normal_form.data()
        };
    }

    /** @brief Sets all energy trackers to zero. */
    void zero() {
        CHECK_CUDA(cudaMemset(normal_pot.data(),     0, normal_pot.size_bytes()));
        CHECK_CUDA(cudaMemset(sliding_pot.data(),    0, sliding_pot.size_bytes()));
        CHECK_CUDA(cudaMemset(rolling_pot.data(),    0, rolling_pot.size_bytes()));
        CHECK_CUDA(cudaMemset(twisting_pot.data(),   0, twisting_pot.size_bytes()));
        CHECK_CUDA(cudaMemset(normal_damp.data(),    0, normal_damp.size_bytes()));
        CHECK_CUDA(cudaMemset(sliding_slip.data(),   0, sliding_slip.size_bytes()));
        CHECK_CUDA(cudaMemset(rolling_slip.data(),   0, rolling_slip.size_bytes()));
        CHECK_CUDA(cudaMemset(twisting_slip.data(),  0, twisting_slip.size_bytes()));
        CHECK_CUDA(cudaMemset(normal_break.data(),   0, normal_break.size_bytes()));
        CHECK_CUDA(cudaMemset(sliding_break.data(),  0, sliding_break.size_bytes()));
        CHECK_CUDA(cudaMemset(rolling_break.data(),  0, rolling_break.size_bytes()));
        CHECK_CUDA(cudaMemset(twisting_break.data(), 0, twisting_break.size_bytes()));
        CHECK_CUDA(cudaMemset(normal_form.data(),    0, normal_form.size_bytes()));
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
    double* sliding_slip;                   // Array - Nmon * sizeof(<T>)
    double* rolling_slip;                   // Array - Nmon * sizeof(<T>)
    double* twisting_slip;                  // Array - Nmon * sizeof(<T>)
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
    HostBuffer<double> sliding_slip;
    HostBuffer<double> rolling_slip;
    HostBuffer<double> twisting_slip;
    HostBuffer<double> normal_break;
    HostBuffer<double> sliding_break;
    HostBuffer<double> rolling_break;
    HostBuffer<double> twisting_break;
    HostBuffer<double> normal_form;

public:
    explicit HostEnergy(size_t Nmon)
        : normal_pot(Nmon),   sliding_pot(Nmon),   rolling_pot(Nmon),   twisting_pot(Nmon)
        , normal_damp(Nmon)
        , sliding_slip(Nmon), rolling_slip(Nmon),  twisting_slip(Nmon)
        , normal_break(Nmon), sliding_break(Nmon), rolling_break(Nmon), twisting_break(Nmon)
        , normal_form(Nmon)
    {}

    /** @brief Generates a view over the pinned host energy trackers for direct data access. */
    HostEnergyView view() {
        return {
            normal_pot.data(),   sliding_pot.data(),   rolling_pot.data(),   twisting_pot.data(),
            normal_damp.data(),
            sliding_slip.data(), rolling_slip.data(),  twisting_slip.data(),
            normal_break.data(), sliding_break.data(), rolling_break.data(), twisting_break.data(),
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
        CHECK_CUDA(cudaMemcpy(normal_pot.data(),     dv.normal_pot,     Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(sliding_pot.data(),    dv.sliding_pot,    Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(rolling_pot.data(),    dv.rolling_pot,    Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(twisting_pot.data(),   dv.twisting_pot,   Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(normal_damp.data(),    dv.normal_damp,    Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(sliding_slip.data(),   dv.sliding_slip,   Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(rolling_slip.data(),   dv.rolling_slip,   Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(twisting_slip.data(),  dv.twisting_slip,  Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(normal_break.data(),   dv.normal_break,   Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(sliding_break.data(),  dv.sliding_break,  Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(rolling_break.data(),  dv.rolling_break,  Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(twisting_break.data(), dv.twisting_break, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(normal_form.data(),    dv.normal_form,    Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    }

    /**
     * @brief Sums the energy trackers over the monomers.
     *
     * @param Nmon The number of monomers in the system.
     * @return The summed energy trackers.
     */
    EnergyRecord totals(size_t Nmon) const {
        auto sum = [Nmon](const HostBuffer<double>& b) {
            double s = 0.0;
            for (size_t i = 0; i < Nmon; i++) s += b.data()[i];
            return s;
        };

        EnergyRecord r;
        r.normal_pot     = sum(normal_pot);
        r.sliding_pot    = sum(sliding_pot);
        r.rolling_pot    = sum(rolling_pot);
        r.twisting_pot   = sum(twisting_pot);
        r.normal_damp    = sum(normal_damp);
        r.sliding_slip   = sum(sliding_slip);
        r.rolling_slip   = sum(rolling_slip);
        r.twisting_slip  = sum(twisting_slip);
        r.normal_break   = sum(normal_break);
        r.sliding_break  = sum(sliding_break);
        r.rolling_break  = sum(rolling_break);
        r.twisting_break = sum(twisting_break);
        r.normal_form    = sum(normal_form);
        return r;
    }
};
