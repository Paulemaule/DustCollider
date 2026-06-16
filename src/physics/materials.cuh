#pragma once

#include "../utils/errors.cuh"
#include "../utils/buffer.cuh"

// TODO: Some material properties are still missing.
// TODO: It might be beneficial to precompute several pairwise quantities that are used often like E*, gamma_ij, ...

/**
 * @brief Non-owning view of device material arrays. Passed directly to kernels.
 */
struct DeviceMaterialsView {
    double* radius;
    double* mass;
    double* moment;
    int*    matID;
    double* density;
    double* surface_energy;
    double* youngs_modulus;
    double* poisson_number;
    double* damping_timescale;
    double* crit_rolling_disp;
};

/**
 * @brief RAII container for the device memory of the monomer material properties.
 *
 * Automatically allocates all arrays on construction and frees them on destruction.
 * To pass the data to the kernel use ::view() to obtain a struct of raw pointers.
 */
class DeviceMaterials {
    DeviceBuffer<double> radius;
    DeviceBuffer<double> mass;
    DeviceBuffer<double> moment;
    DeviceBuffer<int>    matID;
    DeviceBuffer<double> density;
    DeviceBuffer<double> surface_energy;
    DeviceBuffer<double> youngs_modulus;
    DeviceBuffer<double> poisson_number;
    DeviceBuffer<double> damping_timescale;
    DeviceBuffer<double> crit_rolling_disp;

public:
    /**
     * @brief Allocates device memory for the monomer materials of a system with Nmon monomers.
     * @param Nmon The number of monomers in the system.
     */
    explicit DeviceMaterials(size_t Nmon)
        : radius(Nmon),            mass(Nmon),              moment(Nmon)
        , matID(Nmon),             density(Nmon),           surface_energy(Nmon)
        , youngs_modulus(Nmon),    poisson_number(Nmon)
        , damping_timescale(Nmon), crit_rolling_disp(Nmon)
    {}

    /**
     * @brief Generates a view over the device material properties for passing into kernels.
     */
    DeviceMaterialsView view() {
        return {
            radius.data(),            mass.data(),              moment.data(),
            matID.data(),             density.data(),           surface_energy.data(),
            youngs_modulus.data(),    poisson_number.data(),
            damping_timescale.data(), crit_rolling_disp.data()
        };
    }
};

/**
 * @brief Non-owning view of pinned host material arrays.
 */
struct HostMaterialsView {
    double* radius;
    double* mass;
    double* moment;
    int*    matID;
    double* density;
    double* surface_energy;
    double* youngs_modulus;
    double* poisson_number;
    double* damping_timescale;
    double* crit_rolling_disp;
};

/**
 * @brief RAII owner of pinned host material arrays.
 *
 * Allocates all arrays on construction and frees them on destruction.
 * Call view() to obtain a HostMaterialsView of raw pointers.
 */
class HostMaterials {
    HostBuffer<double> radius;
    HostBuffer<double> mass;
    HostBuffer<double> moment;
    HostBuffer<int>    matID;
    HostBuffer<double> density;
    HostBuffer<double> surface_energy;
    HostBuffer<double> youngs_modulus;
    HostBuffer<double> poisson_number;
    HostBuffer<double> damping_timescale;
    HostBuffer<double> crit_rolling_disp;

public:
    explicit HostMaterials(size_t Nmon)
        : radius(Nmon),            mass(Nmon),              moment(Nmon)
        , matID(Nmon),             density(Nmon),           surface_energy(Nmon)
        , youngs_modulus(Nmon),    poisson_number(Nmon)
        , damping_timescale(Nmon), crit_rolling_disp(Nmon)
    {}

    /** @brief Returns a view of raw pointers into the pinned host material arrays. */
    HostMaterialsView view() {
        return {
            radius.data(),            mass.data(),              moment.data(),
            matID.data(),             density.data(),           surface_energy.data(),
            youngs_modulus.data(),    poisson_number.data(),
            damping_timescale.data(), crit_rolling_disp.data()
        };
    }

    /** @brief Copies all material arrays to device memory. */
    void push_to(DeviceMaterials& d, size_t Nmon) {
        DeviceMaterialsView dv = d.view();
        CHECK_CUDA(cudaMemcpy(dv.radius,             radius.data(),            Nmon * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.mass,               mass.data(),              Nmon * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.moment,             moment.data(),            Nmon * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.matID,              matID.data(),             Nmon * sizeof(int),    cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.density,            density.data(),           Nmon * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.surface_energy,     surface_energy.data(),    Nmon * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.youngs_modulus,     youngs_modulus.data(),    Nmon * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.poisson_number,     poisson_number.data(),    Nmon * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.damping_timescale,  damping_timescale.data(), Nmon * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.crit_rolling_disp,  crit_rolling_disp.data(), Nmon * sizeof(double), cudaMemcpyHostToDevice));
    }

    /** @brief Copies all material arrays from device memory into this host struct. */
    void pull_from(DeviceMaterials& d, size_t Nmon) {
        DeviceMaterialsView dv = d.view();
        CHECK_CUDA(cudaMemcpy(radius.data(),            dv.radius,             Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(mass.data(),              dv.mass,               Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(moment.data(),            dv.moment,             Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(matID.data(),             dv.matID,              Nmon * sizeof(int),    cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(density.data(),           dv.density,            Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(surface_energy.data(),    dv.surface_energy,     Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(youngs_modulus.data(),    dv.youngs_modulus,     Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(poisson_number.data(),    dv.poisson_number,     Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(damping_timescale.data(), dv.damping_timescale,  Nmon * sizeof(double), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(crit_rolling_disp.data(), dv.crit_rolling_disp,  Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    }
};