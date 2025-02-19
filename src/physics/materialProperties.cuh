#pragma once

#include "../utils/errors.cuh"

// TODO: Some material properties are still missing.
// TODO: It might be beneficial to precompute several pairwise quantities that are used often like E*, gamma_ij, ...

/**
 * @brief A struct of pointers to host memory of arrays containing static monomer properties.
 */
typedef struct {
    double* radius;                     // Array - Nmon * sizeof(double)
    double* mass;                       // Array - Nmon * sizeof(double)
    double* moment;                     // Array - Nmon * sizeof(double)
    int* matID;                         // Array - Nmon * sizeof(int)
    double* density;                    // Array - Nmon * sizeof(double)
    double* surface_energy;             // Array - Nmon * sizeof(double)
    double* youngs_modulus;             // Array - Nmon * sizeof(double)
    double* poisson_number;             // Array - Nmon * sizeof(double)
    double* damping_timescale;          // Array - Nmon * sizeof(double)
    double* crit_rolling_disp;          // Array - Nmon * sizeof(double)
    //double* magnetic_susceptibiliy;     // Array - Nmon * sizeof(double)
} hostMatProperties;

/**
 * @brief A struct of pointers to device memory of arrays containing static monomer properties.
 */
typedef struct {
    double* radius;                     // Array - Nmon * sizeof(double)
    double* mass;                       // Array - Nmon * sizeof(double)
    double* moment;                     // Array - Nmon * sizeof(double)
    int* matID;                         // Array - Nmon * sizeof(int)
    double* density;                    // Array - Nmon * sizeof(double)
    double* surface_energy;             // Array - Nmon * sizeof(double)
    double* youngs_modulus;             // Array - Nmon * sizeof(double)
    double* poisson_number;             // Array - Nmon * sizeof(double)
    double* damping_timescale;          // Array - Nmon * sizeof(double)
    double* crit_rolling_disp;          // Array - Nmon * sizeof(double)
    //double* magnetic_susceptibiliy;     // Array - Nmon * sizeof(double)
} deviceMatProperties;

/**
 * @brief Allocates memory for a host material properties struct.
 * 
 * @param properties: The host material properties struct for which to allocate memory.
 * @param Nmon: The number of monomers.
 */
void matProperties_allocateHostMemory(hostMatProperties& properties, size_t Nmon) {
    CHECK_CUDA(cudaMallocHost(&properties.radius, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&properties.mass, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&properties.moment, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&properties.matID, Nmon * sizeof(int)));
    CHECK_CUDA(cudaMallocHost(&properties.density, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&properties.surface_energy, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&properties.youngs_modulus, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&properties.poisson_number, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&properties.damping_timescale, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&properties.crit_rolling_disp, Nmon * sizeof(double)));
    //CHECK_CUDA(cudaMallocHost(&properties.magnetic_susceptibiliy, Nmon * sizeof(double)));
}

/**
 * @brief Frees all allocated memory contained by the host material properties struct.
 * 
 * @param properties: The host material properties struct containing the pointers to free.
 */
void matProperties_freeHostMemory(hostMatProperties& properties) {
    CHECK_CUDA(cudaFreeHost(properties.radius));
    CHECK_CUDA(cudaFreeHost(properties.mass));
    CHECK_CUDA(cudaFreeHost(properties.moment));
    CHECK_CUDA(cudaFreeHost(properties.matID));
    CHECK_CUDA(cudaFreeHost(properties.density));
    CHECK_CUDA(cudaFreeHost(properties.surface_energy));
    CHECK_CUDA(cudaFreeHost(properties.youngs_modulus));
    CHECK_CUDA(cudaFreeHost(properties.poisson_number));
    CHECK_CUDA(cudaFreeHost(properties.damping_timescale));
    CHECK_CUDA(cudaFreeHost(properties.crit_rolling_disp));
    //CHECK_CUDA(cudaFreeHost(properties.magnetic_susceptibiliy));

    // Set pointers to nullptr to avoid dangling pointers
    properties.radius = nullptr;
    properties.mass = nullptr;
    properties.moment = nullptr;
    properties.matID = nullptr;
    properties.density = nullptr;
    properties.surface_energy = nullptr;
    properties.youngs_modulus = nullptr;
    properties.poisson_number = nullptr;
    properties.damping_timescale = nullptr;
    properties.crit_rolling_disp = nullptr;
    //properties.magnetic_susceptibiliy = nullptr;
}

/**
 * @brief Allocates memory for a device material properties struct.
 * 
 * @param properties: The device material properties struct for which to allocate memory.
 * @param Nmon: The number of monomers.
 */
void matProperties_allocateDeviceMemory(deviceMatProperties& properties, size_t Nmon) {
    CHECK_CUDA(cudaMalloc(&properties.radius, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&properties.mass, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&properties.moment, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&properties.matID, Nmon * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&properties.density, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&properties.surface_energy, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&properties.youngs_modulus, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&properties.poisson_number, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&properties.damping_timescale, Nmon * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&properties.crit_rolling_disp, Nmon * sizeof(double)));
    //CHECK_CUDA(cudaMalloc(&properties.magnetic_susceptibiliy, Nmon * sizeof(double)));
}

/**
 * @brief Frees all allocated memory contained by the device material properties struct.
 * 
 * @param properties: The device material properties struct containing the pointers to free.
 */
void matProperties_freeDeviceMemory(deviceMatProperties& properties) {
    CHECK_CUDA(cudaFree(properties.radius));
    CHECK_CUDA(cudaFree(properties.mass));
    CHECK_CUDA(cudaFree(properties.moment));
    CHECK_CUDA(cudaFree(properties.matID));
    CHECK_CUDA(cudaFree(properties.density));
    CHECK_CUDA(cudaFree(properties.surface_energy));
    CHECK_CUDA(cudaFree(properties.youngs_modulus));
    CHECK_CUDA(cudaFree(properties.poisson_number));
    CHECK_CUDA(cudaFree(properties.damping_timescale));
    CHECK_CUDA(cudaFree(properties.crit_rolling_disp));
    //CHECK_CUDA(cudaFree(properties.magnetic_susceptibiliy));

    // Set pointers to nullptr to avoid dangling pointers
    properties.radius = nullptr;
    properties.mass = nullptr;
    properties.moment = nullptr;
    properties.matID = nullptr;
    properties.density = nullptr;
    properties.surface_energy = nullptr;
    properties.youngs_modulus = nullptr;
    properties.poisson_number = nullptr;
    properties.damping_timescale = nullptr;
    properties.crit_rolling_disp = nullptr;
    //properties.magnetic_susceptibiliy = nullptr;
}

/**
 * @brief Copies the fields of a host material properties struct to a device material properties struct.
 * 
 * @param host_properties: The host material properties struct from which to copy.
 * @param device_properties: The device material properties struct to which to copy.
 * @param Nmon: The number of monomers.
 */
void matProperties_pushToDevice(const hostMatProperties& host_properties, deviceMatProperties& device_properties, size_t Nmon) {
    CHECK_CUDA(cudaMemcpy(device_properties.radius, host_properties.radius, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.mass, host_properties.mass, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.moment, host_properties.moment, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.matID, host_properties.matID, Nmon * sizeof(int), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.density, host_properties.density, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.surface_energy, host_properties.surface_energy, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.youngs_modulus, host_properties.youngs_modulus, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.poisson_number, host_properties.poisson_number, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.damping_timescale, host_properties.damping_timescale, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(device_properties.crit_rolling_disp, host_properties.damping_timescale, Nmon * sizeof(double), cudaMemcpyHostToDevice));
    //CHECK_CUDA(cudaMemcpy(device_properties.magnetic_susceptibiliy, host_properties.magnetic_susceptibiliy, Nmon * sizeof(double), cudaMemcpyHostToDevice));
}

/**
 * @brief Copies the fields of a device material properties struct to a host material properties struct.
 * 
 * @param device_properties: The device material properties struct from which to copy.
 * @param host_properties: The host material properties struct to which to copy.
 * @param Nmon: The number of monomers.
 */
void matProperties_pullFromDevice(const deviceMatProperties& device_properties, hostMatProperties& host_properties, size_t Nmon) {
    CHECK_CUDA(cudaMemcpy(host_properties.radius, device_properties.radius, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.mass, device_properties.mass, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.moment, device_properties.moment, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.matID, device_properties.matID, Nmon * sizeof(int), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.density, device_properties.density, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.surface_energy, device_properties.surface_energy, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.youngs_modulus, device_properties.youngs_modulus, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.poisson_number, device_properties.poisson_number, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.damping_timescale, device_properties.damping_timescale, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(host_properties.crit_rolling_disp, device_properties.damping_timescale, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
    //CHECK_CUDA(cudaMemcpy(host_properties.magnetic_susceptibiliy, device_properties.magnetic_susceptibiliy, Nmon * sizeof(double), cudaMemcpyDeviceToHost));
}

/**
 * @brief Initializes the fields of a host material properties struct to default values.
 * 
 * @param properties: The host material properties struct to initialize.
 * @param Nmon: The number of monomers.
 */
void matProperties_clearHost(hostMatProperties& properties, size_t Nmon) {
    // FIXME: Memset should not be used for complex types...
    memset(properties.radius, 0, Nmon * sizeof(double));
    memset(properties.mass, 0, Nmon * sizeof(double));
    memset(properties.moment, 0, Nmon * sizeof(double));
    memset(properties.matID, 0, Nmon * sizeof(int));
    memset(properties.density, 0, Nmon * sizeof(double));
    memset(properties.surface_energy, 0, Nmon * sizeof(double));
    memset(properties.youngs_modulus, 0, Nmon * sizeof(double));
    memset(properties.poisson_number, 0, Nmon * sizeof(double));
    memset(properties.damping_timescale, 0, Nmon * sizeof(double));
    memset(properties.crit_rolling_disp, 0, Nmon * sizeof(double));
    //memset(properties.magnetic_susceptibiliy, 0, Nmon * sizeof(double));
}