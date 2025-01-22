#pragma once

#include <iostream>

#include "../makros/printing.cuh"
#include "../makros/errors.cuh"

/**
 * @brief A struct of pointers to host memory of arrays containing system state information.
 */
typedef struct {
    double3* position;                      // Array - Nmon * sizeof(<T>)
    double3* magnetization;                 // Array - Nmon * sizeof(<T>)
    double3* velocity;                      // Array - Nmon * sizeof(<T>)
    double3* omega;                         // Array - Nmon * sizeof(<T>)
    double3* magnetization_change;          // Array - Nmon * sizeof(<T>)
    double3* force;                         // Array - Nmon * sizeof(<T>)
    double3* torque;                        // Array - Nmon * sizeof(<T>)
    double3* contact_compression;           // Matrix - Nmon * Nmon * sizeof(<T>)
    double3* contact_twist;                 // Matrix - Nmon * Nmon * sizeof(<T>)
    double3* contact_pointer;               // Matrix - Nmon * Nmon * sizeof(<T>)
    double3* contact_normal;                // Matrix - Nmon * Nmon * sizeof(<T>)
    double3* contact_rotation;              // Matrix - Nmon * Nmon * sizeof(<T>)
} hostState;

/**
 * @brief A struct of pointers to device memory of arrays containing system state information.
 */
typedef struct {
    double3* position;                      // Array - Nmon * sizeof(<T>)
    double3* magnetization;                 // Array - Nmon * sizeof(<T>)
    double3* velocity;                      // Array - Nmon * sizeof(<T>)
    double3* omega;                         // Array - Nmon * sizeof(<T>)
    double3* magnetization_change;          // Array - Nmon * sizeof(<T>)
    double3* force;                         // Array - Nmon * sizeof(<T>)
    double3* torque;                        // Array - Nmon * sizeof(<T>)
    double3* contact_compression;           // Matrix - Nmon * Nmon * sizeof(<T>)
    double3* contact_twist;                 // Matrix - Nmon * Nmon * sizeof(<T>)
    double3* contact_pointer;               // Matrix - Nmon * Nmon * sizeof(<T>)
    double3* contact_normal;                // Matrix - Nmon * Nmon * sizeof(<T>)
    double3* contact_rotation;              // Matrix - Nmon * Nmon * sizeof(<T>)
} deviceState;

/**
 * @brief Allocates memory for a host state.
 * 
 * @param state: The host state for whos members to allocate memory.
 */
void state_allocateHostMemory(hostState& state, size_t Nmon) {
    // Allocate memory of the appropriate size for all fields of the state.
    // This memory is pinned memory!
    if (cudaMallocHost(& state.position, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for positions.");
    }
    if (cudaMallocHost(& state.magnetization, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for magnetization.");
    }
    if (cudaMallocHost(& state.velocity, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for velocities.");
    }
    if (cudaMallocHost(& state.omega, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for angular velocity.");
    }
    if (cudaMallocHost(& state.magnetization_change, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for change in magnetization.");
    }
    if (cudaMallocHost(& state.force, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for force.");
    }
    if (cudaMallocHost(& state.torque, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for torque.");
    }
    if (cudaMallocHost(& state.contact_compression, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for contact compression lenght.");
    }
    if (cudaMallocHost(& state.contact_twist, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for contact twisting angle.");
    }
    if (cudaMallocHost(& state.contact_pointer, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for contact pointer.");
    }
    if (cudaMallocHost(& state.contact_normal, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for contact normal.");
    }
    if (cudaMallocHost(& state.contact_rotation, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on host for contact pointer rotation.");
    }

    // Check if any CUDA errors occured.
    CUDA_LAST_ERROR_CHECK();
}

/**
 * @brief Allocates memory for a device state.
 * 
 * @param state: The device state for whos members to allocate memory.
 */
void state_allocateDeviceMemory(deviceState& state, size_t Nmon) {
    // Allocate memory of the appropriate size on the device for all fields of the state.
    if (cudaMalloc(& state.position, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for positions.");
    }
    if (cudaMalloc(& state.magnetization, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for magnetization.");
    }
    if (cudaMalloc(& state.velocity, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for velocities.");
    }
    if (cudaMalloc(& state.omega, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for angular velocity.");
    }
    if (cudaMalloc(& state.magnetization_change, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for change in magnetization.");
    }
    if (cudaMalloc(& state.force, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for force.");
    }
    if (cudaMalloc(& state.torque, Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for torque.");
    }
    if (cudaMalloc(& state.contact_compression, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for contact compression lenght.");
    }
    if (cudaMalloc(& state.contact_twist, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for contact twisting angle.");
    }
    if (cudaMalloc(& state.contact_pointer, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for contact pointer.");
    }
    if (cudaMalloc(& state.contact_normal, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for contact normal.");
    }
    if (cudaMalloc(& state.contact_rotation, Nmon * Nmon * sizeof(double3))) {
        PRINT_ERROR("Failed to allocate memory on device for contact pointer rotation.");
    }

    // Check if there were any CUDA errors.
    CUDA_LAST_ERROR_CHECK();
}

/**
 * @brief Copies the fields of a host state to a device state.
 * 
 * @param &host_state: The host state from which to copy.
 * @param &device_state: The device state to which to copy.
 * @param Nmon: The number of monomers.
 */
void state_pushToDevice(hostState& host_state, deviceState& device_state, size_t Nmon) {
    // Copy the contents of the arrays pointed to by the pointers to the device.
    if (cudaMemcpy(device_state.position, host_state.position, 
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push positions to device.");
    }
    if (cudaMemcpy(device_state.magnetization, host_state.magnetization, 
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push magnetization to device.");
    }
    if (cudaMemcpy(device_state.velocity, host_state.velocity,
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push velocity to device.");
    }
    if (cudaMemcpy(device_state.omega, host_state.omega,
            Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push omega to device.");
    }
    if (cudaMemcpy(device_state.magnetization_change, host_state.magnetization_change,
            Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push change in magnetization to device.");
    }
    if (cudaMemcpy(device_state.force, host_state.force,
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push force to device.");
    }
    if (cudaMemcpy(device_state.torque, host_state.torque, 
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push torque to device.");
    }
    if (cudaMemcpy(device_state.contact_compression, host_state.contact_compression, 
            Nmon * Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push  to device.");
    }
    if (cudaMemcpy(device_state.contact_twist, host_state.contact_twist,
            Nmon * Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push  to device.");
    }
    if (cudaMemcpy(device_state.contact_pointer, host_state.contact_pointer,
            Nmon * Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push  to device.");
    }
    if (cudaMemcpy(device_state.contact_normal, host_state.contact_normal,
            Nmon * Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push  to device.");
    }
    if (cudaMemcpy(device_state.contact_rotation, host_state.contact_rotation,
            Nmon * Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyHostToDevice)) {
        PRINT_ERROR("Could not push  to device.");
    }

    // Check if there were any CUDA errors.
    CUDA_LAST_ERROR_CHECK();
}

/**
 * @brief Copies the fields of a device state to a host state.
 * 
 * @param &device_state: The device state from which to copy.
 * @param &host_state: The host state to which to copy.
 * @param Nmon: The number of monomers.
 */
void state_pullFromDevice(deviceState& device_state, hostState& host_state, size_t Nmon) {
    // Copy the contents of the arrays pointed to by the pointers to the host.
    if (cudaMemcpy(host_state.position, device_state.position,
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push positions to device.");
    }
    if (cudaMemcpy(host_state.magnetization, device_state.magnetization, 
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push magnetization to device.");
    }
    if (cudaMemcpy(host_state.velocity, device_state.velocity,
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push velocity to device.");
    }
    if (cudaMemcpy(host_state.omega, device_state.omega,
            Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push omega to device.");
    }
    if (cudaMemcpy(host_state.magnetization_change, device_state.magnetization_change,
            Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push change in magnetization to device.");
    }
    if (cudaMemcpy(host_state.force, device_state.force,
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push force to device.");
    }
    if (cudaMemcpy(host_state.torque, device_state.torque, 
            Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push torque to device.");
    }
    if (cudaMemcpy(host_state.contact_compression, device_state.contact_compression, 
            Nmon * Nmon * sizeof(double3), 
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push  to device.");
    }
    if (cudaMemcpy(host_state.contact_twist, device_state.contact_twist,
            Nmon * Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push  to device.");
    }
    if (cudaMemcpy(host_state.contact_pointer, device_state.contact_pointer,
            Nmon * Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push  to device.");
    }
    if (cudaMemcpy(host_state.contact_normal, device_state.contact_normal,
            Nmon * Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push  to device.");
    }
    if (cudaMemcpy(host_state.contact_rotation, device_state.contact_rotation,
            Nmon * Nmon * sizeof(double3),
            cudaMemcpyKind::cudaMemcpyDeviceToHost)) {
        PRINT_ERROR("Could not push  to device.");
    }

    // Check if there were any CUDA errors.
    CUDA_LAST_ERROR_CHECK();
}

/**
 * @brief Frees all allocated memory contained by the host state.
 * 
 * @param The state contining the pointers that are to be freed.
 */
void state_freeHost(hostState& state) {
    // Free all arrays pointed to by members of the state.
    if(cudaFreeHost(state.position)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.magnetization)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.velocity)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.omega)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.magnetization_change)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.force)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.torque)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.contact_compression)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.contact_twist)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.contact_pointer)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.contact_normal)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFreeHost(state.contact_rotation)) {
        PRINT_ERROR("Could not free host memory.")
    }
    
    // Check if there were any CUDA errors.
    CUDA_LAST_ERROR_CHECK();
}

/**
 * @brief Frees all allocated memory contained by the device state.
 * 
 * @param The state contining the pointers that are to be freed.
 */
void state_freeDevice(deviceState& state) {
    // Free all arrays pointed to by members of the state.
    if(cudaFree(state.position)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.magnetization)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.velocity)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.omega)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.magnetization_change)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.force)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.torque)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.contact_compression)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.contact_twist)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.contact_pointer)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.contact_normal)) {
        PRINT_ERROR("Could not free host memory.")
    }
    if(cudaFree(state.contact_rotation)) {
        PRINT_ERROR("Could not free host memory.")
    }
    
    // Check if there were any CUDA errors.
    CUDA_LAST_ERROR_CHECK();
}
