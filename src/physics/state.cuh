#pragma once

#include <iostream>

#include "/utils/printing.cuh"
#include "/utils/errors.cuh"

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
    CHECK_CUDA(cudaMallocHost(& state.position, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.magnetization, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.velocity, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.omega, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.magnetization_change, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.force, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.torque, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.contact_compression, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.contact_twist, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.contact_pointer, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.contact_normal, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMallocHost(& state.contact_rotation, Nmon * Nmon * sizeof(double3)));
}

/**
 * @brief Allocates memory for a device state.
 * 
 * @param state: The device state for whos members to allocate memory.
 */
void state_allocateDeviceMemory(deviceState& state, size_t Nmon) {
    // Allocate memory of the appropriate size on the device for all fields of the state.
    CHECK_CUDA(cudaMalloc(& state.position, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.magnetization, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.velocity, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.omega, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.magnetization_change, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.force, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.torque, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.contact_compression, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.contact_twist, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.contact_pointer, Nmon * Nmon * sizeof(double3)))
    CHECK_CUDA(cudaMalloc(& state.contact_normal, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(& state.contact_rotation, Nmon * Nmon * sizeof(double3)));
}

/**
 * @brief Frees all allocated memory contained by the host state.
 * 
 * @param The state contining the pointers that are to be freed.
 */
void state_freeHost(hostState& state) {
    // Free all arrays pointed to by members of the state.
    CHECK_CUDA(cudaFreeHost(state.position));
    CHECK_CUDA(cudaFreeHost(state.magnetization));
    CHECK_CUDA(cudaFreeHost(state.velocity));
    CHECK_CUDA(cudaFreeHost(state.omega));
    CHECK_CUDA(cudaFreeHost(state.magnetization_change));
    CHECK_CUDA(cudaFreeHost(state.force));
    CHECK_CUDA(cudaFreeHost(state.torque));
    CHECK_CUDA(cudaFreeHost(state.contact_compression));
    CHECK_CUDA(cudaFreeHost(state.contact_twist));
    CHECK_CUDA(cudaFreeHost(state.contact_pointer));
    CHECK_CUDA(cudaFreeHost(state.contact_normal));
    CHECK_CUDA(cudaFreeHost(state.contact_rotation));
}

/**
 * @brief Frees all allocated memory contained by the device state.
 * 
 * @param The state containing the pointers that are to be freed.
 */
void state_freeDevice(deviceState& state) {
    // Free all arrays pointed to by members of the state.
    CHECK_CUDA(cudaFree(state.position));
    CHECK_CUDA(cudaFree(state.magnetization));
    CHECK_CUDA(cudaFree(state.velocity));
    CHECK_CUDA(cudaFree(state.omega));
    CHECK_CUDA(cudaFree(state.magnetization_change));
    CHECK_CUDA(cudaFree(state.force));
    CHECK_CUDA(cudaFree(state.torque));
    CHECK_CUDA(cudaFree(state.contact_compression));
    CHECK_CUDA(cudaFree(state.contact_twist));
    CHECK_CUDA(cudaFree(state.contact_pointer));
    CHECK_CUDA(cudaFree(state.contact_normal));
    CHECK_CUDA(cudaFree(state.contact_rotation));
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
    CHECK_CUDA(cudaMemcpy(
        device_state.position, host_state.position, 
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.magnetization, host_state.magnetization, 
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.velocity, host_state.velocity,
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.omega, host_state.omega,
        Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.magnetization_change, host_state.magnetization_change,
        Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.force, host_state.force,
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.torque, host_state.torque, 
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.contact_compression, host_state.contact_compression, 
        Nmon * Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.contact_twist, host_state.contact_twist,
        Nmon * Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.contact_pointer, host_state.contact_pointer,
        Nmon * Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.contact_normal, host_state.contact_normal,
        Nmon * Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(
        device_state.contact_rotation, host_state.contact_rotation,
        Nmon * Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyHostToDevice));
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
    CHECK_CUDA(cudaMemcpy(
        host_state.position, device_state.position,
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.magnetization, device_state.magnetization, 
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.velocity, device_state.velocity,
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.omega, device_state.omega,
        Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyDeviceToHost))
    CHECK_CUDA(cudaMemcpy(
        host_state.magnetization_change, device_state.magnetization_change,
        Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.force, device_state.force,
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.torque, device_state.torque, 
        Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.contact_compression, device_state.contact_compression, 
        Nmon * Nmon * sizeof(double3), 
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.contact_twist, device_state.contact_twist,
        Nmon * Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.contact_pointer, device_state.contact_pointer,
        Nmon * Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.contact_normal, device_state.contact_normal,
        Nmon * Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(
        host_state.contact_rotation, device_state.contact_rotation,
        Nmon * Nmon * sizeof(double3),
        cudaMemcpyKind::cudaMemcpyDeviceToHost));
}

/**
 * @brief Sets the values of all fields of a host state to a specified value.
 * 
 * @param &state: The device state whos values are to be set.
 * @param value: The value the fields are to be set to.
 * @param Nmon: The number of monomers.
 */
void state_setValueHost(hostState& state, double value, size_t Nmon) {
    memset(state.position, value, Nmon * sizeof(double3));
    memset(state.magnetization, value, Nmon * sizeof(double3));
    memset(state.velocity, value, Nmon * sizeof(double3));
    memset(state.omega, value, Nmon * sizeof(double3));
    memset(state.magnetization_change, value, Nmon * sizeof(double3));
    memset(state.force, value, Nmon * sizeof(double3));
    memset(state.torque, value, Nmon * sizeof(double3));
    memset(state.contact_compression, value, Nmon * Nmon * sizeof(double3));
    memset(state.contact_twist, value, Nmon * Nmon * sizeof(double3));
    memset(state.contact_pointer, value, Nmon * Nmon * sizeof(double3));
    memset(state.contact_normal, value, Nmon * Nmon * sizeof(double3));
    memset(state.contact_twist, value, Nmon * Nmon * sizeof(double3)); 
}

/**
 * @brief Sets the values of all fields of a device state to a specified value.
 * 
 * @param &state: The device state whos values are to be set.
 * @param value: The value the fields are to be set to.
 * @param Nmon: The number of monomers.
 */
void state_setValueDevice(deviceState& state, double value, size_t Nmon) {
    CHECK_CUDA(cudaMemset(state.position, value, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.magnetization, value, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.velocity, value, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.omega, value, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.magnetization_change, value, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.force, value, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.torque, value, Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.contact_compression, value, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.contact_twist, value, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.contact_pointer, value, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.contact_normal, value, Nmon * Nmon * sizeof(double3)));
    CHECK_CUDA(cudaMemset(state.contact_twist, value, Nmon * Nmon * sizeof(double3)));
}