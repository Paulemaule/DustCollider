#pragma once

#include <algorithm>
#include <cstring>
#include "../utils/buffer.cuh"


/**
 * @brief Non-owning view of device state arrays. Passed directly to kernels.
 */
struct DeviceStateView {
    double3*    position;                      // Array - Nmon * sizeof(<T>)
    double3*    velocity;                      // Array - Nmon * sizeof(<T>)
    double3*    omega;                         // Array - Nmon * sizeof(<T>)
    double3*    force;                         // Array - Nmon * sizeof(<T>)
    double3*    torque;                        // Array - Nmon * sizeof(<T>)
    double*     contact_compression;           // Matrix - Nmon * Nmon * sizeof(<T>)
    double*     contact_twist;                 // Matrix - Nmon * Nmon * sizeof(<T>)
    double3*    contact_pointer;               // Matrix - Nmon * Nmon * sizeof(<T>)
    double3*    contact_normal;                // Matrix - Nmon * Nmon * sizeof(<T>)
    double4*    contact_rotation;              // Matrix - Nmon * Nmon * sizeof(<T>)
};

/**
 * @brief RAII container for the device memory of the system state.
 *
 * Automatically allocates all arrays on construction and frees them on destruction.
 * To pass the data to the kernel use ::view() to obtain a struct of raw pointers.
 */
class DeviceState {
    DeviceBuffer<double3> position;
    DeviceBuffer<double3> velocity;
    DeviceBuffer<double3> omega;
    DeviceBuffer<double3> force;
    DeviceBuffer<double3> torque;
    DeviceBuffer<double>  contact_compression;
    DeviceBuffer<double>  contact_twist;
    DeviceBuffer<double3> contact_pointer;
    DeviceBuffer<double3> contact_normal;
    DeviceBuffer<double4> contact_rotation;

public:
    /**
     * @brief Allocates device memory for the state of a system with Nmon monomers.
     * @param Nmon The number of monomers in the system.
     */
    explicit DeviceState(size_t Nmon) 
        : position(Nmon),            velocity(Nmon),           omega(Nmon)
        , force(Nmon),               torque(Nmon)
        , contact_compression(Nmon * Nmon), contact_twist(Nmon * Nmon)
        , contact_pointer(Nmon * Nmon),     contact_normal(Nmon * Nmon)
        , contact_rotation(Nmon * Nmon)
    {}

    /** @brief Generates a view over the device system state for passing into kernels. */
    DeviceStateView view() {
        return {
            position.data(),            velocity.data(),          omega.data(),
            force.data(),               torque.data(),
            contact_compression.data(), contact_twist.data(),
            contact_pointer.data(),     contact_normal.data(),
            contact_rotation.data()
        };
    }


};

/**
 * @brief Non-owning view of pinned host state arrays.
 */
struct HostStateView {
    double3*    position;                      // Array - Nmon * sizeof(<T>)
    double3*    velocity;                      // Array - Nmon * sizeof(<T>)
    double3*    omega;                         // Array - Nmon * sizeof(<T>)
    double3*    force;                         // Array - Nmon * sizeof(<T>)
    double3*    torque;                        // Array - Nmon * sizeof(<T>)
    double*     contact_compression;           // Matrix - Nmon * Nmon * sizeof(<T>)
    double*     contact_twist;                 // Matrix - Nmon * Nmon * sizeof(<T>)
    double3*    contact_pointer;               // Matrix - Nmon * Nmon * sizeof(<T>)
    double3*    contact_normal;                // Matrix - Nmon * Nmon * sizeof(<T>)
    double4*    contact_rotation;              // Matrix - Nmon * Nmon * sizeof(<T>)
};

/**
 * @brief RAII owner of all pinned host-side state arrays for one time-step buffer.
 *
 * Allocates all arrays on construction and frees them on destruction.
 * Call view() to obtain a HostStateView of raw pointers.
 */
class HostState {
    HostBuffer<double3> position;
    HostBuffer<double3> velocity;
    HostBuffer<double3> omega;
    HostBuffer<double3> force;
    HostBuffer<double3> torque;
    HostBuffer<double>  contact_compression;
    HostBuffer<double>  contact_twist;
    HostBuffer<double3> contact_pointer;
    HostBuffer<double3> contact_normal;
    HostBuffer<double4> contact_rotation;

public:
    explicit HostState(size_t Nmon)
        : position(Nmon),            velocity(Nmon),           omega(Nmon)
        , force(Nmon),               torque(Nmon)
        , contact_compression(Nmon * Nmon), contact_twist(Nmon * Nmon)
        , contact_pointer(Nmon * Nmon),     contact_normal(Nmon * Nmon)
        , contact_rotation(Nmon * Nmon)
    {}

    /** @brief Generates a view over the pinned host system state for direct data access. */
    HostStateView view() {
        return {
            position.data(),            velocity.data(),          omega.data(),
            force.data(),               torque.data(),
            contact_compression.data(), contact_twist.data(),
            contact_pointer.data(),     contact_normal.data(),
            contact_rotation.data()
        };
    }

    /**
     * @brief Copies this host memory system state into a device side system state.
     *  
     * @param d The device state to push the state into.
     * @param Nmon The number of monomers in the system.
     */
    void push_to(DeviceState& d, size_t Nmon) {
        DeviceStateView dv = d.view();
        CHECK_CUDA(cudaMemcpy(dv.position,            position.data(),            Nmon        * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.velocity,            velocity.data(),            Nmon        * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.omega,               omega.data(),               Nmon        * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.force,               force.data(),               Nmon        * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.torque,              torque.data(),              Nmon        * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.contact_compression, contact_compression.data(), Nmon * Nmon * sizeof(double),  cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.contact_twist,       contact_twist.data(),       Nmon * Nmon * sizeof(double),  cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.contact_pointer,     contact_pointer.data(),     Nmon * Nmon * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.contact_normal,      contact_normal.data(),      Nmon * Nmon * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(dv.contact_rotation,    contact_rotation.data(),    Nmon * Nmon * sizeof(double4), cudaMemcpyHostToDevice));
    }

    /** 
     * @brief Copies the system state from device memory into this host memory.
     *  
     * @param d The device state to pull the state from.
     * @param Nmon The number of monomers in the system.
     */
    void pull_from(DeviceState& d, size_t Nmon) {
        DeviceStateView dv = d.view();
        CHECK_CUDA(cudaMemcpy(position.data(),            dv.position,            Nmon        * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(velocity.data(),            dv.velocity,            Nmon        * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(omega.data(),               dv.omega,               Nmon        * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(force.data(),               dv.force,               Nmon        * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(torque.data(),              dv.torque,              Nmon        * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(contact_compression.data(), dv.contact_compression, Nmon * Nmon * sizeof(double),  cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(contact_twist.data(),       dv.contact_twist,       Nmon * Nmon * sizeof(double),  cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(contact_pointer.data(),     dv.contact_pointer,     Nmon * Nmon * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(contact_normal.data(),      dv.contact_normal,      Nmon * Nmon * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(contact_rotation.data(),    dv.contact_rotation,    Nmon * Nmon * sizeof(double4), cudaMemcpyDeviceToHost));
    }
};