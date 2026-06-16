#pragma once

#include <utility>
#include <cuda_runtime.h>
#include "errors.cuh"

/**
 * @brief RAII wrapper for a device memory allocation of N elements of type T.
 *
 * Allocates device memory on construction and automatically frees it on destruction.
 * Ownership of the memory is unique, copying is forbidden and moving transfers ownership.
 *
 * @tparam T The element type of the data in the buffer.
 */
template <typename T>
class DeviceBuffer {
    // Raw pointer to the start of the device memory.
    T* ptr = nullptr;

    // Number of elements of type T in device memory.
    size_t count = 0;

public:
    /**
     * @brief Allocates device memory for n elements of type T.
     * @param n Number of elements to allocate.
     */
    explicit DeviceBuffer(size_t n) {
        count = n;
        CHECK_CUDA(
            cudaMalloc(&ptr, n * sizeof(T))
        );
    }

    /**
     * @brief Frees the device memory.
     * Buffers may end up containing a nullptr which is safe.
     */
    ~DeviceBuffer() { cudaFree(ptr); }

    // Forbid copying via copy constructor
    DeviceBuffer(const DeviceBuffer&) = delete;

    // Forbid copying via copy assignment
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    /**
     * @brief Move constructor. Receives ownership of device memory from o.
     * @param o The buffer to transfer ownership from. o will be left with a safe nullptr.
     */
    DeviceBuffer(DeviceBuffer&& o) noexcept {
        ptr   = std::exchange(o.ptr, nullptr);
        count = std::exchange(o.count, 0);
    }

    /**
     * @brief Move assignment. Frees current memory, then receives ownership from o.
     * @param o The buffer to transfer ownership from. o will be left with a safe nullptr.
     * @return A reference to this buffer.
     */
    DeviceBuffer& operator=(DeviceBuffer&& o) noexcept {
        if (this != &o) {
            cudaFree(ptr);
            ptr   = std::exchange(o.ptr, nullptr);
            count = std::exchange(o.count, 0);
        }
        return *this;
    }

    /** @brief Returns a pointer to the device memory. */
    T* data() { return ptr; }
    /** @brief Returns a const pointer to the device memory. */
    const T* data() const { return ptr; }

    /** @brief Returns the number of elements in this buffer. */
    size_t size() const { return count; }
    /** @brief Returns the total size of the device memory in bytes. */
    size_t size_bytes() const { return count * sizeof(T); }
};

/**
 * @brief RAII wrapper for a pinned host memory allocation of N elements of type T.
 *
 * Allocates pinned host memory on construction and automatically frees it on destruction.
 * Ownership of the memory is unique, copying is forbidden and moving transfers ownership.
 *
 * @tparam T The element type of the data in the buffer.
 */
template <typename T>
class HostBuffer {
    // Raw pointer to the start of the host memory.
    T* ptr = nullptr;

    // Number of elements of type T in host memory.
    size_t count = 0;

public:
    /**
     * @brief Allocates pinned host memory for n elements of type T.
     * @param n Number of elements to allocate.
     */
    explicit HostBuffer(size_t n) {
        count = n;
        CHECK_CUDA(
            cudaMallocHost(&ptr, n * sizeof(T))
        );
    }

    /**
     * @brief Frees the host memory.
     * Buffers may end up containing a nullptr which is safe.
     */
    ~HostBuffer() { cudaFreeHost(ptr); }

    // Forbid copying via copy constructor
    HostBuffer(const HostBuffer&) = delete;

    // Forbid copying via copy assignment
    HostBuffer& operator=(const HostBuffer&) = delete;

    /**
     * @brief Move constructor. Receives ownership of host memory from o.
     * @param o The buffer to transfer ownership from. o will be left with a safe nullptr.
     */
    HostBuffer(HostBuffer&& o) noexcept {
        ptr   = std::exchange(o.ptr, nullptr);
        count = std::exchange(o.count, 0);
    }

    /**
     * @brief Move assignment. Frees current memory, then receives ownership from o.
     * @param o The buffer to transfer ownership from. o will be left with a safe nullptr.
     * @return A reference to this buffer.
     */
    HostBuffer& operator=(HostBuffer&& o) noexcept {
        if (this != &o) {
            cudaFreeHost(ptr);
            ptr   = std::exchange(o.ptr, nullptr);
            count = std::exchange(o.count, 0);
        }
        return *this;
    }

    /** @brief Returns a pointer to the host memory. */
    T* data() { return ptr; }
    /** @brief Returns a const pointer to the host memory. */
    const T* data() const { return ptr; }

    /** @brief Returns the number of elements in this buffer. */
    size_t size() const { return count; }
    /** @brief Returns the total size of the host memory in bytes. */
    size_t size_bytes() const { return count * sizeof(T); }
};
