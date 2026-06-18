#pragma once

#include <cstdio>

#include "logging.cuh"

/////////////////////////

/**
 * @brief This enum is used by the Pipeline to track errors.
 * 
 * A function should return Status::ok if the function executed without problems.
 * If an error was encountered during the execution the function should instead return Status::error.
 */
enum class Status {
    ok,
    error
};

///////////////////////// ERRORS /////////////////////////

// A makro that throws a runtime error with a predefined structure and variable description.
#define PANIC(description) throw std::runtime_error(std::string("ERROR: ") + description + " | " + "Source: " + __FILE__ + ":" + std::to_string(__LINE__) + "\n")

///////////////////////// CUDA API ERRORS /////////////////////////

// A makro that checks for the last cuda error and throws an exception if one is encountered.
#define CUDA_LAST_ERROR_CHECK()                                              \
    do {                                                                     \
        cudaError_t err = cudaGetLastError();                                \
        if (err != cudaError::cudaSuccess) {                                 \
            throw std::runtime_error(std::string("CUDA ERROR: ") +           \
                                     cudaGetErrorString(err) + " | "         \
                                     "Detection: " + __FILE__ + ":" +        \
                                     std::to_string(__LINE__));              \
        }                                                                    \
    } while (0)

// A wrapper makro that will check the return code of a CUDA API call.
#define CHECK_CUDA(ans) { checkErrorCode((ans), __FILE__, __LINE__); }

// A function that will check the return code of a CUDA API call.
inline void checkErrorCode(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaError::cudaSuccess) {
        Logger::error("A CUDA API call failed at {}:{} with Errorcode {} ({})", file, line, (int)code, cudaGetErrorString(code));
        if (abort) exit(code);
   }
}