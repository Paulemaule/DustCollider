#pragma once

///////////////////////// ERRORS /////////////////////////

// A makro that throws a runtime error with a predefined structure and variable description.
#define PANIC(description) throw std::runtime_error(std::string("ERROR: ") + description + " | " + "Source: " + __FILE__ + ":" + std::to_string(__LINE__) + "\n")

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
    