#pragma once

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
      std::printf("A CUDA API failed:\n      \'%s\' at:\n      %s:%d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}