#pragma once

///////////////////////// LOGGING CONFIG /////////////////////////

#define VERBOSITY 3

///////////////////////// CUDA CONFIG /////////////////////////

// A makro that enables code execution on the device.
#define RUN_ON_GPU

// A makro for the number of threads per block for CUDA-Kernel execution.
#define BLOCK_SIZE 256

///////////////////////// SYSTEM CONFIG /////////////////////////

#ifdef _WIN32
    /**
     * A makro for the path seperator \\.
     */
    #define SEP '\\'
#elif __linux__
    /** 
     * A makro for the path seperator /.
     */
    #define SEP '/'
#endif

///////////////////////// ALGORITHM CONFIG /////////////////////////

// Makro for a very negative value, used in minimum number algorithm.
#define MIN_DEF  1e100

// Makro for a very negative value, used in maximum number algorithm.
#define MAX_DEF -1e100