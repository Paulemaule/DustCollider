#pragma once

///////////////////////// BUILD CONFIG /////////////////////////

///////////////////////// LOGGING CONFIG /////////////////////////

#define VERBOSITY 3

// How many progress reports will be printed in total.
#define PROGRESS_LOG_AMMOUNT 5
// The number of iterations that are skipped before the first progress report.
#define PROGRESS_LOG_OFFSET 5
// The weight of the rolling average algorithm used to determine the time per iteration.
#define ROLLING_AVERAGE_WEIGHT 0.01

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