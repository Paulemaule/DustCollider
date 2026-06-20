#pragma once

///////////////////////// LOGGING CONFIG /////////////////////////

// How many progress reports will be printed in total.
#define PROGRESS_LOG_AMOUNT 5
// The number of iterations that are skipped before the first progress report.
#define PROGRESS_LOG_OFFSET 5
// The weight of the rolling average algorithm used to determine the time per iteration.
#define ROLLING_AVERAGE_WEIGHT 0.01

///////////////////////// CUDA CONFIG /////////////////////////

// A macro for the number of threads per block for CUDA-Kernel execution.
#define BLOCK_SIZE 128

///////////////////////// SYSTEM CONFIG /////////////////////////

#ifdef _WIN32
    /**
     * A macro for the path separator \\.
     */
    #define SEP '\\'
#elif __linux__
    /** 
     * A macro for the path separator /.
     */
    #define SEP '/'
#endif

///////////////////////// ALGORITHM CONFIG /////////////////////////

// Macro for a very positive value, used in minimum number algorithm.
#define MIN_DEF  1e100

// Macro for a very negative value, used in maximum number algorithm.
#define MAX_DEF -1e100