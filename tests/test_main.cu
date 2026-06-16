/**
 * @file test_main.cu
 * @brief Entry point for the unit test suite.
 *
 * All test files are included directly into this single translation unit to
 * avoid multiple-definition linker errors that would arise from the
 * non-inline functions defined in .cuh headers.
 *
 * Adding a new test module
 * ------------------------
 * 1. Create tests/test_definitions/unit/test_<module>.cu with a void test_<module>() function.
 * 2. Add  #include "test_definitions/unit/test_<module>.cu"  below.
 * 3. Add  RUN_SUITE("<module>", test_<module>);  in main().
 */

#include "framework/testkit.h"

// --- test module includes go here ---
#include "test_definitions/unit/test_vectors.cu"
#include "test_definitions/unit/test_math.cu"
#include "test_definitions/unit/test_cluster.cu"
#include "test_definitions/unit/test_parser.cu"
#include "test_definitions/unit/test_aggregate.cu"
#include "test_definitions/unit/test_pipeline.cu"
#include "test_definitions/unit/test_buffer.cu"

int main() {
    // --- RUN_SUITE calls go here ---
    RUN_SUITE("vector utilities",     test_vectors);
    RUN_SUITE("math utilities",       test_math);
    RUN_SUITE("cluster detection",    test_cluster);
    RUN_SUITE("command file parser",  test_parser);
    RUN_SUITE("aggregate file I/O",   test_aggregate);
    RUN_SUITE("pipeline setup",       test_pipeline);
    RUN_SUITE("GPU memory wrappers",  test_buffer);

    return testkit::summary();
}
