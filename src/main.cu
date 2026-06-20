/**
 * @file main.cu
 * @brief Entry point for the dustCollider code.
 * 
 * Orchestrates the three phases of the code:
 *  1) The Setup        - Simulation setup from input files
 *  2) The Simulation   - Run the time integration loop
 *  3) The Output       - Write stored snapshots to disk
 * 
 * Error handling is inconsistent at this point.
 * Sometimes errors are thrown, sometimes the code returns -1.
 */

#include <chrono>

#include "utils/logging.cuh"
#include "utils/errors.cuh"

#include "simulationSetup/simulationConfig.cuh"
#include "simulationSetup/simulationSetup.cuh"
#include "simulator.cuh"

int main(const int argc, const char** argv)
{
    auto start = std::chrono::high_resolution_clock::now();

    // Run the simulation setup
    Pipeline setup;

    Status _s = setup.run(argc, argv);
    if (_s != Status::ok) {
        Logger::error("An error occurred during simulation setup. Terminating.");
        return -1;
    }

    // Initialize the simulation
    Simulator simulator(setup.getConfig());

    // Start the simulation loop
    simulator.run();

    // Write the output files
    simulator.write_output();

    // Calculate and print the final runtime.
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    char buffer[14];
    ns_to_time_string(elapsed.count(), buffer, 14);
    Logger::print("Total runtime : {} .", buffer);
    
    Logger::lineBreak();
    Logger::header("DONE");

    return 0;
}