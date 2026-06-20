#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "../utils/logging.cuh"
#include "../utils/errors.cuh"

/**
 * @brief A struct that contains header information about an aggregate snapshot.
 */
struct AggregateHeader {
    int     Nmon;                           // The number of monomer in the aggregate [m].
    double  external_radius;                // The external radius of the aggregate [m].
    double  effective_radius;               // The effective radius of the aggregate [m].
};

/**
 * @brief A struct containing all relevant information for a snapshot of a single aggregate.
 */
struct AggregateMonomers {
    std::vector<double3> positions;         // The positions of the monomers [m].
    std::vector<double> radii;              // The radii of the monomers [m].
    std::vector<int> material_ids;          // The material IDs of the monomers.
};

/**
 * @brief This class handles reading and writing aggregate snapshots to and from files.
 */
class Aggregate {
public:
    // Class destructor.
    ~Aggregate() {}; 

    // These fields hold the aggregates information
    AggregateHeader header;                 // Meta information about the aggregate.
    AggregateMonomers monomers;             // Information about the individual aggregate monomers.

    /**
     * @brief Construct an Aggregate object from a file.
     * 
     * This function will read an aggregate file at the specified path and construct 
     * an Aggregate object.
     */
    static Aggregate from_file(const std::string& aggregate_file_path) {
        // The instance of Aggregate that will be filled and returned by this function.
        Aggregate out = Aggregate();

        std::ifstream aggregate_file(aggregate_file_path);

        // Check ifstream health
        if ( !aggregate_file.is_open() ) {
            Logger::error("Could not open aggregate file {}", aggregate_file_path);
            throw std::runtime_error("Failed to open aggregate file.");
        }

        std::string line;
        int line_counter = -1;                  // Start at -1 due to line_counter++ at the start of each loop iteration.

        // Iterate over the lines of the aggregate file.
        while ( std::getline(aggregate_file, line) ) {
            line_counter++;
            std::vector<double> line_contents;

            // Line 0: header — Nmon, external_radius [nm], effective_radius [nm]
            if ( line_counter == 0 ) {
                if ( parse_line_values(line, line_contents) != Status::ok ) {
                    Logger::error("Could not read aggregate header line.");
                    throw std::runtime_error("Could not read aggregate header line.");
                }

                if ( line_contents.size() < 3 ) {
                    Logger::error("Aggregate header needs at least 3 values, not {}", line_contents.size());
                    throw std::runtime_error("Aggregate header too short.");
                }

                out.header.Nmon             = static_cast<int>(line_contents[0]);
                out.header.external_radius  = line_contents[1] * 1e-9;
                out.header.effective_radius = line_contents[2] * 1e-9;
            }

            // Lines 1-5: metadata / comments — skip.

            // Lines 5+: monomer data — x y z ? radius ? mat_id  (all in nm, mat_id 1-indexed)
            if ( line_counter > 4 ) {
                if ( parse_line_values(line, line_contents) != Status::ok ){
                    Logger::error("Could not read aggregate monomer line. '{}'", line);
                    throw std::runtime_error("Could not read aggregate monomer line.");
                }

                if ( line_contents.size() != 7 ) {
                    Logger::error("Aggregate monomer line '{}' has {} values, expected 7.",
                        line_counter, line_contents.size());
                    throw std::runtime_error("Aggregate monomer line has wrong number of values.");
                }

                double3 position;
                position.x = 1e-9 * line_contents[0];
                position.y = 1e-9 * line_contents[1];
                position.z = 1e-9 * line_contents[2];

                double radius      = 1e-9 * line_contents[4];
                int    material_id = static_cast<int>(line_contents[6]) - 1;   // convert to 0-indexed

                out.monomers.positions.push_back(position);
                out.monomers.radii.push_back(radius);
                out.monomers.material_ids.push_back(material_id);
            }
        }

        return out;
    }

    /**
     * @brief A function that would construct an instances of Aggregate from the system state.
     * TODO: Implement
     */
    static Aggregate from_state() {
        Aggregate out = Aggregate();

        throw std::runtime_error("Not implemented yet.");

        return out;
    }

    /**
     * @brief A function that would store an instance of Aggregate to the disk.
     * TODO: Implement
     */
    Status to_file(std::string& aggregate_file_path) {
        throw std::runtime_error("Not implemented yet.");

        return Status::error;
    }
    
private:
    // Constructors for this class should not be publicly accessible.
    Aggregate() {};

    /**
     * @brief Splits a whitespace-separated line into its double values.
     *
     * @param[in] line    Input string with numeric tokens.
     * @param[out] out    Vector filled with the parsed values; cleared before use.
     *
     * @return Status::ok    All tokens converted successfully.
     * @return Status::error A token could not be converted to double.
     */
    static Status parse_line_values(const std::string& line, std::vector<double>& out) {
        std::istringstream iss(line);
        out.clear();

        for ( std::string t; iss >> t; ) {
            try {
                out.push_back(std::stod(t));
            } catch (...) {
                Logger::error("Could not convert '{}' to double in aggregate file.", t);
                return Status::error;
            }
        }

        return Status::ok;
    }
};