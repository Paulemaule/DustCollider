#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>

#include "../utils/logging.cuh"
#include "../utils/errors.cuh"
#include "simulationConfig.cuh"

class CommandFile{
public:
    const std::string command_file_path;

    /**
     * @brief Constructor of the CommandFile class.
     * 
     * @param path: The path of the command file.
     */
    CommandFile(const std::string& path) : command_file_path(path) {}

    ~CommandFile() {}

    /**
     * @brief Parses the command file and fills a SimulationConfig struct with the run parameters it finds.
     * 
     * This function will read the command file specified when creating this CommandFile object and parse 
     * the run parameters. The run parameters are then stored in a SimulationConfig object.
     * The command file is read line by line and each line is sanitized before attempting to match a tag to the line.
     * Anything after a hashtag in a line is considered a comment and will not be parsed.
     * 
     * @param *out_config: Reference to a SimulationConfig object that will be filled with the parameters
     * read from the command file.
     */
    Status parse(SimulationConfig& out_config) {
        // Open the command file as an ifstream object
        std::ifstream command_file(command_file_path);

        // Check ifstream health
        if ( !command_file.is_open() ) {
            Logger::error("Could not open command file '{}'.", command_file_path);
            return Status::error;
        }

        // Parse the command file line by line
        std::string line;
        int line_counter = -1;

        while ( std::getline(command_file, line) ) {
            // Increment line counter
            line_counter++;

            // Sanitize the line
            line = sanitize_line(line);
            
            // If the line is empty after sanitation skip it
            if ( line.empty() ) continue;

            // Extract values from the line and place them in the SimulationConfig object
            Status _s = extract_values(line, out_config);
            if ( _s != Status::ok ) {
                return _s;          // Propagate any error that occurred during parameter extraction
            }
        }

        // Return Status:ok by default.
        return Status::ok;
    }

private:
    /**
     * @brief Extracts the first substring enclosed by "" from a line of the command file.
     * 
     * @param line: The line of the command file from which a substring is to be extracted.
     * 
     * @returns The substring found in the line. If no substring is found, return an empty string.
     */
    std::string extract_substring(std::string& line) {
        std::string::size_type pos1 = 0, pos2 = 0;
        std::string substring = "";
        int len = -1;

        // Search for the first "..." pattern and extract substring.
        if (line.find_first_of("\"") != std::string::npos)
        {
            pos1 = line.find("\"");
            pos2 = line.find("\"", pos1 + 1);

            // If the line contains only one '"' the function will simply return an empty string as the substring
            if ( pos2 == std::string::npos ) 
                return substring;

            len = int(pos2 - pos1 + 1);

            if (len < 0) 
                return substring;

            substring = line.substr(pos1, len);
            line.erase(pos1, len);
        }

        // Remove '"' from the extracted substring.
        while (substring.find("\"") != std::string::npos)
        {
            pos1 = substring.find("\"");
            substring.erase(pos1, 1);
        }

        return substring;
    }

    /**
     * @brief Sanitize and normalize a single line from the command file.
     *
     * This function transforms a raw input line into a normalized form suitable
     * for tag-based parsing. The following steps are taken:
     *
     * - Temporarily extracts quoted substrings of the form '"..."' to protect their
     *      contents from modification.
     * - Inserts whitespace around structural characters.
     * - Removes or replaces unwanted characters.
     * - Remove excessive whitespaces, by trimming consecutive whitespaces as well as 
     *      leading and trailing ones.
     * - Replaces commas with dots to enforce a consistent decimal separator.
     * - Removes inline comments starting with '#' or '!'. If a comment marker appears
     *      at the beginning of the line, any previously extracted quoted substring is
     *      discarded.
     * - Re-attaches a previously extracted quoted substring (if valid) at the end of
     *   the line.
     *
     * Empty or comment-only lines are reduced to an empty string.
     *
     * @param[in,out] line
     *     The input line to be sanitized. The contents are modified in place during
     *     processing.
     *
     * @return std::string
     *     The sanitized and normalized version of the input line. This is typically
     *     identical to `line` after modification, but is returned for convenience and
     *     chaining.
     *
     * @note
     *     Encapsulated substrings might cause issues if misused.
     *     They are meant only for encapsulating paths.
     */
    std::string sanitize_line(std::string& line) {
        std::string::size_type pos = 0;

        // Search for and remove substrings of the form "..."
        // They are added back later
        std::string extracted_substring = extract_substring(line);

        if (line.size() == 0)
            return line;

        // Replace '>' with '> '
        if (line.find(">") != std::string::npos)
        {
            pos = line.find(">");
            line.replace(pos, 1, "> ");
        }

        // Replace '=' with ' = '
        if (line.find("=") != std::string::npos)
        {
            pos = line.find("=");
            line.replace(pos, 1, " = ");
        }

        // Delete ';'
        while (line.find(";") != std::string::npos)
        {
            pos = line.find(";");
            line.replace(pos, 1, " ");
        }

        // Delete '?'
        while (line.find("?") != std::string::npos)
        {
            pos = line.find("?");
            line.replace(pos, 1, " ");
        }

        // Delete '*'
        while (line.find("*") != std::string::npos)
        {
            pos = line.find("*");
            line.replace(pos, 1, " ");
        }

        // Delete '\t'
        while (line.find('\t') != std::string::npos)
        {
            pos = line.find('\t');
            line.replace(pos, 1, " ");
        }

        // Delete ' \r\n'
        while (line.find(" \r\n") != std::string::npos)
        {
            pos = line.find(" \r\n");
            line.replace(pos, 3, " ");
        }

        // Delete ' \r'
        while (line.find(" \r") != std::string::npos)
        {
            pos = line.find(" \r");
            line.replace(pos, 2, " ");
        }

        // Delete ' \n'
        while (line.find(" \n") != std::string::npos)
        {
            pos = line.find(" \n");
            line.replace(pos, 2, " ");
        }

        // Delete '/r/n'
        while (line.find("\r\n") != std::string::npos)
        {
            pos = line.find("\r\n");
            line.replace(pos, 2, " ");
        }

        // Delete '/r'
        while (line.find("\r") != std::string::npos)
        {
            pos = line.find("\r");
            line.replace(pos, 1, " ");
        }

        // Delete '\n'
        while (line.find("\n") != std::string::npos)
        {
            pos = line.find("\n");
            line.replace(pos, 1, " ");
        }

        // Replace all '  ' with ' ' (double space with space)
        while (line.find("  ") != std::string::npos)
        {
            pos = line.find("  ");
            line.replace(pos, 2, " ");
        }

        // Replace all ',' with '.'
        while (line.find(",") != std::string::npos)
        {
            pos = line.find(",");
            line.replace(pos, 1, ".");
        }

        if (line == " ")
            line = "";

        // Remove trailing ' '
        if (line.size() > 0)
        {
            while (line.c_str()[line.size() - 1] == ' ')
            {
                pos = line.find_last_of(' ');
                line.erase(pos, 1);
            }
        }

        // Remove leading ' '
        if (line.size() > 0) {
            while ( line.c_str()[0] == ' ' )
            {
                line.erase(0, 1);
            }
        }

        // If the line contains '#', everything after it will be removed.
        // If the first character is '#', extracted substrings ("...") will be removed.
        /**
         * TODO: The handling of substrings is not robust enough.
         */
        if (line.find_first_of("#") != std::string::npos)
        {
            pos = line.find("#");

            if (pos == 0) {
                extracted_substring.clear(); 
            }

            line.erase(pos, line.length() - pos);
        }

        // Same as with '#' comments.
        if (line.find_first_of("!") != std::string::npos)
        {
            pos = line.find("!");

            // If the '#' was found at the start of the line any extracted substring is considered invalid
            if (pos == 0)
                extracted_substring = "";

            line.erase(pos, line.length() - pos);
        }

        // If there was a substring of the type "..." add it back to the line.
        if ( extracted_substring.size() > 0 )
            line += " \"" + extracted_substring + "\"";


        // Return the sanitized line.
        return line;
    }

    /**
     * @brief Extract a command tag and its associated value from a sanitized command line.
     *
     * This function interprets a single, already-sanitized line from the command file.
     * The lines are expected to have the form:
     *     <tag> value
     *
     * The extracted tag and value are matched against known configuration
     * options and written into the provided SimulationConfig struct.
     *
     * If the expected tag format cannot be found, the function reports an error
     * and aborts parsing.
     *
     * @param[in,out] line
     *     A sanitized command line. The contents are not modified.
     *
     * @param[out] out_config
     *     Reference to the SimulationConfig struct into which the run parameters are written.
     *
     * @return Status::ok
     *     If a valid tag–value pair was found and successfully processed.
     *
     * @return Status::error
     *     If the line does not contain a valid command tag in the form `<...>`, or
     *     if command matching fails.
     */
    Status extract_values(std::string& line, SimulationConfig& out_config) {
        // Positions of the tag start and end
        std::string::size_type l_pos;
        std::string::size_type r_pos;

        std::string command_tag;
        std::string command_value;

        // Search for a tag in the line.
        l_pos = line.find("<");
        r_pos = line.find(">", l_pos);

        if ( l_pos != std::string::npos && r_pos != std::string::npos && r_pos > l_pos) {
            // Extract the command tag.
            command_tag = line.substr(l_pos + 1, r_pos - l_pos - 1);
            // Extract the command value. The command value starts at the first non whitespace character after the closing bracket '>'.
            command_value = line.substr(r_pos + 1);
            std::string::size_type value_start = command_value.find_first_not_of(' ');
            command_value = (value_start == std::string::npos) ? "" : command_value.substr(value_start);
        } else {
            Logger::error("Cannot find a tag in command file line: '{}'", line.c_str());
            return Status::error;
        }

        // Match the command tag and value into
        Status _s = process_tag_value(command_tag, command_value, out_config);
        if ( _s != Status::ok ) {
            return _s;          // Propagate errors.
        }

        // Return Status::ok by default.
        return Status::ok;
    }

    /**
     * @brief Apply a parsed '<tag>' and its associated value to a SimulationConfig.
     * 
     * @param[in] tag
     *     The command tag name (without angle brackets), e.g. '"pos_A"'.
     *     The string is not modified.
     *
     * @param[in] value
     *     The raw value string associated with the tag. The string is not modified.
     *
     * @param[out] out_config
     *     Configuration struct to be updated according to the parsed tag/value.
     *
     * @return Status::ok
     *     If the tag is recognized and the value was successfully parsed and stored.
     *
     * @return Status::error
     *     If the tag is unknown, the value has an invalid format (wrong number of
     *     components, empty path, failed conversion), or if a called function reports an error.
     */
    Status process_tag_value(std::string& tag, std::string& value, SimulationConfig& out_config) {
        // The output path of the simulation
        if ( tag == "path_results" ) {
            std::string path = extract_substring(value);

            // Check for errors.
            if ( path.size() == 0 ) {
                Logger::error("No path found in <path_results>.");
                return Status::error;
            }

            // Store results and return Ok status.
            out_config.output.path = path;
            return Status::ok;
        }

        // Aggregate parameter parsing
        if ( tag.find("aggregate_") != std::string::npos ) {
            // Extract the aggregate key: "aggregate_A_path" -> "A"
            size_t first  = tag.find('_');
            size_t second = tag.find('_', first + 1);
            std::string aggregate_key = tag.substr(first + 1, second - first - 1);

            // Find an existing entry or create a new one; work through a pointer so
            // modifications are written directly into the vector element.
            AggregateConfig* cfg = nullptr;
            for ( AggregateConfig& a : out_config.aggregates ) {
                if ( a.name == aggregate_key ) { cfg = &a; break; }
            }
            if ( cfg == nullptr ) {
                out_config.aggregates.push_back(AggregateConfig{});
                out_config.aggregates.back().name = aggregate_key;
                cfg = &out_config.aggregates.back();
            }

            if ( tag.find("_path") != std::string::npos ) {
                std::string path = extract_substring(value);
                if ( path.empty() ) {
                    Logger::error("No path found in <{}>.", tag);
                    return Status::error;
                }
                cfg->path = path;
                return Status::ok;
            }

            if ( tag.find("_pos") != std::string::npos ) {
                std::vector<double> components;
                Status _s = split_values(value, components);
                if ( _s != Status::ok ) return _s;
                if ( components.size() != 3 ) {
                    Logger::error("<{}> needs to be three dimensional.", tag);
                    return Status::error;
                }
                cfg->position = { components[0], components[1], components[2] };
                return Status::ok;
            }

            if ( tag.find("_vel") != std::string::npos ) {
                std::vector<double> components;
                Status _s = split_values(value, components);
                if ( _s != Status::ok ) return _s;
                if ( components.size() != 3 ) {
                    Logger::error("<{}> needs to be three dimensional.", tag);
                    return Status::error;
                }
                cfg->velocity = { components[0], components[1], components[2] };
                return Status::ok;
            }

            if ( tag.find("_ang") != std::string::npos ) {
                std::vector<double> components;
                Status _s = split_values(value, components);
                if ( _s != Status::ok ) return _s;
                if ( components.size() != 3 ) {
                    Logger::error("<{}> needs to be three dimensional.", tag);
                    return Status::error;
                }
                cfg->angular = { components[0], components[1], components[2] };
                return Status::ok;
            }
        }

        // Simulation parameters
        if (tag == "N_iter") {
            int N_iter;

            // Attempt to convert value to integer
            try {
                N_iter = std::stoi(value);
            } catch (...) {
                Logger::error("Could not convert <N_iter> '{}' to integer.", value);
                return Status::error;
            }

            // Store results and return Ok status.
            out_config.N_iter = N_iter;
            return Status::ok;
        }

        if (tag == "N_save") {
            int N_save;

            // Attempt to convert value to integer
            try {
                N_save = std::stoi(value);
            } catch (...) {
                Logger::error("Could not convert <N_save> '{}' to integer.", value);
                return Status::error;
            }

            // Store results and return Ok status.
            out_config.output.N_save = N_save;
            return Status::ok;
        }

        if (tag == "save_ovito") {
            bool save;

            Status _s = to_bool(value, save);
            if ( _s != Status::ok ) {
                return _s;
            }
            
            // Store results and return Ok status.
            out_config.output.ovito = save;
            return Status::ok;
        }

        if (tag == "save_position") {
            bool save;

            Status _s = to_bool(value, save);
            if ( _s != Status::ok ) {
                return _s;
            }
            
            // Store results and return Ok status.
            out_config.output.position = save;
            return Status::ok;
        }

        if (tag == "save_velocity") {
            bool save;

            Status _s = to_bool(value, save);
            if ( _s != Status::ok ) {
                return _s;
            }
            
            // Store results and return Ok status.
            out_config.output.velocity = save;
            return Status::ok;
        }

        if (tag == "save_angular") {
            bool save;

            Status _s = to_bool(value, save);
            if ( _s != Status::ok ) {
                return _s;
            }
            
            // Store results and return Ok status.
            out_config.output.angular = save;
            return Status::ok;
        }

        if (tag == "save_force") {
            bool save;

            Status _s = to_bool(value, save);
            if ( _s != Status::ok ) {
                return _s;
            }
            
            // Store results and return Ok status.
            out_config.output.force = save;
            return Status::ok;
        }

        if (tag == "save_torque") {
            bool save;

            Status _s = to_bool(value, save);
            if ( _s != Status::ok ) {
                return _s;
            }
            
            // Store results and return Ok status.
            out_config.output.torque = save;
            return Status::ok;
        }

        if (tag == "save_energy") {
            bool save;

            Status _s = to_bool(value, save);
            if ( _s != Status::ok ) {
                return _s;
            }
            
            // Store results and return Ok status.
            out_config.output.energy = save;
            return Status::ok;
        }

        // Material parameter parsing.
        // Command file format: <material id="N">  "name"  gamma E nu rho xi tvis [tss tsl Msat chi Tc]
        // After sanitize_line the tag becomes "material id = " and the value has the name first,
        // numeric params in the middle, and the id string at the end (moved there by the sanitizer).
        if ( tag.rfind("material", 0) == 0 ) {
            // Extract material name (first quoted token in value)
            std::string name = extract_substring(value);
            if ( name.empty() ) {
                Logger::error("No material name found in <{}>.", tag);
                return Status::error;
            }

            // Extract material ID (second quoted token; the sanitizer moved it to the end of value)
            std::string id_str = extract_substring(value);
            if ( id_str.empty() ) {
                Logger::error("No material ID found in <{}>.", tag);
                return Status::error;
            }
            
            int mat_id;
            try {
                mat_id = std::stoi(id_str);
            } catch (...) {
                Logger::error("Could not convert material ID '{}' to integer.", id_str);
                return Status::error;
            }

            // Parse the numeric material parameters
            std::vector<double> params;
            Status _s = split_values(value, params);
            if ( _s != Status::ok ) return _s;

            MaterialEntry mat;
            mat.name = name;
            if ( params.size() == 6 ) {
                // Non-magnetic: gamma [J/m²], E [Pa], nu [-], rho [kg/m³], xi [m], tvis [s]
                mat.gamma = params[0];
                mat.E     = params[1];
                mat.nu    = params[2];
                mat.rho   = params[3];
                mat.xi    = params[4];
                mat.tvis  = params[5];
            } else if ( params.size() == 11 ) {
                // Magnetic: as above plus tss [s], tsl [s], Msat [A/m], chi [-], Tc [K]
                mat.gamma = params[0];
                mat.E     = params[1];
                mat.nu    = params[2];
                mat.rho   = params[3];
                mat.xi    = params[4];
                mat.tvis  = params[5];
                mat.tss   = params[6];
                mat.tsl   = params[7];
                mat.Msat  = params[8];
                mat.chi   = params[9];
                mat.Tc    = params[10];
            } else {
                Logger::error("Material '{}' has {} parameters, expected 6 (non-magnetic) or 11 (magnetic).",
                    name, params.size());
                return Status::error;
            }

            // IDs are 1-indexed in the command file; store at (id - 1)
            int idx = mat_id - 1;
            if ( idx < 0 ) {
                Logger::error("Material ID must be >= 1, got {}.", mat_id);
                return Status::error;
            }
            if ( idx >= static_cast<int>(out_config.materials.size()) ) {
                out_config.materials.resize(idx + 1);
            }
            out_config.materials[idx] = mat;

            return Status::ok;
        }

        if (tag == "B_ext") {
            std::vector<double> components;
            Status _s = split_values(value, components);
            if ( _s != Status::ok ) return _s;
            if ( components.size() != 3 ) {
                Logger::error("<B_ext> needs to be three dimensional.");
                return Status::error;
            }
            out_config.B_ext = { components[0], components[1], components[2] };
            return Status::ok;
        }

        if (tag == "T_dust") {
            try {
                out_config.T_dust = std::stod(value);
            } catch (...) {
                Logger::error("Could not convert <T_dust> to double.");
                return Status::error;
            }
            return Status::ok;
        }

        if (tag == "timestep") {
            try {
                out_config.timestep = std::stod(value);
            } catch (...) {
                Logger::error("Could not convert <timestep> to double.");
                return Status::error;
            }
            return Status::ok;
        }

        // Tag matches no known pattern
        Logger::error("Unknown tag '{}' with value '{}'", tag.c_str(), value.c_str());
        return Status::error;
    }

    /**
     * @brief Splits a whitespace-separated list of numeric values and convert them to doubles.
     *
     * @param[in] values
     *     Input string containing one or more numeric values separated by any number of
     *     whitespaces.
     *
     * @param[out] out
     *     Vector that will be filled with the parsed numeric values. Any previous
     *     contents are discarded.
     *
     * @return Status::ok
     *     If all tokens were successfully converted to doubles.
     *
     * @return Status::error
     *     If any token cannot be converted to a double.
     *
     * @note
     *     This function does not care about the dimensionality of the final vector.
     *     Care needs to be taken to check the dimension from the caller side.
     */
    Status split_values(std::string& values, std::vector<double>& out) {
        std::istringstream iss(values);
        out.clear();

        for ( std::string t; iss >> t; ) {
            try {
                out.push_back(std::stod(t));
            } catch (...) {
                Logger::error("Error when attemting to convert the string '{}' to double.", t);
                return Status::error;
            }
        }

        return Status::ok;
    }

    /**
     * @brief Convert a string to a boolean value.
     *
     * The conversion from string to boolean only accepts the follwoing values:
     *      {'true', 'True', '1'} and {'false', 'False', '0'}
     *
     * Any other input is considered invalid and results in an error.
     *
     * @param[in] s
     *     Input string to be interpreted as a boolean value.
     *
     * @param[out] out
     *     Boolean variable that will receive the converted value on success.
     *     Its value is undefined if the function returns `Status::error`.
     *
     * @return Status::ok
     *     If the input string matches a supported boolean literal.
     *
     * @return Status::error
     *     If the input string cannot be interpreted as a boolean.
     */
    Status to_bool(const std::string& s, bool& out) {
        if (s == "1" || s == "true"  || s == "True") {
            out = true;
            return Status::ok;
        }
        else if (s == "0" || s == "false" || s == "False") {
            out = false;
            return Status::ok;
        }
        else {
            out = false;
            Logger::error("The string '{}' could not be converted to a boolean value. 'False' was set by default.", s);
            return Status::error;
        }
    }
};