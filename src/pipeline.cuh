#pragma once

#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#endif

#include "utils/printing.cuh"
#include "utils/vector.cuh"

/**
 * @brief A class for managing the simulation pipeline.
 * 
 * This class contains fields that represent the run parameters and functions to prepare, process, summarize and save the data.
 */
class CPipeline
{
public:
    CPipeline()
    {
        // Initializing the run parameters.
        pos_A.x = 0;
        pos_A.y = 0;
        pos_A.z = 0;

        pos_B.x = 0;
        pos_B.y = 0;
        pos_B.z = 0;

        vel_A.x = 0;
        vel_A.y = 0;
        vel_A.z = 0;

        vel_B.x = 0;
        vel_B.y = 0;
        vel_B.z = 0;

        ang_A.x = 0;
        ang_A.y = 0;
        ang_A.z = 0;

        ang_B.x = 0;
        ang_B.y = 0;
        ang_B.z = 0;

        N_iter = 0;
        N_save = 0;

        time_start = 0;
        time_stop = 0;
        time_step = 0;

        a_eff_A = 0;
        a_eff_B = 0;

        Nmon_A = 0;
        Nmon_B = 0;

        a_out_A = 0;
        a_out_B = 0;

        a_mon_min = 0;
        a_mon_max = 0;

        save_pos = false;
        save_vel = false;
        save_force = false;
        save_torque = false;
        save_cluster = false;
        save_omega = false;
        save_mag = false;
        save_ovito = false;

        Bext.x = 0;
        Bext.y = 0;
        Bext.z = 0;

        T_dust = 15;

        // Assume initially that the materials are non magnetic.
        mat_type = MAT_TYPE_NONE;

        // Initializing values for min/max parameter determination.
        min_gamma = MIN_DEF;
        min_E = MIN_DEF;
        min_nu = MIN_DEF;
        min_rho = MIN_DEF;
        min_xi = MIN_DEF;
        min_Tvis = MIN_DEF;

        min_tss = MIN_DEF;
        min_tsl = MIN_DEF;
        min_Msat = MIN_DEF;
        min_chi = MIN_DEF;
        min_Tc = MIN_DEF;

        max_gamma = MAX_DEF;
        max_E = MAX_DEF;
        max_nu = MAX_DEF;
        max_rho = MAX_DEF;
        max_xi = MAX_DEF;
        max_Tvis = MAX_DEF;

        max_tss = MAX_DEF;
        max_tsl = MAX_DEF;
        max_Msat = MAX_DEF;
        max_chi = MAX_DEF;
        max_Tc = MAX_DEF;
    };

    ~CPipeline() {}

    /**
     * @brief This function parses the command line input to extract the command file location.
     * 
     * @param argc: The count of command line arguments.
     * @param argv: The command line arguments.
     * @return True if the function was successful, False otherwise.
     */
    bool init(int argc, const char** argv)
    {
        PRINT_HEADLINE();
        PRINT_CLR_LINE();

        #ifdef RELEASE
        std::cout << "Compiled in Release build." << std::endl;
        PRINT_LOG("Parsing command line input", 2);
        if (argc != 2)
        {
            PRINT_ERROR("Wrong number of command line inputs. Only the command file location is required.");
            return false;
        }
        cmd_filename = argv[1];
        #else
        std::cout << "Compiled in Debug build." << std::endl;
        cmd_filename = DEBUG_CMD_FILE;
        PRINT_LOG("Assuming cmd file location is:\n      " + cmd_filename, 2)
        #endif

        return true;
    }

    /**
     * @brief A function that sanitizes a string for command file reading.
     * 
     * @param &line: A reference to the string that is to be formatted.
     */
    void formatLine(string& line) {
        string::size_type pos = 0;

        // Search for and remove substrings of the form "...".
        string tmp_str = seperateString(line);

        if (line.size() == 0)
            return;

        // Replace '>' with '> '
        if (line.find(">") != string::npos)
        {
            pos = line.find(">");
            line.replace(pos, 1, "> ");
        }

        // Replace '=' with ' = '
        if (line.find("=") != string::npos)
        {
            pos = line.find("=");
            line.replace(pos, 1, " = ");
        }

        // Delete ';'
        while (line.find(";") != string::npos)
        {
            pos = line.find(";");
            line.replace(pos, 1, " ");
        }

        // Delete '?'
        while (line.find("?") != string::npos)
        {
            pos = line.find("?");
            line.replace(pos, 1, " ");
        }

        // Delete '*'
        while (line.find("*") != string::npos)
        {
            pos = line.find("*");
            line.replace(pos, 1, " ");
        }

        // Delete '\t'
        while (line.find('\t') != string::npos)
        {
            pos = line.find('\t');
            line.replace(pos, 1, " ");
        }

        // Delete ' \r\n'
        while (line.find(" \r\n") != string::npos)
        {
            pos = line.find(" \r\n");
            line.replace(pos, 3, " ");
        }

        // Delete ' \r'
        while (line.find(" \r") != string::npos)
        {
            pos = line.find(" \r");
            line.replace(pos, 2, " ");
        }

        // Delete ' \n'
        while (line.find(" \n") != string::npos)
        {
            pos = line.find(" \n");
            line.replace(pos, 2, " ");
        }

        // Delete '/r/n'
        while (line.find("\r\n") != string::npos)
        {
            pos = line.find("\r\n");
            line.replace(pos, 2, " ");
        }

        // Delete '/r'
        while (line.find("\r") != string::npos)
        {
            pos = line.find("\r");
            line.replace(pos, 1, " ");
        }

        // Delete '\n'
        while (line.find("\n") != string::npos)
        {
            pos = line.find("\n");
            line.replace(pos, 1, " ");
        }

        // Replace '  ' with ' ' (double space with space)
        while (line.find("  ") != string::npos)
        {
            pos = line.find("  ");
            line.replace(pos, 2, " ");
        }

        // Replace ',' with '.'
        while (line.find(",") != string::npos)
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
        while (line.c_str()[0] == ' ')
        {
            line.erase(0, 1);
        }

        // If the line contains '#', everything after it will be removed.
        // If the first character is '#', extracted substrings ("...") will be removed.
        if (line.find_first_of("#") != string::npos)
        {
            pos = line.find("#");
            if (pos == 0)
                tmp_str = "";

            line.erase(pos, line.length() - pos);
        }

        // Same as with '#' comments.
        if (line.find_first_of("!") != string::npos)
        {
            pos = line.find("!");
            if (pos == 0)
                tmp_str = "";

            line.erase(pos, line.length() - pos);
        }

        // If there was a substring of the type "..." add it back to the line.
        if (tmp_str.size() != 0)
            line += " \"" + tmp_str + "\"";
    };

    /**
     * @brief Extracts the first substring of the type "..." from a string.
     * 
     * @param &str: A reference to the string that is to be extracted from. The substring will be removed from this.
     * @return The substring. If no such string is found return ''.
     */
    string seperateString(string& str)
    {
        string::size_type pos1 = 0, pos2 = 0;
        string ret = "";
        int len = -1;

        // Search for the first "..." pattern and extract substring.
        if (str.find_first_of("\"") != string::npos)
        {
            pos1 = str.find("\"");
            pos2 = str.find("\"", pos1 + 1);

            len = int(pos2 - pos1 + 1);

            if (len < 0)
                return ret;

            ret = str.substr(pos1, len);
            str.erase(pos1, len);
        }

        // Remove '"' from the extracted substring.
        while (ret.find("\"") != string::npos)
        {
            pos1 = ret.find("\"");
            ret.erase(pos1, 1);
        }

        return ret;
    };

    /**
     * @brief Split a string at each space and convert all substring into floats.
     * 
     * @param &str: The string from which numbers are to be extracted. Remains untouched.
     * @return A list of the floats. Strings that are not numbers are converted to 0.0
     */
    dlist parseValues(string& str)
    {
        int pos;
        dlist values;
        string v;

        formatLine(str);

        // If no pattern is found return an empty container.
        if (str.size() == 0)
            return values;

        // Split the string at ' ' and convert the substring into floats which are added to `values`.
        while (str.find(" ") != string::npos)
        {
            pos = int(str.find(" "));
            v = str.substr(0, pos);
            str.erase(0, pos + 1);
            values.push_back(atof(v.c_str()));
        }

        // When the string contains only one substring.
        values.push_back(atof(str.c_str()));

        return values;
    };

    /**
     * @brief Load the command file and parse the line one by one.
     * 
     * @return True on success, False on failure.
     */
    bool parse()
    {
        PRINT_LOG("Reading the command file.", 2);

        ifstream reader(cmd_filename.c_str());

        string line, cmd;
        string::size_type pos = 0;
        int line_counter = 0;

        if (reader.fail())
        {
            cout << "\nERROR: Cannot open command file:\n\t" << cmd_filename << endl;
            return false;
        }

        // Iterate over the lines in the command file.
        while (getline(reader, line))
        {
            line_counter++;

            formatLine(line);

            // When the line is empty skip
            if (line.compare("") == 0)
                continue;

            // When the line starts with '#' skip
            if (line.c_str()[0] == '#')
                continue;

            // Searching for the command file pattern '<tag> value'
            pos = line.find(">");

            if (pos != string::npos)
            {
                // Extract the tag and value from the pattern.
                cmd = line.substr(0, pos + 1);
                line.erase(0, pos + 2);
            }
            else
            {
                printf("ERROR: Cannot find a tag in command file line %d\n      \'%s\'", line_counter, line.c_str());
                return false;
            }

            // If the tag is emtpy, skip
            if (cmd.compare("") == 0)
                continue;

            // If the tag is a space, skip
            if (cmd.compare(" ") == 0)
                continue;

            // Parse the line and set run parameters accordingly.
            if (!parseLine(cmd, line))
            {
                printf("ERROR: Cannot parse the command \'%s\' in command file line %d\n      \'%s\'", cmd.c_str(), line_counter, line.c_str());
                return false;
            }
        }

        return true;
    };

    /**
     * @brief Checks the run parameter validity.
     * 
     * Checks the run parameters for validity. Issues some warnings for unconventional parameters. Creates the ouput directories.
     * 
     * @return True on success, False on failure.
     */
    bool checkParameters()
    {
        PRINT_LOG("Validating run parameters.", 2);

        int len = int(lst_matID.size());

        if (len == 0)
        {
            cout << "ERROR: No materials defined!  \n";
            return false;
        }

        if (N_save <= 0)
        {
            cout << "ERROR: Command \"<N_save>\" needs to be larger than zero!  \n";
            return false;
        }

        if (lst_matID[0] != 0)
        {
            cout << "ERROR: Lowest material ID has to be 1!  \n";
            return false;
        }

        for (int i = 0; i < len; i++)
        {
            if (lst_matID[i] != i)
            {
                cout << "ERROR: Material IDs not incremented by 1!  \n";
                return false;
            }
        }

        // Check each of the material parameters.
        for (int i = 0; i < len; i++)
        {
            if ((abs(lst_Msat[i]) == 0 && lst_Tc[i] > 0) || (abs(lst_Msat[i]) > 0 && lst_Tc[i] == 0))
            {
                cout << "ERROR: Spont. magnetization and Curie temperature connot both be zero!  \n";
                return false;
            }
        }

        if (mat_type == MAT_TYPE_NONE)
        {
            if (vec_lenght(Bext) > 0)
            {
                cout << "WARNING: None of the material can be magnetized but an external B-field is defined! \n";
            }
        }

        if (mat_type == MAT_TYPE_MAG)
        {
            if (vec_lenght(Bext) == 0)
            {
                cout << "WARNING: Some material can be magnetized but no external B-field is defined! \n";
            }
        }

        PRINT_LOG("Creating output directories.", 2);

        // Create the output directories.
        if (path_results.length() > 5)
        {
            if (!createPath(path_results))
                return false;

            if (!createPath(path_binary))
                return false;

            if (!createPath(path_ovito))
                return false;

            if (!createPath(path_plots))
                return false;
        }
        else
        {
            std::cout << "ERROR: Invalid path for simulation results!" << std::endl << std::flush;
            return false;
        }

        if (path_A.length() <= 5)
        {
            std::cout << "ERROR: Invalid path for aggregate A!" << std::endl << std::flush;
            return false;
        }

        if (path_B.length() <= 5)
        {
            std::cout << "ERROR: Invalid path for aggregate B!" << std::endl << std::flush;
            return false;
        }

        //TODO: Use the functions here?
        double len_pos_A = sqrt(pos_A.x * pos_A.x + pos_A.y * pos_A.y + pos_A.z * pos_A.z); // Distance of aggregate A from the origin.
        double len_pos_B = sqrt(pos_B.x * pos_B.x + pos_B.y * pos_B.y + pos_B.z * pos_B.z); // Distance of aggregate B from the origin.

        if (len_pos_A + len_pos_B == 0)
        {
            std::cout << "ERROR: Both aggregate postions are at the center! \n";
            return false;
        }

        if (T_dust <= 0 && T_dust != -1)
        {
            std::cout << "ERROR: Dust temperature needs to be larger than zero! \n";
            return false;
        }

        // TODO: Check if this is implemented.
        if (T_dust == -1)
        {
            cout << "WARNING: Impact of dust temperature is ignored! \n";
        }

        double len_vel_A = sqrt(vel_A.x * vel_A.x + vel_A.y * vel_A.y + vel_A.z * vel_A.z); // The initial absolute velocity of aggregate A.
        double len_vel_B = sqrt(vel_B.x * vel_B.x + vel_B.y * vel_B.y + vel_B.z * vel_B.z); // The initial absolute velocity of aggregate B.

        if (len_vel_A + len_vel_B == 0)
        {
            cout << "WARNING: Both aggregate velocities are zero! \n";
        }

        double len_ang_A = sqrt(ang_A.x * ang_A.x + ang_A.y * ang_A.y + ang_A.z * ang_A.z); // The initial absolute angular velocity of aggregate A.
        double len_ang_B = sqrt(ang_B.x * ang_B.x + ang_B.y * ang_B.y + ang_B.z * ang_B.z); // The initial absolute angular velocity of aggregate B.

        if (len_ang_A + len_ang_B == 0)
        {
            cout << "WARNING: Non rotating aggregates! \n";
        }

        if (N_iter <= 0)
        {
            PRINT_ERROR("Invalid number of iterations!");
            return false;
        }

        if (time_start < 0)
        {
            cout << "WARNING: Invalid starting time! \n";
        }

        if (time_stop <= 0)
        {
            cout << "WARNING: Invalid stopping time! \n";
        }

        // TODO: check time step?
        return true;
    };

    /**
     * @brief Uses a parameter name and value pair to set a run parameter.
     * 
     * @param cmd: The name of the parameter.
     * @param data: The value of the parameter.
     * @return True if successful, False if failure.
     */
    bool parseLine(string cmd, string data)
    {
        if (cmd.compare("<path_results>") == 0)
        {
            path_results = seperateString(data);

            path_binary = path_results + "binary" + SEP;
            path_ovito = path_results + "ovito" + SEP;
            path_plots = path_results + "plots" + SEP;

            return true;
        }

        if (cmd.compare("<path_A>") == 0)
        {
            path_A = seperateString(data);
            return true;
        }

        if (cmd.compare("<path_B>") == 0)
        {
            path_B = seperateString(data);
            return true;
        }

        if (cmd.compare("<pos_A>") == 0)
        {
            dlist values = parseValues(data);

            if (values.size() != 3)
            {
                PRINT_ERROR(std::string("Wrong ammount of coordinates in command \'") + cmd + "\'!\n");
                return false;
            }

            pos_A.x = values[0];
            pos_A.y = values[1];
            pos_A.z = values[2];

            return true;
        }

        if (cmd.compare("<pos_B>") == 0)
        {
            dlist values = parseValues(data);

            if (values.size() != 3)
            {
                PRINT_ERROR(std::string("Wrong ammount of coordinates in command \'") + cmd + "\'!\n");
                return false;
            }

            pos_B.x = values[0];
            pos_B.y = values[1];
            pos_B.z = values[2];

            return true;
        }

        if (cmd.compare("<B_ext>") == 0)
        {
            dlist values = parseValues(data);

            if (values.size() != 3)
            {
                PRINT_ERROR(std::string("Wrong ammount of coordinates in command \'") + cmd + "\'!\n");
                return false;
            }

            Bext.x = values[0];
            Bext.y = values[1];
            Bext.z = values[2];

            return true;
        }

        if (cmd.compare("<vel_A>") == 0)
        {
            dlist values = parseValues(data);

            if (values.size() != 3)
            {
                PRINT_ERROR(std::string("Wrong ammount of coordinates in command \'") + cmd + "\'!\n");
                return false;
            }

            vel_A.x = values[0];
            vel_A.y = values[1];
            vel_A.z = values[2];

            return true;
        }

        if (cmd.compare("<vel_B>") == 0)
        {
            dlist values = parseValues(data);

            if (values.size() != 3)
            {
                PRINT_ERROR(std::string("Wrong ammount of coordinates in command \'") + cmd + "\'!\n");
                return false;
            }

            vel_B.x = values[0];
            vel_B.y = values[1];
            vel_B.z = values[2];

            return true;
        }

        if (cmd.compare("<ang_A>") == 0)
        {
            dlist values = parseValues(data);

            if (values.size() != 3)
            {
                PRINT_ERROR(std::string("Wrong ammount of coordinates in command \'") + cmd + "\'!\n");
                return false;
            }

            ang_A.x = values[0];
            ang_A.y = values[1];
            ang_A.z = values[2];

            return true;
        }

        if (cmd.compare("<ang_B>") == 0)
        {
            dlist values = parseValues(data);

            if (values.size() != 3)
            {
                PRINT_ERROR(std::string("Wrong ammount of coordinates in command \'") + cmd + "\'!\n");
                return false;
            }

            ang_B.x = values[0];
            ang_B.y = values[1];
            ang_B.z = values[2];

            return true;
        }

        if (cmd.compare("<T_dust>") == 0)
        {
            double t = double(atof(data.c_str()));
            T_dust = t;

            return true;
        }

        if (cmd.compare("<time_start>") == 0)
        {
            double t = double(atof(data.c_str()));
            time_start = t;

            return true;
        }

        if (cmd.compare("<time_stop>") == 0)
        {
            double t = double(atof(data.c_str()));
            time_stop = t;

            return true;
        }

        if (cmd.compare("<time_step>") == 0)
        {
            double t = double(atof(data.c_str()));
            time_stop = t;

            return true;
        }

        if (cmd.compare("<N_iter>") == 0)
        {
            ullong i = ullong(atof(data.c_str()));
            N_iter = i;

            return true;
        }

        if (cmd.compare("<N_save>") == 0)
        {
            ullong i = ullong(atof(data.c_str()));
            N_save = i;

            return true;
        }

        if (cmd.compare("<save_pos>") == 0)
        {
            int s = int(atof(data.c_str()));

            if (s == 1)
                save_pos = true;

            return true;
        }

        if (cmd.compare("<save_vel>") == 0)
        {
            int s = int(atof(data.c_str()));

            if (s == 1)
                save_vel = true;

            return true;
        }

        if (cmd.compare("<save_force>") == 0)
        {
            int s = int(atof(data.c_str()));

            if (s == 1)
                save_force = true;

            return true;
        }

        if (cmd.compare("<save_torque>") == 0)
        {
            int s = int(atof(data.c_str()));

            if (s == 1)
                save_torque = true;

            return true;
        }

        if (cmd.compare("<save_cluster>") == 0)
        {
            int s = int(atof(data.c_str()));

            if (s == 1)
                save_cluster = true;

            return true;
        }

        if (cmd.compare("<save_omega>") == 0)
        {
            int s = int(atof(data.c_str()));

            if (s == 1)
                save_omega = true;

            return true;
        }

        if (cmd.compare("<save_mag>") == 0)
        {
            int s = int(atof(data.c_str()));

            if (s == 1)
                save_mag = true;

            return true;
        }

        if (cmd.compare("<save_ovito>") == 0)
        {
            int s = int(atof(data.c_str()));

            if (s == 1)
                save_ovito = true;

            return true;
        }

        if (cmd.compare("<material id = >") == 0)
        {
            string str_m = seperateString(data);
            string str_id = seperateString(data);
            int ID = int(atof(str_id.c_str()));

            dlist values = parseValues(data);
            int length = int(values.size());

            if (values.size() == 6)
            {
                lst_matName.push_back(str_m);
                lst_matID.push_back(ID - 1);

                lst_gamma.push_back(values[0]);
                lst_E.push_back(values[1]);
                lst_nu.push_back(values[2]);
                lst_rho.push_back(values[3]);
                lst_xi.push_back(values[4]);
                lst_Tvis.push_back(values[5]);

                lst_tss.push_back(0);
                lst_tsl.push_back(0);
                lst_Msat.push_back(0);
                lst_chi.push_back(0);
                lst_Tc.push_back(0);

                // TODO: Add warning that magnetic values have not been read

                return true;
            }

            if (values.size() == 11)
            {   
                // Magnetic materials are defined.
                mat_type = MAT_TYPE_MAG;

                lst_matName.push_back(str_m);
                lst_matID.push_back(ID - 1);

                lst_gamma.push_back(values[0]);
                lst_E.push_back(values[1]);
                lst_nu.push_back(values[2]);
                lst_rho.push_back(values[3]);
                lst_xi.push_back(values[4]);
                lst_Tvis.push_back(values[5]);

                lst_tss.push_back(values[6]);
                lst_tsl.push_back(values[7]);
                lst_Msat.push_back(values[8]);
                lst_chi.push_back(values[9]);
                lst_Tc.push_back(values[10]);

                return true;
            }

            PRINT_ERROR("Wrong ammount of material constants in command file!");
            return false;
        }

        return false;
    };

    /**
     * @brief Reads the aggregates from files and calculates the initial state.
     * 
     * @param *&pos: An array of the momomer positions will get placed here.
     * @param *&vel: An array of the monomer velocities will be placed here.
     * @param *&omega_tot: An array of the monomer angular velocities will be placed here.
     * @param *&mag: An array of the monomer magnetizations will be placed here.
     * @param *&amon: An array of the monomer radii will be placed here.
     * @param *&mass: An array of the monomer masses will be placed here.
     * @param *&moment: An array of the monomer moment of inertias will be placed here.
     * @param *&matIDs: An array of the momomer material IDs will be placed here.
     * @param *&Nmon: The total number of monomers will be placed here.
     * @returns True if successful, False if failure.
     */
    bool prepareData(
            vec3D*& pos,
            vec3D*& vel,
            vec3D*& omega_tot,
            vec3D*& mag,
            double*& amon,
            double*& mass,
            double*& moment,
            int*& matIDs,
            int& Nmon
        ) {
        // Initialize readers for aggregate files.
        ifstream reader_A, reader_B;

        vlist lst_pos_A, lst_pos_B;
        ilist lst_matID_A, lst_matID_B;
        dlist lst_amon_A, lst_amon_B;

        string line_A, line_B;
        int line_counter_A = 0, line_counter_B = 0;

        agg_filename_A = path_A;
        agg_filename_B = path_B;

        // Read the monomer data for aggregate A and aggregate B.
        reader_A.open(agg_filename_A.c_str());

        if (reader_A.fail())
        {
            PRINT_ERROR("Cannot read aggregate A from \n      " + agg_filename_A);
            return false;
        }

        reader_B.open(agg_filename_B.c_str());

        if (reader_B.fail())
        {
            PRINT_ERROR("Cannot read aggregate B from \n      " + agg_filename_B);
            return false;
        }

        // Initialize values for min and max value determination.
        a_mon_min = 1e200;
        a_mon_max = 0;

        // Read aggregate files.
        PRINT_LOG("Reading aggregate A from file.", 2);
        while (getline(reader_A, line_A))
        {
            line_counter_A++;

            // Use the first line to determine the number of monomers.
            if (line_counter_A == 1)
            {
                dlist values = parseValues(line_A);

                Nmon_A = int(values[0]);
                a_eff_A = 1e-9 * values[2];
            }

            // Read the monomer data from line 5 onward.
            if (line_counter_A > 5)
            {
                dlist values = parseValues(line_A);

                // Skip empty lines.
                if (values.size() == 0)
                    continue;
                
                if (line_counter_A % 50 == 0)
                    cout << "Reading aggregate A: " << 100.0 * float(line_counter_A) / float(Nmon_A) << "                 \r" << flush;

                // Read monomer positions from the first three numbers.
                vec3D tmp_pos;
                tmp_pos.x = 1e-9 * values[0];
                tmp_pos.y = 1e-9 * values[1];
                tmp_pos.z = 1e-9 * values[2];
                
                lst_pos_A.push_back(tmp_pos);

                // Read monomer radius from 5th number
                double a_mon = 1e-9 * values[4];
                lst_amon_A.push_back(a_mon);

                // Determine min and max monomer radii
                if (a_mon_min > a_mon)
                    a_mon_min = a_mon;

                if (a_mon_max < a_mon)
                    a_mon_max = a_mon;

                // Determine the maximum for the outer most surface from the origin of aggregate A
                double distance = sqrt(tmp_pos.x * tmp_pos.x + tmp_pos.y * tmp_pos.y + tmp_pos.z * tmp_pos.z) + a_mon;

                if (a_out_A < distance)
                    a_out_A = distance;

                // Read the material ID from the 7th number
                int mat_id = int(values[6] - 1);

                if (!isMatID(mat_id))
                {
                    cout << "ERROR: Material id " << mat_id << " in aggregate A line " << line_counter_A << " does not match the IDs in the command file!   \n" << flush;
                    return false;
                }

                lst_matID_A.push_back(mat_id);
            }
        }

        PRINT_LOG("Reading aggregate B from file.", 2);
        while (getline(reader_B, line_B))
        {
            line_counter_B++;

            if (line_counter_B == 1)
            {
                dlist values = parseValues(line_B);

                Nmon_B = int(values[0]);
                a_eff_B = 1e-9 * values[2];
            }

            if (line_counter_B > 5)
            {
                dlist values = parseValues(line_B);

                if (values.size() == 0)
                    continue;

                if (line_counter_B % 50 == 0)
                    cout << "Reading aggregate B: " << 100.0 * float(line_counter_B) / float(Nmon_B) << "                 \r" << flush;

                vec3D tmp_pos;

                tmp_pos.x = 1e-9 * values[0];
                tmp_pos.y = 1e-9 * values[1];
                tmp_pos.z = 1e-9 * values[2];
                lst_pos_B.push_back(tmp_pos);

                double a_mon = 1e-9 * values[4];
                lst_amon_B.push_back(a_mon);

                if (a_mon_min > a_mon)
                    a_mon_min = a_mon;

                if (a_mon_max < a_mon)
                    a_mon_max = a_mon;

                double distance = sqrt(tmp_pos.x * tmp_pos.x + tmp_pos.y * tmp_pos.y + tmp_pos.z * tmp_pos.z) + a_mon;

                if (a_out_B < distance)
                    a_out_B = distance;

                int mat_id = int(values[6] - 1);

                if (!isMatID(mat_id))
                {
                    cout << "ERROR: Material id " << mat_id << " in aggregate B line " << line_counter_B << " does not match the IDs in the command file!   \n" << flush;
                    return false;
                }

                lst_matID_B.push_back(mat_id);
            }
        }

        Nmon = Nmon_A + Nmon_B;

        // Create arrays to contain position, velocity, magnetization, angular velocities, material IDs, monomer radii
        PRINT_LOG("Calculating initial state.", 2);
        
        pos = new vec3D[Nmon];
        vel = new vec3D[Nmon];
        mag = new vec3D[Nmon];
        omega_tot = new vec3D[Nmon];
        matIDs = new int[Nmon];
        amon = new double[Nmon];

        moment = new double[Nmon];
        mass = new double[Nmon];

        // Iterate over the 
        for (int i = 0; i < Nmon_A; i++)
        {
            // Determine monomer position in experimental frame.
            pos[i].x = lst_pos_A[i].x + pos_A.x;
            pos[i].y = lst_pos_A[i].y + pos_A.y;
            pos[i].z = lst_pos_A[i].z + pos_A.z;

            amon[i] = lst_amon_A[i];

            // Determine the velocity of the monomers due to aggregate rotation.
            vec3D r;
            r.x = 0; // FIXME: This should be a bug?
            r.y = lst_pos_A[i].y;
            r.z = lst_pos_A[i].z;

            vec3D vel_tan = cpu_vec3D_cross(ang_A, r);

            vel[i].x = vel_A.x + vel_tan.x;
            vel[i].y = vel_A.y + vel_tan.y;
            vel[i].z = vel_A.z + vel_tan.z;

            // The angular momentum of the monomers is equal to the aggregates angular momentum.
            omega_tot[i] = ang_A;

            // FIXME: Is mat_id not redundant?
            int mat_id = lst_matID_A[i];
            matIDs[i] = mat_id;

            double rho = lst_rho[mat_id];

            // Calculate moment of inertia and mass of the monomer.
            mass[i] = 4. / 3. * PI * rho * amon[i] * amon[i] * amon[i];
            moment[i] = 2. / 5. * mass[i] * amon[i] * amon[i];

            // TODO: Code from the physics engine should be used here...
            // The magnetic susceptibility.
            double chi = lst_chi[mat_id];

            // Set initial magnetization depending on run setup
            if (abs(chi) != 0)
            {
                if (abs(chi) > LIMIT_FER)
                {
                    // Ferromagnetic materials
                    double len_Bext = vec_lenght(Bext);

                    if (len_Bext > 0.)
                    {
                        // Saturation magnetization in direction of external field.
                        // FIXME: Magnetization initialization
                        //mag[i].x = lst_Msat[mat_id] * Bext.x / len_Bext;
                        //mag[i].y = lst_Msat[mat_id] * Bext.y / len_Bext;
                        //mag[i].z = lst_Msat[mat_id] * Bext.z / len_Bext;

                        mag[i].x = 0.984807753 * lst_Msat[mat_id];
                        mag[i].y = 0;
                        mag[i].z = 0.173648178 * lst_Msat[mat_id];
                    }
                    else
                    {   
                        double len_omega = cpu_vec3D_length(omega_tot[i]);
                        if (len_omega > 0.)
                        {
                            // Saturation magnetization in direction of angular momentum.
                            mag[i].x = lst_Msat[mat_id] * omega_tot[i].x / len_omega;
                            mag[i].y = lst_Msat[mat_id] * omega_tot[i].y / len_omega;
                            mag[i].z = lst_Msat[mat_id] * omega_tot[i].z / len_omega;
                        }
                        else
                        {
                            // Saturation magnetization in direction of x axis.
                            mag[i].x = lst_Msat[mat_id];
                            mag[i].y = 0.;
                            mag[i].z = 0.;
                        }
                    }
                }
                else
                {
                    // Para + Diamagnetic materials
                    vec3D M_ind, M_Bar;
                    double chi_fac = chi / (chi + 1.);

                    M_ind.x = chi_fac * Bext.x / mu0;
                    M_ind.y = chi_fac * Bext.y / mu0;
                    M_ind.z = chi_fac * Bext.z / mu0;

                    M_Bar.x = chi * omega_tot[i].x / PROD_BARR;
                    M_Bar.y = chi * omega_tot[i].y / PROD_BARR;
                    M_Bar.z = chi * omega_tot[i].z / PROD_BARR;

                    // Set the magnetization of the monomer to the induced and Barnett moment.
                    mag[i].x = M_ind.x + M_Bar.x;
                    mag[i].y = M_ind.y + M_Bar.y;
                    mag[i].z = M_ind.z + M_Bar.z;

                    double len_mag = cpu_vec3D_length(mag[i]);

                    // Clip the magnetization to saturation.
                    if (len_mag > lst_Msat[mat_id])
                    {
                        mag[i].x = lst_Msat[mat_id] * mag[i].x / len_mag;
                        mag[i].y = lst_Msat[mat_id] * mag[i].y / len_mag;
                        mag[i].z = lst_Msat[mat_id] * mag[i].z / len_mag;
                    }

                    // FIXME: This needs to go.
                    mag[i].x = 0.984807753 * lst_Msat[mat_id];
                    mag[i].y = 0;
                    mag[i].z = 0.173648178 * lst_Msat[mat_id];
                }
            }
            else
            {
                // Nonmagnetic materials.
                mag[i].x = 0;
                mag[i].y = 0;
                mag[i].z = 0;
            }
        }

        // FIXME: This should just be a single loop, the logic should not be duplicated.
        // Iterate over the monomers of aggregate B
        for (int i = Nmon_A; i < Nmon; i++)
        {
            int index = i - Nmon_A;
            pos[i].x = lst_pos_B[index].x + pos_B.x;
            pos[i].y = lst_pos_B[index].y + pos_B.y;
            pos[i].z = lst_pos_B[index].z + pos_B.z;

            amon[i] = lst_amon_B[index];

            vec3D r;

            r.x = 0;
            r.y = lst_pos_B[index].y;
            r.z = lst_pos_B[index].z;

            vec3D vel_tan = cpu_vec3D_cross(ang_B, r);

            vel[i].x = vel_B.x + vel_tan.x;
            vel[i].y = vel_B.y + vel_tan.y;
            vel[i].z = vel_B.z + vel_tan.z;

            omega_tot[i] = ang_B;

            int mat_id = lst_matID_B[index];
            matIDs[i] = mat_id;

            double rho = lst_rho[mat_id];

            mass[i] = 4. / 3. * PI * rho * amon[i] * amon[i] * amon[i];
            moment[i] = 2. / 5. * mass[i] * amon[i] * amon[i];

            double chi = lst_chi[mat_id];

            if (abs(chi) != 0)
            {
                if (abs(chi) > LIMIT_FER) //ferromagnetic
                {
                    double len_Bext = vec_lenght(Bext);

                    if (len_Bext > 0.) //set initial direction to Bext
                    {
                        //mag[i].x = lst_Msat[mat_id] * Bext.x / len_Bext;
                        //mag[i].y = lst_Msat[mat_id] * Bext.y / len_Bext;
                        //mag[i].z = lst_Msat[mat_id] * Bext.z / len_Bext;

                        // FIXME: This needs to go
                        mag[i].x = 0.984807753 * lst_Msat[mat_id];
                        mag[i].y = 0;
                        mag[i].z = 0.173648178 * lst_Msat[mat_id];
                    }
                    else
                    {
                        double len_omega = cpu_vec3D_length(omega_tot[i]);
                        if (len_omega > 0.) //set initial direction to omega
                        {
                            mag[i].x = lst_Msat[mat_id] * omega_tot[i].x / len_omega;
                            mag[i].y = lst_Msat[mat_id] * omega_tot[i].y / len_omega;
                            mag[i].z = lst_Msat[mat_id] * omega_tot[i].z / len_omega;
                        }
                        else //set initial direction to x-direction
                        {
                            mag[i].x = lst_Msat[mat_id];
                            mag[i].y = 0.;
                            mag[i].z = 0.;
                        }
                    }
                }
                else
                {
                    vec3D M_ind, M_Bar;
                    double chi_fac = chi / (chi + 1.);

                    M_ind.x = chi_fac * Bext.x / mu0;
                    M_ind.y = chi_fac * Bext.y / mu0;
                    M_ind.z = chi_fac * Bext.z / mu0;

                    M_Bar.x = chi * omega_tot[i].x / PROD_BARR;
                    M_Bar.y = chi * omega_tot[i].y / PROD_BARR;
                    M_Bar.z = chi * omega_tot[i].z / PROD_BARR;

                    mag[i].x = M_ind.x + M_Bar.x;
                    mag[i].y = M_ind.y + M_Bar.y;
                    mag[i].z = M_ind.z + M_Bar.z;

                    double len_mag = cpu_vec3D_length(mag[i]);

                    if (len_mag > lst_Msat[mat_id])
                    {
                        mag[i].x = lst_Msat[mat_id] * mag[i].x / len_mag;
                        mag[i].y = lst_Msat[mat_id] * mag[i].y / len_mag;
                        mag[i].z = lst_Msat[mat_id] * mag[i].z / len_mag;
                    }

                    mag[i].x = 0.984807753 * lst_Msat[mat_id];
                    mag[i].y = 0;
                    mag[i].z = 0.173648178 * lst_Msat[mat_id];
                }
            }
            else
            {
                mag[i].x = 0;
                mag[i].y = 0;
                mag[i].z = 0;
            }
        }

        // Find the smallest timestep.
        PRINT_LOG("Determining simulation timestep.", 2);
        if (time_step == 0)
        {
            time_step = 1e200;

            // Iterate over monomer pairs.
            for (int i = 0; i < Nmon; i++)
            {
                for (int j = 0; j < Nmon; j++)
                {
                    if (i == j)
                        continue;

                    int mat_id_A = matIDs[i];
                    int mat_id_B = matIDs[j];

                    double mass_A = mass[i];
                    double mass_B = mass[j];
                    double mass = (mass_A * mass_B) / (mass_A + mass_B);

                    double a_mon_A = amon[i];
                    double a_mon_B = amon[j];
                    double R = (a_mon_A * a_mon_B) / (a_mon_A + a_mon_B);

                    double nu_A = lst_nu[mat_id_A];
                    double nu_B = lst_nu[mat_id_B];

                    double E_A = lst_E[mat_id_A];
                    double E_B = lst_E[mat_id_B];

                    double Es = (1 - nu_A * nu_A) / E_A + (1 - nu_B * nu_B) / E_B;
                    Es = 1. / Es;

                    double gamma_A = lst_gamma[mat_id_A];
                    double gamma_B = lst_gamma[mat_id_B];
                    double gamma = gamma_A + gamma_B - 2. / (1. / gamma_A + 1. / gamma_B);

                    double r0 = pow(9 * PI * gamma * R * R / Es, 1. / 3.);
                    double delta_c = 0.5 * r0 * r0 / (R * pow(6., 1. / 3.));
                    double F_c = 3.0 * PI * gamma * R;

                    double tc = sqrt(mass * delta_c / F_c);

                    // TODO: Should this stay commented out?
                    //double tc = 0.95 * (pow(R, 7. / 6.) * sqrt(lst_rho[mat_id_A])) / (pow(gamma, 1. / 6.) * pow(Es, 1. / 3.));

                    if (time_step > tc)
                        time_step = tc;
                }
            }

            int Nmat = int(lst_matName.size());

            for (int i = 0; i < Nmat; i++)
            {
                double tss = lst_tss[i];
                double tsl = lst_tsl[i];

                if (tss > 0)
                {
                    if (time_step > tss)
                        time_step = tss;
                }

                if (tsl > 0)
                {
                    if (time_step > tsl)
                        time_step = tsl;
                }
            }

            // Set the timestep to 1/200 of the minimum physical timescale of the system.
            time_step = 0.005 * time_step;
        }

        PRINT_CLR_LINE();

        return true;
    };

    bool prepareData_(
            double3*& pos,
            double3*& vel,
            double3*& omega_tot,
            double3*& mag,
            double*& amon,
            double*& mass,
            double*& moment,
            int*& matIDs,
            int& Nmon
        ) {
        // Initialize readers for aggregate files.
        ifstream reader_A, reader_B;

        std::vector<double3> lst_pos_A, lst_pos_B;
        std::vector<int> lst_matID_A, lst_matID_B;
        std::vector<double> lst_amon_A, lst_amon_B;

        string line_A, line_B;
        int line_counter_A = 0, line_counter_B = 0;

        agg_filename_A = path_A;
        agg_filename_B = path_B;

        // Read the monomer data for aggregate A and aggregate B.
        reader_A.open(agg_filename_A.c_str());

        if (reader_A.fail())
        {
            PRINT_ERROR("Cannot read aggregate A from \n      " + agg_filename_A);
            return false;
        }

        reader_B.open(agg_filename_B.c_str());

        if (reader_B.fail())
        {
            PRINT_ERROR("Cannot read aggregate B from \n      " + agg_filename_B);
            return false;
        }

        // Initialize values for min and max value determination.
        a_mon_min = 1e200;
        a_mon_max = 0;

        // Read aggregate files.
        PRINT_LOG("Reading aggregate A from file.", 2);
        while (getline(reader_A, line_A))
        {
            line_counter_A++;

            // Use the first line to determine the number of monomers.
            if (line_counter_A == 1)
            {
                dlist values = parseValues(line_A);

                Nmon_A = int(values[0]);
                a_eff_A = 1e-9 * values[2];
            }

            // Read the monomer data from line 5 onward.
            if (line_counter_A > 5)
            {
                dlist values = parseValues(line_A);

                // Skip empty lines.
                if (values.size() == 0)
                    continue;
                
                if (line_counter_A % 50 == 0)
                    cout << "Reading aggregate A: " << 100.0 * float(line_counter_A) / float(Nmon_A) << "                 \r" << flush;

                // Read monomer positions from the first three numbers.
                double3 tmp_pos;
                tmp_pos.x = 1e-9 * values[0];
                tmp_pos.y = 1e-9 * values[1];
                tmp_pos.z = 1e-9 * values[2];
                
                lst_pos_A.push_back(tmp_pos);

                // Read monomer radius from 5th number
                double a_mon = 1e-9 * values[4];
                lst_amon_A.push_back(a_mon);

                // Determine min and max monomer radii
                if (a_mon_min > a_mon)
                    a_mon_min = a_mon;

                if (a_mon_max < a_mon)
                    a_mon_max = a_mon;

                // Determine the maximum for the outer most surface from the origin of aggregate A
                double distance = sqrt(tmp_pos.x * tmp_pos.x + tmp_pos.y * tmp_pos.y + tmp_pos.z * tmp_pos.z) + a_mon;

                if (a_out_A < distance)
                    a_out_A = distance;

                // Read the material ID from the 7th number
                int mat_id = int(values[6] - 1);

                if (!isMatID(mat_id))
                {
                    cout << "ERROR: Material id " << mat_id << " in aggregate A line " << line_counter_A << " does not match the IDs in the command file!   \n" << flush;
                    return false;
                }

                lst_matID_A.push_back(mat_id);
            }
        }

        PRINT_LOG("Reading aggregate B from file.", 2);
        while (getline(reader_B, line_B))
        {
            line_counter_B++;

            if (line_counter_B == 1)
            {
                dlist values = parseValues(line_B);

                Nmon_B = int(values[0]);
                a_eff_B = 1e-9 * values[2];
            }

            if (line_counter_B > 5)
            {
                dlist values = parseValues(line_B);

                if (values.size() == 0)
                    continue;

                if (line_counter_B % 50 == 0)
                    cout << "Reading aggregate B: " << 100.0 * float(line_counter_B) / float(Nmon_B) << "                 \r" << flush;

                double3 tmp_pos;

                tmp_pos.x = 1e-9 * values[0];
                tmp_pos.y = 1e-9 * values[1];
                tmp_pos.z = 1e-9 * values[2];
                lst_pos_B.push_back(tmp_pos);

                double a_mon = 1e-9 * values[4];
                lst_amon_B.push_back(a_mon);

                if (a_mon_min > a_mon)
                    a_mon_min = a_mon;

                if (a_mon_max < a_mon)
                    a_mon_max = a_mon;

                double distance = sqrt(tmp_pos.x * tmp_pos.x + tmp_pos.y * tmp_pos.y + tmp_pos.z * tmp_pos.z) + a_mon;

                if (a_out_B < distance)
                    a_out_B = distance;

                int mat_id = int(values[6] - 1);

                if (!isMatID(mat_id))
                {
                    cout << "ERROR: Material id " << mat_id << " in aggregate B line " << line_counter_B << " does not match the IDs in the command file!   \n" << flush;
                    return false;
                }

                lst_matID_B.push_back(mat_id);
            }
        }

        Nmon = Nmon_A + Nmon_B;

        // Create arrays to contain initial position, velocity, magnetization, angular velocities, material IDs, monomer radii
        PRINT_LOG("Calculating initial state.", 2);
        
        pos = (double3*) malloc(Nmon * sizeof(double3));
        vel = (double3*) malloc(Nmon * sizeof(double3));
        mag = (double3*) malloc(Nmon * sizeof(double3));
        omega_tot = (double3*) malloc(Nmon * sizeof(double3));
        matIDs = (int*) malloc(Nmon * sizeof(int));
        amon = (double*) malloc(Nmon * sizeof(double));
        moment = (double*) malloc(Nmon * sizeof(double));
        mass = (double*) malloc(Nmon * sizeof(double));

        // Iterate over the monomers in cluster A
        for (int i = 0; i < Nmon_A; i++)
        {
            // Determine monomer position in experimental frame.
            pos[i].x = lst_pos_A[i].x + pos_A.x;
            pos[i].y = lst_pos_A[i].y + pos_A.y;
            pos[i].z = lst_pos_A[i].z + pos_A.z;

            amon[i] = lst_amon_A[i];

            // Determine the velocity of the monomers due to aggregate rotation.
            double3 r;
            r.x = 0; // FIXME: This should be a bug?
            r.y = lst_pos_A[i].y;
            r.z = lst_pos_A[i].z;

            // TODO: Change CPipeline::ang_A to double3, or better: why ist ang_A even a class field?
            double3 ang_A_ = make_double3(ang_A.x, ang_A.y, ang_B.z);

            double3 vel_tan = vec_cross(ang_A_, r);

            vel[i].x = vel_A.x + vel_tan.x;
            vel[i].y = vel_A.y + vel_tan.y;
            vel[i].z = vel_A.z + vel_tan.z;

            // The angular momentum of the monomers is equal to the aggregates angular momentum.
            omega_tot[i] = ang_A_;

            // FIXME: Is mat_id not redundant?
            int mat_id = lst_matID_A[i];
            matIDs[i] = mat_id;

            double rho = lst_rho[mat_id];

            // Calculate moment of inertia and mass of the monomer.
            mass[i] = 4. / 3. * PI * rho * amon[i] * amon[i] * amon[i];
            moment[i] = 2. / 5. * mass[i] * amon[i] * amon[i];

            // TODO: Code from the physics engine should be used here...
            // The magnetic susceptibility.
            double chi = lst_chi[mat_id];

            // Set initial magnetization depending on run setup
            if (abs(chi) != 0)
            {
                if (abs(chi) > LIMIT_FER)
                {
                    // Ferromagnetic materials
                    double len_Bext = vec_lenght(Bext);

                    if (len_Bext > 0.)
                    {
                        // Saturation magnetization in direction of external field.
                        // FIXME: Magnetization initialization
                        //mag[i].x = lst_Msat[mat_id] * Bext.x / len_Bext;
                        //mag[i].y = lst_Msat[mat_id] * Bext.y / len_Bext;
                        //mag[i].z = lst_Msat[mat_id] * Bext.z / len_Bext;

                        mag[i].x = 0.984807753 * lst_Msat[mat_id];
                        mag[i].y = 0;
                        mag[i].z = 0.173648178 * lst_Msat[mat_id];
                    }
                    else
                    {   
                        double len_omega = vec_lenght(omega_tot[i]);
                        if (len_omega > 0.)
                        {
                            // Saturation magnetization in direction of angular momentum.
                            mag[i].x = lst_Msat[mat_id] * omega_tot[i].x / len_omega;
                            mag[i].y = lst_Msat[mat_id] * omega_tot[i].y / len_omega;
                            mag[i].z = lst_Msat[mat_id] * omega_tot[i].z / len_omega;
                        }
                        else
                        {
                            // Saturation magnetization in direction of x axis.
                            mag[i].x = lst_Msat[mat_id];
                            mag[i].y = 0.;
                            mag[i].z = 0.;
                        }
                    }
                }
                else
                {
                    // Para + Diamagnetic materials
                    vec3D M_ind, M_Bar;
                    double chi_fac = chi / (chi + 1.);

                    M_ind.x = chi_fac * Bext.x / mu0;
                    M_ind.y = chi_fac * Bext.y / mu0;
                    M_ind.z = chi_fac * Bext.z / mu0;

                    M_Bar.x = chi * omega_tot[i].x / PROD_BARR;
                    M_Bar.y = chi * omega_tot[i].y / PROD_BARR;
                    M_Bar.z = chi * omega_tot[i].z / PROD_BARR;

                    // Set the magnetization of the monomer to the induced and Barnett moment.
                    mag[i].x = M_ind.x + M_Bar.x;
                    mag[i].y = M_ind.y + M_Bar.y;
                    mag[i].z = M_ind.z + M_Bar.z;

                    double len_mag = vec_lenght(mag[i]);

                    // Clip the magnetization to saturation.
                    if (len_mag > lst_Msat[mat_id])
                    {
                        mag[i].x = lst_Msat[mat_id] * mag[i].x / len_mag;
                        mag[i].y = lst_Msat[mat_id] * mag[i].y / len_mag;
                        mag[i].z = lst_Msat[mat_id] * mag[i].z / len_mag;
                    }

                    // FIXME: This needs to go.
                    mag[i].x = 0.984807753 * lst_Msat[mat_id];
                    mag[i].y = 0;
                    mag[i].z = 0.173648178 * lst_Msat[mat_id];
                }
            }
            else
            {
                // Nonmagnetic materials.
                mag[i].x = 0;
                mag[i].y = 0;
                mag[i].z = 0;
            }
        }

        // FIXME: This should just be a single loop, the logic should not be duplicated.
        // Iterate over the monomers of aggregate B
        for (int i = Nmon_A; i < Nmon; i++)
        {
            int index = i - Nmon_A;
            pos[i].x = lst_pos_B[index].x + pos_B.x;
            pos[i].y = lst_pos_B[index].y + pos_B.y;
            pos[i].z = lst_pos_B[index].z + pos_B.z;

            amon[i] = lst_amon_B[index];

            vec3D r;

            r.x = 0;
            r.y = lst_pos_B[index].y;
            r.z = lst_pos_B[index].z;

            vec3D vel_tan = cpu_vec3D_cross(ang_B, r);

            vel[i].x = vel_B.x + vel_tan.x;
            vel[i].y = vel_B.y + vel_tan.y;
            vel[i].z = vel_B.z + vel_tan.z;

            // TODO: Change CPipeline::ang_A to double3, or better: why ist ang_A even a class field?
            double3 ang_B_ = make_double3(ang_B.x, ang_B.y, ang_B.z);

            omega_tot[i] = ang_B_;

            int mat_id = lst_matID_B[index];
            matIDs[i] = mat_id;

            double rho = lst_rho[mat_id];

            mass[i] = 4. / 3. * PI * rho * amon[i] * amon[i] * amon[i];
            moment[i] = 2. / 5. * mass[i] * amon[i] * amon[i];

            double chi = lst_chi[mat_id];

            if (abs(chi) != 0)
            {
                if (abs(chi) > LIMIT_FER) //ferromagnetic
                {
                    double len_Bext = vec_lenght(Bext);

                    if (len_Bext > 0.) //set initial direction to Bext
                    {
                        //mag[i].x = lst_Msat[mat_id] * Bext.x / len_Bext;
                        //mag[i].y = lst_Msat[mat_id] * Bext.y / len_Bext;
                        //mag[i].z = lst_Msat[mat_id] * Bext.z / len_Bext;

                        // FIXME: This needs to go
                        mag[i].x = 0.984807753 * lst_Msat[mat_id];
                        mag[i].y = 0;
                        mag[i].z = 0.173648178 * lst_Msat[mat_id];
                    }
                    else
                    {
                        double len_omega = vec_lenght(omega_tot[i]);
                        if (len_omega > 0.) //set initial direction to omega
                        {
                            mag[i].x = lst_Msat[mat_id] * omega_tot[i].x / len_omega;
                            mag[i].y = lst_Msat[mat_id] * omega_tot[i].y / len_omega;
                            mag[i].z = lst_Msat[mat_id] * omega_tot[i].z / len_omega;
                        }
                        else //set initial direction to x-direction
                        {
                            mag[i].x = lst_Msat[mat_id];
                            mag[i].y = 0.;
                            mag[i].z = 0.;
                        }
                    }
                }
                else
                {
                    vec3D M_ind, M_Bar;
                    double chi_fac = chi / (chi + 1.);

                    M_ind.x = chi_fac * Bext.x / mu0;
                    M_ind.y = chi_fac * Bext.y / mu0;
                    M_ind.z = chi_fac * Bext.z / mu0;

                    M_Bar.x = chi * omega_tot[i].x / PROD_BARR;
                    M_Bar.y = chi * omega_tot[i].y / PROD_BARR;
                    M_Bar.z = chi * omega_tot[i].z / PROD_BARR;

                    mag[i].x = M_ind.x + M_Bar.x;
                    mag[i].y = M_ind.y + M_Bar.y;
                    mag[i].z = M_ind.z + M_Bar.z;

                    double len_mag = vec_lenght(mag[i]);

                    if (len_mag > lst_Msat[mat_id])
                    {
                        mag[i].x = lst_Msat[mat_id] * mag[i].x / len_mag;
                        mag[i].y = lst_Msat[mat_id] * mag[i].y / len_mag;
                        mag[i].z = lst_Msat[mat_id] * mag[i].z / len_mag;
                    }

                    mag[i].x = 0.984807753 * lst_Msat[mat_id];
                    mag[i].y = 0;
                    mag[i].z = 0.173648178 * lst_Msat[mat_id];
                }
            }
            else
            {
                mag[i].x = 0;
                mag[i].y = 0;
                mag[i].z = 0;
            }
        }

        // Find the smallest timestep.
        PRINT_LOG("Determining simulation timestep.", 2);
        if (time_step == 0)
        {
            time_step = 1e200;

            // Iterate over monomer pairs.
            for (int i = 0; i < Nmon; i++)
            {
                for (int j = 0; j < Nmon; j++)
                {
                    if (i == j)
                        continue;

                    int mat_id_A = matIDs[i];
                    int mat_id_B = matIDs[j];

                    double mass_A = mass[i];
                    double mass_B = mass[j];
                    double mass = (mass_A * mass_B) / (mass_A + mass_B);

                    double a_mon_A = amon[i];
                    double a_mon_B = amon[j];
                    double R = (a_mon_A * a_mon_B) / (a_mon_A + a_mon_B);

                    double nu_A = lst_nu[mat_id_A];
                    double nu_B = lst_nu[mat_id_B];

                    double E_A = lst_E[mat_id_A];
                    double E_B = lst_E[mat_id_B];

                    double Es = (1 - nu_A * nu_A) / E_A + (1 - nu_B * nu_B) / E_B;
                    Es = 1. / Es;

                    double gamma_A = lst_gamma[mat_id_A];
                    double gamma_B = lst_gamma[mat_id_B];
                    double gamma = gamma_A + gamma_B - 2. / (1. / gamma_A + 1. / gamma_B);

                    double r0 = pow(9 * PI * gamma * R * R / Es, 1. / 3.);
                    double delta_c = 0.5 * r0 * r0 / (R * pow(6., 1. / 3.));
                    double F_c = 3.0 * PI * gamma * R;

                    double tc = sqrt(mass * delta_c / F_c);

                    // TODO: Should this stay commented out?
                    //double tc = 0.95 * (pow(R, 7. / 6.) * sqrt(lst_rho[mat_id_A])) / (pow(gamma, 1. / 6.) * pow(Es, 1. / 3.));

                    if (time_step > tc)
                        time_step = tc;
                }
            }

            int Nmat = int(lst_matName.size());

            for (int i = 0; i < Nmat; i++)
            {
                double tss = lst_tss[i];
                double tsl = lst_tsl[i];

                if (tss > 0)
                {
                    if (time_step > tss)
                        time_step = tss;
                }

                if (tsl > 0)
                {
                    if (time_step > tsl)
                        time_step = tsl;
                }
            }

            // Set the timestep to 1/200 of the minimum physical timescale of the system.
            time_step = 0.005 * time_step;
        }

        PRINT_CLR_LINE();

        return true;
    };
    
    /**
     * @brief Processes the material parameters and places them in arrays.
     * 
     * @param *&mat An array of materials will be placed here.
     * @param &int The ammount of materials will be placed here.
     */
    void prepareMaterial(material*& mat, int& Nmat)
    {
        PRINT_LOG("Preparing the material parameters.", 2);

        Nmat = int(lst_matName.size());
        mat = new material[Nmat];

        for (int i = 0; i < Nmat; i++)
        {
            // Transfer the material parameters into the material struct.
            mat[i].gamma = lst_gamma[i];
            mat[i].E = lst_E[i];
            mat[i].nu = lst_nu[i];
            mat[i].rho = lst_rho[i];
            mat[i].xi = lst_xi[i];
            mat[i].tvis = lst_Tvis[i];

            mat[i].tss = lst_tss[i];
            mat[i].tsl = lst_tsl[i];
            mat[i].Msat = lst_Msat[i];
            mat[i].chi = lst_chi[i];
            mat[i].Tc = lst_Tc[i];

            // Adjust material parameters based on dust temperature.
            if (T_dust != -1) // When T_dust is -1, dust temperature impact is ignored.
            {
                double corr;

                // Curie's law
                corr = CHI_20 / T_dust;
                //mat[i].chi *= corr; // FIXME: Turn on the temperature impact.

                if (abs(mat[i].Msat) > 0 && mat[i].Tc > 0)
                {
                    //Bloch T^3/2 law
                    corr = 1.0 - pow(T_dust / mat[i].Tc, 1.5);

                    if (corr < 0.0)
                        corr = 0.0;

                    //mat[i].Msat *= corr; // FIXME: Turn on the temperature impact
                }

                // reduced surface energy Bogdan+ 2020
                corr = SLOPE * T_dust + INTERCEPT;
                //mat[i].gamma *= corr; // FIXME: Turn on the energy correction!
            }

            // Determine minimum and maximum values for the material parameters.
            if (min_gamma > mat[i].gamma)
                min_gamma = mat[i].gamma;

            if (max_gamma < mat[i].gamma)
                max_gamma = mat[i].gamma;

            if (min_E > mat[i].E)
                min_E = mat[i].E;

            if (max_E < mat[i].E)
                max_E = mat[i].E;

            if (min_nu > mat[i].nu)
                min_nu = mat[i].nu;

            if (max_nu < mat[i].nu)
                max_nu = mat[i].nu;

            if (min_rho > mat[i].rho)
                min_rho = mat[i].rho;

            if (max_rho < mat[i].rho)
                max_rho = mat[i].rho;

            if (min_xi > mat[i].xi)
                min_xi = mat[i].xi;

            if (max_xi < mat[i].xi)
                max_xi = mat[i].xi;

            if (min_Tvis > mat[i].tvis)
                min_Tvis = mat[i].tvis;

            if (max_Tvis < mat[i].tvis)
                max_Tvis = mat[i].tvis;

            if (mat[i].tss != 0)
            {
                if (min_tss > mat[i].tss)
                    min_tss = mat[i].tss;

                if (max_tss < mat[i].tss)
                    max_tss = mat[i].tss;
            }

            if (mat[i].tsl != 0)
            {
                if (min_tsl > mat[i].tsl)
                    min_tsl = mat[i].tsl;

                if (max_tsl < mat[i].tsl)
                    max_tsl = mat[i].tsl;
            }

            if (mat[i].Msat != 0)
            {
                if (min_Msat > mat[i].Msat)
                    min_Msat = mat[i].Msat;

                if (max_Msat < mat[i].Msat)
                    max_Msat = mat[i].Msat;
            }

            if (mat[i].chi != 0)
            {
                if (min_chi > mat[i].chi)
                    min_chi = mat[i].chi;

                if (max_chi < mat[i].chi)
                    max_chi = mat[i].chi;
            }

            if (mat[i].Tc != 0)
            {
                if (min_Tc > mat[i].Tc)
                    min_Tc = mat[i].Tc;

                if (max_Tc < mat[i].Tc)
                    max_Tc = mat[i].Tc;
            }
        }
    };

    /**
     * @brief Checks if a material of the specified ID exists.
     * 
     * @param id: The material ID that is to be checked.
     */
    bool isMatID(int id)
    {
        int len = int(lst_matID.size());

        for (int i = 0; i < len; i++)
        {
            if (lst_matID[i] == i)
                return true;
        }

        return false;
    };

    // TODO: Make Nmon and Nsto consants.
    /**
     * @brief Writes the provided data into an ovito dump file.
     * 
     * @param *pos: Monomer positions.
     * @param *vel: Monomer velocities.
     * @param *force: Total forces on the monomers.
     * @param *torque: Total torques on the monomers.
     * @param *omega: Monomer angular momenta.
     * @param *mag: Monomer magnetizations.
     * @param *cluserIDs: Monomer cluster IDs.
     * @param *amon: Monomer radii.
     * @param *matID: Monomer material IDs.
     * @param *Nmon: Number of monomers.
     * @param *Nsto: // TODO: determine what this is.
     */
    bool writeAllOVITO(const double3* pos, const double3* vel, const double3* force, const double3* torque, const double3* omega, const double3* mag, const int* cluserIDs, const double* amon, const int* matID, int Nmon, int Nsto)
    {
        char str_tmp[1024];
        char str_end[1024];

        double b_size = 0;

        for (int i = 0; i < Nsto; i++)
        {
            int index = i % Nmon;
            double a = amon[index];
            double distance = vec_lenght(pos[i]) + a;

            if (b_size < distance)
                b_size = distance;
        }

        b_size *= 1.15e9;

        if (b_size > 10000)
            b_size = 10000;

        int steps = Nsto / Nmon;

        double max_vel = 0;
        double max_force = 0;
        double max_torque = 0;
        double max_omega = 0;
        double max_mag = 0;

        for (int i = 0; i < Nmon; i++)
        {
            if (vel != 0)
            {
                double len_vel = vec_lenght(vel[i]);
                if (max_vel < len_vel)
                    max_vel = sqrt(len_vel);
            }

            if (force != 0)
            {
                double len_force = vec_lenght(force[i]);
                if (max_force < len_force)
                    max_force = sqrt(len_force);
            }

            if (torque != 0)
            {
                double len_torque = vec_lenght(torque[i]);
                if (max_torque < len_torque)
                    max_torque = sqrt(len_torque);
            }

            if (omega != 0)
            {
                double len_omega = vec_lenght(omega[i]);
                if (max_omega < len_omega)
                    max_omega = sqrt(len_omega);
            }

            if (mag != 0)
            {
                double len_mag = vec_lenght(mag[i]);
                if (max_mag < len_mag)
                    max_mag = sqrt(len_mag);
            }
        }

        for (int i = 0; i < steps; i++)
        {
            #ifdef _WIN32
            strcpy_s(str_tmp, "t_%05lu.dump");
            sprintf_s(str_end, str_tmp, i);
            #elif __linux__
            strcpy(str_tmp, "t_%05lu.dump");
            sprintf(str_end, str_tmp, i);
            #endif

            string str_file = path_ovito + str_end;

            ofstream writer(str_file.c_str());

            if (writer.fail())
            {
                PRINT_ERROR("Cannot open file:\n      " + str_file);
                return false;
            }

            writer << "ITEM: TIMESTEP\n";
            writer << i << "\n";
            writer << "ITEM: NUMBER OF ATOMS\n";
            writer << Nmon << "\n";
            writer << "ITEM: BOX BOUNDS pp pp pp\n";
            
            #ifdef _WIN32
            strcpy_s(str_tmp, "%.4f %.4f\n");
            sprintf_s(str_end, str_tmp, -b_size, b_size);
            #elif __linux__
            strcpy(str_tmp, "%.4f %.4f\n");
            sprintf(str_end, str_tmp, -b_size, b_size);
            #endif

            writer << str_end;
            writer << str_end;
            writer << str_end;

            writer << "ITEM: ATOMS id mol type x y z vx vy vz fx fy fz tqx tqy tqz omegax omegay omegaz mux muy muz radius\n";

            for (int j = 0; j < Nmon; j++)
            {
                double x = 1.0e9 * pos[i * Nmon + j].x;
                double y = 1.0e9 * pos[i * Nmon + j].y;
                double z = 1.0e9 * pos[i * Nmon + j].z;

                int cl_id = 0;
                double vx = 0, vy = 0, vz = 0;
                double fx = 0, fy = 0, fz = 0;
                double tx = 0, ty = 0, tz = 0;

                double ox = 0, oy = 0, oz = 0;
                double mx = 0, my = 0, mz = 0;

                if (cluserIDs != 0)
                    cl_id = cluserIDs[i * Nmon + j];

                if (vel != 0)
                {
                    double len_vel = vec_lenght(vel[i * Nmon + j]);
                    if (len_vel > 0)
                    {
                        vx = vel[i * Nmon + j].x / len_vel;
                        vy = vel[i * Nmon + j].y / len_vel;
                        vz = vel[i * Nmon + j].z / len_vel;
                    }
                }

                if (force != 0)
                {
                    double len_force = vec_lenght(force[i * Nmon + j]);
                    if (len_force > 0)
                    {
                        fx = force[i * Nmon + j].x / len_force * pow(len_force, 1.0 / 8.0);
                        fy = force[i * Nmon + j].y / len_force * pow(len_force, 1.0 / 8.0);
                        fz = force[i * Nmon + j].z / len_force * pow(len_force, 1.0 / 8.0);
                    }
                }

                if (torque != 0)
                {
                    double len_torque = vec_lenght(torque[i * Nmon + j]);
                    if (len_torque > 0)
                    {
                        double3 torque_single = torque[i * Nmon + j];
                        vec_normalize(torque_single);
                        tx = torque_single.x;
                        ty = torque_single.y;
                        tz = torque_single.z;
                        //tx = torque[i * Nmon + j].x / len_torque;
                        //ty = torque[i * Nmon + j].y / len_torque;
                        //tz = torque[i * Nmon + j].z / len_torque;
                    }
                }

                if (omega != 0)
                {
                    double len_omega = vec_lenght(omega[i * Nmon + j]);
                    if (len_omega > 0)
                    {
                        //ox = omega[i * Nmon + j].x / len_omega;
                        //oy = omega[i * Nmon + j].y / len_omega;
                        //oz = omega[i * Nmon + j].z / len_omega;
                        double3 omega_single = omega[i * Nmon + j];
                        vec_normalize(omega_single);
                        ox = omega_single.x;
                        oy = omega_single.y;
                        oz = omega_single.z;
                    }
                }

                if (mag != 0)
                {
                    double len_mag = vec_lenght(mag[i * Nmon + j]);
                    if (len_mag != 0)
                    {
                        // Rescaling the magnetization to make it usefull.
                        mx = mag[i * Nmon + j].x / len_mag * pow(len_mag, 1.0 / 8.0);
                        my = mag[i * Nmon + j].y / len_mag * pow(len_mag, 1.0 / 8.0);
                        mz = mag[i * Nmon + j].z / len_mag * pow(len_mag, 1.0 / 8.0);
                    }
                }

                int mat_id = matID[j] + 1;

                double r = 1.0e9 * amon[j];
                
                #ifdef _WIN32
                strcpy_s(str_tmp, "%d %d %d %.5f %.5f %.5f %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5f\n");
                sprintf_s(str_end, str_tmp, j, cl_id, mat_id, x, y, z, vx, vy, vz, fx, fy, fz, tx, ty, tz, ox, oy, oz, mx, my, mz, r);
                #elif __linux__
                strcpy(str_tmp, "%d %d %d %.5f %.5f %.5f %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5f\n");
                sprintf(str_end, str_tmp, j, cl_id, mat_id, x, y, z, vx, vy, vz, fx, fy, fz, tx, ty, tz, ox, oy, oz, mx, my, mz, r);
                #endif
                
                writer << str_end;
            }
            writer.close();
        }

        return true;
    };

    /**
     * @brief Prints a summary of the run parameters.
     */
    void printParameters()
    {
        PRINT_TITLE("OVERVIEW OF RUN PARAMETERS");
        PRINT_CLR_LINE();
        
        std::cout << "  - command file :\n\t" << cmd_filename << "\n\n";

        std::cout << "  - Nr. of iterations: " << N_iter << "\n";
        std::cout << "  - start time   : " << time_start << " sec\n";
        std::cout << "  - stop time    : " << time_stop << " sec\n";
        std::cout << "  - time step    : " << time_step << " sec\n\n";

        std::cout << "  - position A: (" << pos_A.x << ", " << pos_A.y << ", " << pos_A.z << ") m\n";
        std::cout << "  - position B: (" << pos_B.x << ", " << pos_B.y << ", " << pos_B.z << ") m\n\n";

        std::cout << "  - velocity A: (" << vel_A.x << ", " << vel_A.y << " " << vel_A.z << ") m/s\n";
        std::cout << "  - velocity B: (" << vel_B.x << ", " << vel_B.y << " " << vel_B.z << ") m/s\n\n";

        std::cout << "  - ang. velocity A: (" << ang_A.x << ", " << ang_A.y << ", " << ang_A.z << ") m\n";
        std::cout << "  - ang. velocity B: (" << ang_B.x << ", " << ang_B.y << ", " << ang_B.z << ") m\n\n";

        if (Bext.x + Bext.y + Bext.z == 0)
            cout << "  - ext. mag. field: none\n";
        else
            cout << "  - ext. mag. field: (" << Bext.x << ", " << Bext.y << ", " << Bext.z << ") T\n";

        cout << "  - dust temp.     : " << T_dust << " K \n\n";

        if (save_pos || save_vel || save_force || save_torque || save_cluster)
        {
            cout << "  - Output files every " << N_save << "-th time step to: \n";

            if (save_pos)
                cout << "    - positions :\n\t" << path_binary << "pos.dat\n";

            if (save_vel)
                cout << "    - velocities:\n\t" << path_binary << "vel.dat\n";

            if (save_force)
                cout << "    - forces    :\n\t" << path_binary << "force.dat\n";

            if (save_torque)
                cout << "    - torques   :\n\t" << path_binary << "torque.dat\n";

            if (save_cluster)
                cout << "    - clusters  :\n\t" << path_binary << "cluster.dat\n";
        }

        PRINT_CLR_LINE();
        PRINT_TITLE("OVERVIEW OF MATERIAL PARAMETERS")
        PRINT_CLR_LINE();

        printf("  - surface energy   : [%6.2e , %6.2e] J m^-1\n", min_gamma, max_gamma);
        printf("  - Young's modulus  : [%6.2e , %6.2e] Pa\n", min_E, max_E);
        printf("  - Poisson's number : [%6.2e , %6.2e] \n", min_nu, max_nu);
        printf("  - density          : [%6.2e , %6.2e] kg m^-3\n", min_rho, max_rho);
        printf("  - crit. roll.      : [%6.2e , %6.2e] m\n", min_xi, max_xi);
        printf("  - visc. damp time  : [%6.2e , %6.2e] sec\n", min_Tvis, max_Tvis);

        cout << "\n";
        if (mat_type == MAT_TYPE_MAG)
        {
            if (min_tss != MIN_DEF && max_tss != MAX_DEF)
                cout << "  - spin-spin rel. time   : " << min_tss << " - " << max_tss << " sec    \n";

            if (min_tsl != MIN_DEF && max_tsl != MAX_DEF)
                cout << "  - spin-lattice rel. time: " << min_tsl << " - " << max_tsl << " sec    \n";

            if (min_Msat != MIN_DEF && max_Msat != MAX_DEF)
                cout << "  - sat. magnetization    : " << min_Msat << " - " << max_Msat << " A m^-1    \n";

            if (min_chi != MIN_DEF && max_chi != MAX_DEF)
                cout << "  - mag. suceptibility    : " << min_chi << " - " << max_chi << "     \n";

            if (min_Tc != MIN_DEF && max_Tc != MAX_DEF)
                cout << "  - Curie temperature     : " << min_Tc << " - " << max_Tc << " K    \n";

            cout << "\n";
        }

        PRINT_CLR_LINE();
        PRINT_TITLE("OVERVIEW OF AGGREGATE PARAMETERS");
        PRINT_CLR_LINE();
        
        cout << "  - Nr. of monomers: " << Nmon_A + Nmon_B << "\n";
        cout << "  - monomer radius : " << a_mon_min << " - " << a_mon_max << " m\n\n";

        cout << "  - file A:\n\t" << agg_filename_A << "\n";
        cout << "  - Nr. of monomers  A: " << Nmon_A << "\n";
        cout << "  - effective radius A: " << a_eff_A << " m\n";
        cout << "  - outer radius     A: " << a_out_A << " m\n\n";

        cout << "  - file B:\n\t" << agg_filename_B << "\n";
        cout << "  - Nr. of monomers  B: " << Nmon_B << "\n";
        cout << "  - effective radius B: " << a_eff_B << " m\n";
        cout << "  - outer radius     B: " << a_out_B << " m\n";

        PRINT_CLR_LINE();
    };

    // Defining getter functions.

    bool savePos() { return save_pos; }
    bool saveVel() { return save_vel; }
    bool saveForce() { return save_force; }
    bool saveTorque() { return save_torque; }
    bool saveCluster() { return save_cluster; }
    bool saveOvito() { return save_ovito; }

    bool saveOmega() { return save_omega; }
    bool saveMag() { return save_mag; }

    double3 getBext() { return Bext; }

    /**
     * @brief Writes a header file that contains important information for the run.
     * 
     * @return True on success, False on failure.
     */
    bool writeHeader()
    {
        string str_file = path_binary;
        str_file += "header.txt";

        int Nmat = int(lst_gamma.size());

        ofstream writer(str_file.c_str());

        if (writer.fail())
        {
            cout << "\nERROR: Cannot write header file:\n\t" << str_file << endl;
            return false;
        }

        writer << "#N_iter\tNmon_A\tNmon_B\tN_save\tN_mat\ttime_step [s]\n";
        writer << N_iter << "\t" << Nmon_A << "\t" << Nmon_B << "\t" << N_save << "\t" << Nmat << "\t" << time_step << "\n";

        writer.close();
        return true;
    };

    /**
     * @brief Saves a list of vectors into the binary output directory.
     * 
     * @param name_file: The name of the file.
     * @param *data: The array of vectors that is to be saved.
     * @param N: // TODO: Determine what exactly this is.
     * 
     * @return True on success, False on failure.
     */
    bool writeBinaryVec(string name_file, const double3* data, ullong N)
    {
        string path_tmp = path_binary + name_file;
        ullong len_array = N * sizeof(double3);

        ofstream bin_writer(path_tmp.c_str(), ios::binary);

        if (bin_writer.fail())
        {
            cout << "\nERROR: Cannot write binary file:\n\t";
            cout << path_tmp << "\n";
            return false;
        }

        PRINT_LOG(std::string("Writing binary file to:\n     ") + path_tmp, 2);

        bin_writer.write((const char*)data, len_array);

        return true;
    };

    /**
     * @brief Saves a list of doubles into the ouput directory.
     * 
     * @param name_file: The name of the file.
     * @param *data: An array of doubles that is to be saved.
     * @param N: // TODO: Determine what exactly this is.
     */
    bool writeBinaryDouble(string name_file, const double* data, ullong N)
    {
        string path_tmp = path_binary + name_file;
        ullong len_array = N * sizeof(double);

        ofstream bin_writer(path_tmp.c_str(), ios::out | ios::binary);

        if (bin_writer.fail())
        {
            PRINT_ERROR("Cannot write binary file:\n      " + path_tmp);
            return false;
        }

        PRINT_LOG(std::string("Writing binary file to:\n     ") + path_tmp, 2);

        bin_writer.write((const char*)data, len_array);

        return true;
    };

    /**
     * @brief Saves a list of integers into the ouput directory.
     * 
     * @param name_file: The name of the file.
     * @param *data: An array of doubles that is to be saved.
     * @param N: // TODO: Determine what exactly this is.
     */
    bool writeBinaryInt(string name_file, const int* data, ullong N)
    {
        string path_tmp = path_binary + name_file;
        ullong len_array = N * sizeof(int);

        ofstream bin_writer(path_tmp.c_str(), ios::out | ios::binary);

        if (bin_writer.fail())
        {
            cout << "\nERROR: Cannot write binary file:\n\t";
            cout << path_tmp << "\n";
            return false;
        }

        PRINT_LOG(std::string("Writing binary file to:\n     ") + path_tmp, 2);

        bin_writer.write((const char*)data, len_array);

        return true;
    };

    //getter and setter
    ullong getNIter() { return N_iter; }
    ullong getNSave() { return N_save; }

    double getTimeStep() { return time_step; };

private:
    /**
     * @brief Creates a directory.
     * 
     * @param path: The path of the directory.
     * @return True if success, False if failure.
     */
    bool createPath(string path) 
    {
        #ifdef _WIN32
        // Windows-specific 
        if (_mkdir(path.c_str()) == 0)
        {
            cout << "Directory created successfully on Windows:\n\t" << path << std::endl;
            return true;
        }
        else
        {
            if (errno == EEXIST)
            {
                cout << "Directory exists:\n\t" << path << std::endl;
                return true;
            }
            else
            {
                if (errno == ENOENT)
                {
                    cout << "ERROR: Path does not exist:\n\t" << path << std::endl;
                    return false;

                    /*stringstream str_strem;
                    string segment;
                    vector<std::string> seglist;

                    str_strem << path;

                    while (getline(str_strem, segment, SEP))
                    {
                        seglist.push_back(segment);
                    }

                    //for(int i=0;i< seglist.size(), i++)*/


                }
                else
                {
                    cout << "Failed to create directory on Windows! " << errno << std::endl;
                    cout << "Directory:\n\t" << path << std::endl;

                    return false;
                }
            }
        }
        #elif __linux__
        // This function creates the directory with full permissions for everyone. // TODO: Restrict permissions?
        if (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) == 0)
        {
            PRINT_LOG("Creating directory:\n      " + path, 3);
            return true;
        } 
        else 
        {
            if (errno == EEXIST)
            {
                PRINT_LOG("Directory allready exists:\n      " + path, 2);
                return true;
            }
            else
            {
                std::cout << "Could not create directory." << std::endl << std::flush;
                //PRINT_ERROR("Failed to create directory: " + strerror(errno));
            }
        }
        #else
        PRINT_ERROR("Unsupported operating system.")
        #endif
        
        return false;
    };

    // Field definition starts here:
    string cmd_filename; // The path to the command file.
    string path_results; // The path to the output directory.

    string path_binary; // The path to the binary output directory.
    string path_ovito; // The path to the ovito output directory.
    string path_plots; // The path to the plot directory.

    string path_A; // Temporary variable for the path to aggregate A.
    string path_B; // Temporary variable for the path to aggregate B.

    string agg_filename_A; // Path to aggregate A.
    string agg_filename_B; // Path to aggregate B.

    vec3D pos_A; // The position of aggregate A.
    vec3D pos_B; // The position of aggregate B.

    vec3D vel_A; // The initial velocity of aggregate A.
    vec3D vel_B; // The initial velocity of aggregate B.

    vec3D ang_A; // The initial angular velocity of aggregate A.
    vec3D ang_B; // The initial angular velocity of aggregate B.

    ullong N_iter; // The number of total iterations.
    ullong N_save; // The save interval.

    double time_start;
    double time_stop;
    double time_step; // The timestep size of the simulation.

    strlist lst_matName; // List of material names.
    ilist lst_matID;  // List of material ids.
    dlist lst_gamma;  // List of material surface energies
    dlist lst_E;      // List of material Young's moduli
    dlist lst_nu;     // List of material Poisson numbers
    dlist lst_rho;    // List of material densities.
    dlist lst_xi;     // List of material critical rolling lengths
    dlist lst_Tvis;   // List of materialviscous dumping time 

    dlist lst_tss;     // List of material spin-spin time scale
    dlist lst_tsl;     // List of material spin-lattice time scale
    dlist lst_Msat;    // List of material spont. magnetization
    dlist lst_chi;     // List of material mag. susceptibillity
    dlist lst_Tc;      // List of material Curie temperature

    // TODO: Check if this is used to control the behavior of the code or if its just for printing.
    int mat_type; // MAT_TYPE_NONE if no material has magnetic properties defined. MAT_TYPE_MAG if at least on material has magnetic properties defined.

    double a_eff_A; // Effective radius of aggregate A.
    double a_eff_B; // Effective radius of aggregate B.

    double a_out_A; // Outer radius of aggregate A (???).
    double a_out_B; // Outer radius of aggregate B (???)

    double a_mon_min; // ???
    double a_mon_max; // ???

    int Nmon_A; // Number of monomers in aggregate A.
    int Nmon_B; // Number of monomers in aggregate B.

    bool save_pos;
    bool save_vel;
    bool save_force;
    bool save_torque;
    bool save_cluster;
    bool save_omega;
    bool save_mag;
    bool save_ovito;

    double3 Bext; // The external magnetic field.
    double T_dust; // The temperature of the dust.

    // These are only here for pretty printing...
    double min_gamma;
    double min_E;
    double min_nu;
    double min_rho;
    double min_xi;
    double min_Tvis;

    double min_tss;
    double min_tsl;
    double min_Msat;
    double min_chi;
    double min_Tc;

    double max_gamma;
    double max_E;
    double max_nu;
    double max_rho;
    double max_xi;
    double max_Tvis;

    double max_tss;
    double max_tsl;
    double max_Msat;
    double max_chi;
    double max_Tc;
};
