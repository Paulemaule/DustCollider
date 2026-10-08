#pragma once

#include "commandFile.cuh"
#include "aggregate.cuh"

#include "../utils/logging.cuh"
#include "../utils/errors.cuh"
#include "../utils/vector.cuh"
#include "../physics/integrator_utils.cuh"
#include "../utils/constant.cuh"

#include <climits>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <string>
/**
 * TODO: Replace Macro based logging with dedicated Logging tool
 */
#include "../utils/printing.cuh"

/**
 * TODO: Documentation
 */
class Pipeline {
public:
    Pipeline() {
    }

    ~Pipeline() {
    }

    const SimulationConfig& getConfig() const {
        return run_config; 
    }

    /**
     * @brief Runs the setup for the simulation.
     * 
     * This function will 
     *  1) Parse the command line input.
     *  2) Read the command file and extract the run parameters.
     *  3) Read the aggregate files to extract monomer information.
     *  4) Build the initial system state.
     *  5) Calculate system parameters (timestep)
     *  6) Resolve the run length and snapshot schedule into iteration counts.
     *  7) Run a sanity check on the simulation config.
     * 
     * @param argc
     * @param argv
     */
    Status run(int argc, const char** argv) {
        // Check and log the version and build of the code.
        std::string title = std::string("DUST COLLIDER - ") + VERSION;

        Logger::title(title);
        Logger::lineBreak();

        #if defined(RELEASE)
            Logger::log("Code compiled in Release build");
        #elif defined(DEBUG)
            Logger::warn("Code compiled in Debug build");
        #elif defined(TEST)
            Logger::warn("Code compiled in Test build");
        #else
            Logger::error("Compiler did not pass build information.");
            return Status::error;
        #endif

        Logger::log("Run started at {}.", current_time_string());

        Logger::lineBreak();

        Logger::header("SETUP");
        Logger::lineBreak();

        // Parse the command line arguments
        std::string command_file_path;
        Status _s = parse_commandline(argc, argv, command_file_path);
        if ( _s != Status::ok ) return _s;

        // Parse the command file into the member
        _s = parse_commandfile(command_file_path, run_config);
        if ( _s != Status::ok ) return _s;

        // Load aggregate files
        std::vector<Aggregate> aggregates;
        for ( const AggregateConfig& agg_cfg : run_config.aggregates ) {
            aggregates.push_back(Aggregate::from_file(agg_cfg.path));
        }

        // Calculate the initial system state from the aggregates
        _s = build_initial_state(aggregates);
        if ( _s != Status::ok ) return _s;

        // Calculate system parameters
        _s = calculate_system_properties();
        if ( _s != Status::ok ) return _s;

        // Convert the run length and snapshot schedule into iterations. Needs the timestep.
        _s = resolve_run_schedule();
        if ( _s != Status::ok ) return _s;

        Logger::lineBreak();
        
        // Pretty print the entire simulation config
        print_config(run_config);

        // Sanity checks on the simulation config
        Logger::separator();
        Logger::lineBreak();

        _s = check_simulation_config(run_config);
        if ( _s != Status::ok ) return _s;

        return Status::ok;
    }

private:
    SimulationConfig run_config;

    /**
     * @brief Populates run_config.initial_state from loaded aggregate data.
     */
    Status build_initial_state(const std::vector<Aggregate>& aggregates) {
        Logger::log("Prepare initial state from aggregate data.");

        size_t Nmon = 0;
        for ( const Aggregate& agg : aggregates ) Nmon += agg.header.Nmon;

        run_config.initial_state.positions.resize(Nmon);
        run_config.initial_state.velocities.resize(Nmon);
        run_config.initial_state.omegas.resize(Nmon);
        run_config.initial_state.radii.resize(Nmon);
        run_config.initial_state.masses.resize(Nmon);
        run_config.initial_state.moments.resize(Nmon);
        run_config.initial_state.material_ids.resize(Nmon);

        size_t i = 0;
        for ( size_t j = 0; j < aggregates.size(); j++ ) {
            const Aggregate&       agg     = aggregates[j];
            const AggregateConfig& agg_cfg = run_config.aggregates[j];

            for ( size_t k = 0; k < (size_t)agg.header.Nmon; k++ ) {
                double3 position;
                position.x = agg.monomers.positions[k].x + agg_cfg.position.x;
                position.y = agg.monomers.positions[k].y + agg_cfg.position.y;
                position.z = agg.monomers.positions[k].z + agg_cfg.position.z;

                double3 vel_tang = vec_cross(agg_cfg.angular, agg.monomers.positions[k]);
                double3 velocity;
                velocity.x = agg_cfg.velocity.x + vel_tang.x;
                velocity.y = agg_cfg.velocity.y + vel_tang.y;
                velocity.z = agg_cfg.velocity.z + vel_tang.z;

                double r   = agg.monomers.radii[k];
                int    mid = agg.monomers.material_ids[k];
                double rho = run_config.materials[mid].rho;
                double m   = (4.0 / 3.0) * PI * rho * r * r * r;

                run_config.initial_state.positions[i]    = position;
                run_config.initial_state.velocities[i]   = velocity;
                run_config.initial_state.omegas[i]       = agg_cfg.angular;
                run_config.initial_state.radii[i]        = r;
                run_config.initial_state.masses[i]       = m;
                run_config.initial_state.moments[i]      = (2.0 / 5.0) * m * r * r;
                run_config.initial_state.material_ids[i] = mid;

                i++;
            }
        }

        return Status::ok;
    }

    /**
     * @brief Auto-calculates several system properties.
     *
     * This includes the smallest dynamical timescale tau_min and, unless <time_step> was set, the timestep.
     * tau_min is calculated in both cases, as the run schedule is reported in units of it.
     */
    Status calculate_system_properties() {
        Logger::log("Calculate simulation timestep.");

        const bool skip_timestep = (run_config.timestep != 0.0);
        if ( skip_timestep ) {
            Logger::warn("<time_step> was set in command file. Timestep was set to {} s.", run_config.timestep);
        }

        const size_t Nmon = run_config.initial_state.radii.size();
        double tau_min  = 1e200; // The smallest dynamical timescale of the system.
        double min_Fc_R = 1e200; // The smallest contact torque scale F_c * R_red [N*m].

        // Iterate over all monomer pairs to determine the systems extremal dynamical parameters
        for ( size_t i = 0; i < Nmon; i++ ) {
            for ( size_t j = i + 1; j < Nmon; j++ ) {
                // Find the materials of the two monomers
                const MaterialEntry& mat_i = run_config.materials[run_config.initial_state.material_ids[i]];
                const MaterialEntry& mat_j = run_config.materials[run_config.initial_state.material_ids[j]];

                // Read the properties of the monomers
                double m_i   = run_config.initial_state.masses[i];
                double m_j   = run_config.initial_state.masses[j];

                double r_i   = run_config.initial_state.radii[i];
                double r_j   = run_config.initial_state.radii[j];

                // Calculate pair properties. These MUST come from integrator_utils.cuh rather than be
                // re-derived here, or the timestep silently drifts away from the kernels that use them.
                double M            = (m_i * m_j) / (m_i + m_j);
                double R            = get_R(r_i, r_j);

                double Es           = get_E_s(mat_i.E, mat_j.E, mat_i.nu, mat_j.nu);
                double gamma_ij     = get_gamma(mat_i.gamma, mat_j.gamma);

                double a0           = get_a_0(gamma_ij, R, Es);
                double delta_N_c    = get_delta_N_crit(a0, R);
                double F_c          = 3.0 * PI * gamma_ij * R;

                // The timescale(s) of the system
                double tau_N        = sqrt(M * delta_N_c / F_c);

                // TODO: Calculating the timestep only from the normal interaction only works if the normal direction actually dominates the other interactions. This is not necessarily given for all material parameters.
                if ( tau_N < tau_min  ) tau_min  = tau_N;
                if ( F_c * R < min_Fc_R ) min_Fc_R = F_c * R;
            }
        }

        // Without a monomer pair there is no contact timescale.
        if ( Nmon >= 2 ) {
            run_config.tau_min = tau_min;
        } else {
            Logger::warn("The timestep could not be calculated for a simulation with less then 2 monomers!");
        }

        if ( !skip_timestep ) {
            run_config.timestep = 0.005 * tau_min;
        }

        return Status::ok;
    }

    // TODO There is a structural conflict baked into SimulationConfig that becomes clearest here. SimulationConfig is meant to be a struct that is handed to the Simulation and containing only info relevant to the Simulation. Currently it also holds all the parameters set in the command file. It would make more sense to split the current struct into two structs one for each purpose (SimulationParameters and Input).
    /**
     * @brief Resolves the run length and the snapshot schedule into iteration counts.
     *
     * The run length is set by exactly one of <N_iter> or <t_end> [s].
     * The snapshot schedule is set by exactly one of <N_save>, <t_save> [s] or <N_snap>. 
     * This function converts any pairing of these tags into a valid run lenght and 
     * snapshot schedule while taking the timestep of the simulation into account.
     *
     * Snapshots lie on a uniform grid of N_save iterations, starting with the initial state at t = 0.
     * N_iter is rounded up to a whole number of snapshot intervals, so the final state is always stored.
     * On return N_iter and N_save hold the resolved values, while t_end, t_save and N_snap keep the
     * values from the command file. // TODO: See the comment on top of this funtion for why this is a bad idea.
     */
    Status resolve_run_schedule() {
        Logger::log("Resolve run length and snapshot schedule.");

        SimulationConfig& cfg = run_config;
        const double      dt  = cfg.timestep;

        // Exactly one tag of each group must be set. Like <time_step>, a tag counts as set if it is non-zero.
        const int n_length = (cfg.N_iter != 0) + (cfg.t_end != 0.0);
        if ( n_length == 0 ) {
            Logger::error("The run length is not set. Use either <N_iter> or <t_end>.");
            return Status::error;
        }
        if ( n_length > 1 ) {
            Logger::error("<N_iter> and <t_end> are mutually exclusive, only one may be set.");
            return Status::error;
        }

        const int n_snapshot = (cfg.output.N_save != 0) + (cfg.output.t_save != 0.0) + (cfg.output.N_snap != 0);
        if ( n_snapshot == 0 ) {
            Logger::error("The snapshot schedule is not set. Use one of <N_save>, <t_save> or <N_snap>.");
            return Status::error;
        }
        if ( n_snapshot > 1 ) {
            Logger::error("<N_save>, <t_save> and <N_snap> are mutually exclusive, but {} of them are set.", n_snapshot);
            return Status::error;
        }

        // The minimum number of iterations the run has to cover. Extracted from <N_iter> of <t_end>, whichever is set.
        int64_t N_min = 0;

        if ( cfg.N_iter != 0 ) { // When the number of iterations is explicitly set
            if ( cfg.N_iter < 1 ) {
                Logger::error("<N_iter> = {} must be > 0.", cfg.N_iter);
                return Status::error;
            }
            N_min = cfg.N_iter;
        } else { // When the number of iterations is via the lenght of the simulation
            if ( !is_finite(cfg.t_end) || cfg.t_end <= 0.0 ) {
                Logger::error("<t_end> = {} s must be > 0.", cfg.t_end);
                return Status::error;
            }
            const double x = ceil_tol(cfg.t_end / dt);
            if ( !(x >= 1.0 && x <= double(INT_MAX)) ) {
                Logger::error("<t_end> = {:.4e} s corresponds to {:.4e} iterations of dt = {:.4e} s, but must be between 1 and {} iterations.",
                              cfg.t_end, x, dt, INT_MAX);
                return Status::error;
            }
            N_min = int64_t(x);
        }

        // The number of iterations between two snapshots and the number of snapshot intervals.
        int64_t G = 0;
        int64_t I = 0;

        if ( cfg.output.N_save != 0 ) { // If the number of iterations between snapshots is set
            if ( cfg.output.N_save < 1 ) {
                Logger::error("<N_save> = {} must be > 0.", cfg.output.N_save);
                return Status::error;
            }
            G = cfg.output.N_save;
            I = ceil_div(N_min, G);
        } else if ( cfg.output.t_save != 0.0 ) { // If the time interval between snapshots is set
            if ( !is_finite(cfg.output.t_save) || cfg.output.t_save <= 0.0 ) {
                Logger::error("<t_save> = {} s must be > 0.", cfg.output.t_save);
                return Status::error;
            }
            // Calculate the number of interations based on the timestep and time interval between snapshots
            const double ratio = cfg.output.t_save / dt;
            const double x     = ceil_tol(ratio); // Round up with tolerance
            if ( !(ratio * (1.0 + 1e-9) >= 1.0 && x <= double(INT_MAX)) ) {
                Logger::error("<t_save> = {:.4e} s corresponds to {:.4e} iterations of dt = {:.4e} s, but must be between 1 and {} iterations.",
                              cfg.output.t_save, ratio, dt, INT_MAX);
                return Status::error;
            }
            G = int64_t(x);
            I = ceil_div(N_min, G);
        } else { // If the number of snapshots is set
            if ( cfg.output.N_snap < 2 ) {
                Logger::error("<N_snap> = {} must be >= 2, as the initial and the final state are always stored.", cfg.output.N_snap);
                return Status::error;
            }
            I = cfg.output.N_snap - 1;
            G = ceil_div(N_min, I);
        }
        
        const int64_t N_iter = G * I;
        if ( N_iter > INT_MAX ) {
            Logger::error("The run needs {} iterations ({} snapshot intervals of {} iterations), more than the maximum of {}.",
                          N_iter, I, G, INT_MAX);
            return Status::error;
        }

        // Extend the runtime to allow for the last snapshot to be stored
        if ( N_iter > N_min ) {
            Logger::warn("The run was extended from {} to {} iterations (+{:.4e} s) to end on a snapshot.",
                         N_min, N_iter, double(N_iter - N_min) * dt);
        }

        cfg.N_iter        = int(N_iter);
        cfg.output.N_save = int(G);

        return Status::ok;
    }

    /**
     * @brief Parse the command line.
     */
    Status parse_commandline(int argc, const char** argv, std::string& command_file_path) {
        Logger::log("Parsing command line input.");

        if (argc != 2) {
            Logger::error("Wrong number ({} !=) of command line inputs. Command line should be 'dust_collider <command_file_path>'", argc);
            return Status::error;
        }

        command_file_path = argv[1];
        return Status::ok;
    }

    /**
     * @brief Parses the command file and return the filled out simulation config.
     */
    Status parse_commandfile(const std::string& command_file_path, SimulationConfig& out_config) {
        CommandFile command_file(command_file_path);
        Status _s = command_file.parse(out_config);
        if ( _s != Status::ok ) return _s;

        // Resolve relative paths against the command file's directory
        std::filesystem::path command_file_parent = std::filesystem::path(command_file_path).parent_path();

        // Helper function that resolves file system paths inplace
        auto resolve = [&](std::string& p) {
            if (!p.empty() && !std::filesystem::path(p).is_absolute())
                p = std::filesystem::weakly_canonical(command_file_parent / p).string();
        };

        // Resolve out path
        resolve(out_config.output.path);

        // Resolve aggregate paths
        for (auto& agg : out_config.aggregates)
            resolve(agg.path);

        return Status::ok;
    }

    static inline bool is_finite(double x) {
        return std::isfinite(x);
    }

    static inline bool is_finite3(const double3& v) {
        return is_finite(v.x) && is_finite(v.y) && is_finite(v.z);
    }

    // Integer division rounding up, for a >= 0 and b > 0.
    static inline int64_t ceil_div(int64_t a, int64_t b) {
        return (a + b - 1) / b;
    }

    // Rounds a time in units of dt up to whole iterations with tolerance to avoid issues from flop noise.
    static inline double ceil_tol(double x) {
        return std::ceil(x * (1.0 - 1e-9));
    }

    static inline bool needs_output(const SimulationConfig& c) {
        return c.output.ovito   || c.output.position || c.output.velocity ||
               c.output.angular || c.output.force    || c.output.torque;
    }

    static inline Status check_aggregate_config(const AggregateConfig& a) {
        if (a.path.empty()) {
            Logger::error("Aggregate '{}'s path '{}' is empty.", a.name, a.path);
            return Status::error;
        }
        if (!is_finite3(a.position)) {
            Logger::error("Aggregate '{}'s position '({:.2e}, {:.2e}, {:.2e})' contains non-finite values.", 
                            a.name, a.position.x, a.position.y, a.position.z);
            return Status::error;
        }
        if (!is_finite3(a.velocity)) {
            Logger::error("Aggregate '{}'s velocity '({:.2e}, {:.2e}, {:.2e})' contains non-finite values.", 
                            a.name, a.velocity.x, a.velocity.y, a.velocity.z);
            return Status::error;
        }
        if (!is_finite3(a.angular)) {
            Logger::error("Aggregate '{}'s angular velocity '({:.2e}, {:.2e}, {:.2e})' contains non-finite values.",
                            a.name, a.angular.x, a.angular.y, a.angular.z);
            return Status::error;
        }
        return Status::ok;
    }

    /**
     * @brief Check the simulation config for errors or unusual setups.
     */
    Status check_simulation_config(const SimulationConfig& cfg) {
        Logger::log("Validating simulation config.");

        // N_iter and N_save are validated in resolve_run_schedule(), which runs before the config is printed.

        if (cfg.aggregates.empty()) {
            Logger::error("No aggregates defined in command file.");
            return Status::error;
        }

        for ( const AggregateConfig& a : cfg.aggregates ) {
            if ( check_aggregate_config(a) != Status::ok ) return Status::error;
        }

        if ( needs_output(cfg) && cfg.output.path.empty() ) {
            Logger::error("Output requested (save_* flags set) but <path_results> is empty.");
            return Status::error;
        }

        // Material checks
        if (cfg.materials.empty()) {
            Logger::error("No materials defined.");
            return Status::error;
        }

        for (int i = 0; i < (int)cfg.materials.size(); i++) {
            const MaterialEntry& mat = cfg.materials[i];

            // Msat and Tc must both be zero (non-magnetic) or both non-zero (magnetic).
            if ((mat.Msat == 0.0) != (mat.Tc == 0.0)) {
                Logger::error("Material '{}': Msat and Tc must both be zero or both be non-zero.", mat.name);
                return Status::error;
            }
        }

        // Warn if B_ext is set but no material is magnetic, or vice versa.
        const bool has_bext = vec_length_sq(cfg.B_ext) > 0.0;
        bool has_mag_mat = false;
        for (const MaterialEntry& mat : cfg.materials) {
            if (mat.Msat != 0.0 || mat.chi != 0.0) { has_mag_mat = true; break; }
        }
        if (has_bext && !has_mag_mat)
            Logger::warn("<B_ext> is set but no material has magnetic properties.");
        if (!has_bext && has_mag_mat)
            Logger::warn("A magnetic material is defined but <B_ext> is zero.");

        // Dust temperature
        if (cfg.T_dust == -1.0) {
            Logger::warn("T_dust = -1: temperature corrections are disabled.");
        } else if (cfg.T_dust <= 0.0) {
            Logger::error("T_dust must be > 0 (or -1 to disable temperature corrections).");
            return Status::error;
        }

        // Warn if all aggregate positions are at the origin
        bool any_nonzero_pos = false;
        for (const AggregateConfig& a : cfg.aggregates) {
            if (vec_length_sq(a.position) > 0.0) { any_nonzero_pos = true; break; }
        }
        if (!any_nonzero_pos) {
            Logger::warn("All aggregates are centred at the origin.");
        }

        // Warn if all velocities are zero — simulation will be static
        bool any_motion = false;
        for ( const AggregateConfig& a : cfg.aggregates ) {
            if ( vec_length_sq(a.angular) + vec_length_sq(a.velocity) > 0.0 )
                any_motion = true;
        }
        if ( !any_motion ) {
            Logger::warn("The simulation is static — all aggregate velocities are zero.");
        }

        return Status::ok;
    }

    /**
     * @brief Pretty prints the simulation configuration to console.
     * 
     * @param cfg The configuration to print.
     */
    static void print_config(const SimulationConfig& cfg) {
        // Run parameters
        Logger::header("SIMULATION CONFIG");
        Logger::lineBreak();

        // Expresses a time in units of the smallest contact timescale, if it is available.
        auto in_tau = [&](double t) -> std::string {
            return (cfg.tau_min > 0.0) ? std::format(" = {:.4g} tau_{{dyn,min}}", t / cfg.tau_min) : std::string();
        };

        const double t_run = double(cfg.N_iter) * cfg.timestep;

        Logger::print("   Nmon:      {}", cfg.initial_state.positions.size());
        Logger::print("   timestep:  {:.4e} s{}", cfg.timestep, in_tau(cfg.timestep));
        if (cfg.tau_min > 0.0) {
            Logger::print("   tau_{{N,min}}: {:.4e} s", cfg.tau_min);
        } else {
            Logger::print("   tau_{{N,min}}: n/a (fewer than two monomers)");
        }
        if (cfg.t_end != 0.0) {
            Logger::print("   N_iter:    {} (from <t_end> = {:.4e} s)", cfg.N_iter, cfg.t_end);
        } else {
            Logger::print("   N_iter:    {} (from <N_iter>)", cfg.N_iter);
        }
        Logger::print("      =>      simulated time {:.4e} s{}", t_run, in_tau(t_run));
        Logger::print("   T_dust:    {:.2g} K", cfg.T_dust);
        if (cfg.T_dust == -1.0) {
            Logger::print("      T_dust:   disabled");
        }
        Logger::print("   B_ext:     ({:.3e}, {:.3e}, {:.3e}) T",
               cfg.B_ext.x, cfg.B_ext.y, cfg.B_ext.z);

        Logger::lineBreak();

        // Aggregates
        Logger::print("Aggregate:");
        Logger::lineBreak();
        for (const AggregateConfig& a : cfg.aggregates) {
            Logger::print("   [{}] from {}", a.name.c_str(), a.path.c_str());
            Logger::print("      pos: ({:.3e}, {:.3e}, {:.3e}) m",
                   a.position.x, a.position.y, a.position.z);
            Logger::print("      vel: ({:.3e}, {:.3e}, {:.3e}) m/s",
                   a.velocity.x, a.velocity.y, a.velocity.z);
            Logger::print("      ang: ({:.3e}, {:.3e}, {:.3e}) rad/s",
                   a.angular.x, a.angular.y, a.angular.z);
        }

        Logger::lineBreak();

        // Materials
        Logger::print("Materials:");
        Logger::lineBreak();
        Logger::print("   - Mechanical properties -");
        Logger::print("   {:<2} | {:<20} | {:<9} | {:<9} | {:<9} | {:<9} | {:<9} | {:<9}",
            "ID", "name", "gamma", "E", "nu", "rho", "xi", "tvis");
        Logger::print("   {:<2} | {:<20} | {:<9} | {:<9} | {:<9} | {:<9} | {:<9} | {:<9}",
            "", "", "J/m^2", "Pa", "", "kg/m^3", "m", "s");

        Logger::lineBreak();

        bool magnetic_print = false;
        for ( int i = 0; i < (int)cfg.materials.size(); i++ ) {
            const MaterialEntry& m = cfg.materials[i];
            Logger::print("   {:<2} | {:>20} | {:>9.2e} | {:>9.2e} | {:>9.2e} | {:>9.2e} | {:>9.2e} | {:>9.2e}", 
                i, m.name, m.gamma, m.E, m.nu, m.rho, m.xi, m.tvis);

            // Check if any of the materials have magnetic properties defined
            if ( m.Msat != 0.0 || m.chi != 0.0 ) {
                magnetic_print = true;
            }
        }

        Logger::lineBreak();

        if (magnetic_print) {
            Logger::print("   - Magnetic properties -");
            Logger::print("   {:<2} | {:<20} | {:<9} | {:<9} | {:<9} | {:<9} | {:<9}",
                "ID", "name", "tss", "tsl", "Msat", "chi", "Tc");
            Logger::print("   {:<2} | {:<20} | {:<9} | {:<9} | {:<9} | {:<9} | {:<9}",
                "", "", "s", "s", "", "A/m", "", "K");
            Logger::lineBreak();

            for ( int i = 0; i < (int)cfg.materials.size(); i++ ) {
                const MaterialEntry& m = cfg.materials[i];
                if ( m.Msat != 0.0 || m.chi != 0.0 ) {
                    Logger::print("   {:<2} | {:<20} | {:<9.2e} | {:<9.2e} | {:<9.2e} | {:<9.2e} | {:<9.2e}",
                        i, m.name, m.tss, m.tsl, m.Msat, m.chi, m.Tc);
                } else {
                    Logger::print("   {:<2} | {:<20} | XXX",
                        i, m.name);
                }
            }

            Logger::lineBreak();
        }

        // Output
        Logger::print("Simulation Output");
        Logger::lineBreak();
        const double t_save = double(cfg.output.N_save) * cfg.timestep;

        Logger::print("   Path:      {}", cfg.output.path.c_str());
        if (cfg.output.N_snap != 0) {
            Logger::print("   N_save:    {} (from <N_snap> = {})", cfg.output.N_save, cfg.output.N_snap);
        } else if (cfg.output.t_save != 0.0) {
            Logger::print("   N_save:    {} (from <t_save> = {:.4e} s)", cfg.output.N_save, cfg.output.t_save);
        } else {
            Logger::print("   N_save:    {} (from <N_save>)", cfg.output.N_save);
        }
        Logger::print("      =>      a snapshot every {:.4e} s{}", t_save, in_tau(t_save));
        Logger::print("      =>      {} snapshots from t = 0 to {:.4e} s will be saved",
            cfg.N_iter / cfg.output.N_save + 1, t_run);
        Logger::print("   ovito={:<3}  pos={:<3}  vel={:<3}  ang={:<3}  force={:<3}  torque={:<3}",
            cfg.output.ovito    ? "yes" : "no",
            cfg.output.position ? "yes" : "no",
            cfg.output.velocity ? "yes" : "no",
            cfg.output.angular  ? "yes" : "no",
            cfg.output.force    ? "yes" : "no",
            cfg.output.torque   ? "yes" : "no");

        Logger::lineBreak();
    }
};
