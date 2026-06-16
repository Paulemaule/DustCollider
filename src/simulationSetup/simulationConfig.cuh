#pragma once

#include <vector>
#include <string>

#include "aggregate.cuh"

/**
 * @brief Initial conditions for a single aggregate in the simulation.
 *
 * @param name      Key extracted from the command file tag, e.g. "A", "spin".
 * @param path      Path to the aggregate file.
 * @param position  Initial position offset [m].
 * @param velocity  Initial translational velocity [m/s].
 * @param angular   Initial angular velocity [rad/s].
 */
struct AggregateConfig {
    std::string  name;
    std::string  path{};
    double3      position{};
    double3      velocity{};
    double3      angular{};
};

/**
 * @brief Parameters for a single material as defined in the command file.
 *
 * Non-magnetic materials use only the contact mechanics fields (gamma through tvis).
 * Magnetic materials additionally define tss, tsl, Msat, chi and Tc; the remaining
 * fields are left at their default value of 0.
 *
 * @param name   Display name of the material (e.g. "forsterite").
 * @param gamma  Surface energy [J/m²].
 * @param E      Young's modulus [Pa].
 * @param nu     Poisson's ratio [-].
 * @param rho    Density [kg/m³].
 * @param xi     Critical rolling displacement [m].
 * @param tvis   Viscous damping timescale [s].
 * @param tss    Spin-spin relaxation time [s].
 * @param tsl    Spin-lattice relaxation time [s].
 * @param Msat   Saturation magnetization [A/m].
 * @param chi    Magnetic susceptibility [-].
 * @param Tc     Curie temperature [K].
 */
struct MaterialEntry {
    std::string  name;
    double  gamma = 0.0;
    double  E     = 0.0;
    double  nu    = 0.0;
    double  rho   = 0.0;
    double  xi    = 0.0;
    double  tvis  = 0.0;
    double  tss   = 0.0;
    double  tsl   = 0.0;
    double  Msat  = 0.0;
    double  chi   = 0.0;
    double  Tc    = 0.0;
};

/**
 * @brief Output configuration for the simulation.
 *
 * @param path      Path to the output directory.
 * @param N_save    Interval (in iterations) at which data is written to disk.
 * @param ovito     Write trajectory output in OVITO-compatible format.
 * @param position  Store monomer positions to disk.
 * @param velocity  Store monomer velocities to disk.
 * @param angular   Store monomer angular velocities to disk.
 * @param force     Store monomer forces to disk.
 * @param torque    Store monomer torques to disk.
 * @param energy    Store potential and dissipated energy to disk.
 */
struct OutputConfig {
    std::string  path{};
    int          N_save   = 0;
    bool         ovito    = false;
    bool         position = false;
    bool         velocity = false;
    bool         angular  = false;
    bool         force    = false;
    bool         torque   = false;
    bool         energy   = false;
};

/**
 * @brief Initial state of the aggregate model calculated directly from the aggregate files and command file inputs.
 *
 * Uses plain std::vector storage (no pinned or device memory).
 *
 * @param positions     Initial position of each monomer [m].
 * @param velocities    Initial velocity of each monomer [m/s].
 * @param omegas        Initial angular velocity of each monomer [rad/s].
 * @param radii         Radius of each monomer [m].
 * @param material_ids  0-indexed material ID of each monomer.
 */
struct InitialState {
    std::vector<double3>  positions;
    std::vector<double3>  velocities;
    std::vector<double3>  omegas;
    std::vector<double>   radii;
    std::vector<double>   masses;
    std::vector<double>   moments;
    std::vector<int>      material_ids;
};

/**
 * @brief Complete simulation setup as defined in the command file.
 *
 * The fields of this struct describe the setup of the simulation.
 *
 * @param aggregates      Initial conditions for each aggregate, keyed by name.
 * @param materials       Material parameters indexed by (command-file id - 1).
 * @param output          Output path and save flags.
 * @param N_iter          Total number of simulation iterations.
 * @param timestep        Simulation timestep [s]; 0 means auto-calculate from material properties.
 * @param B_ext           External magnetic field [T].
 * @param T_dust          Dust temperature [K]; -1 disables temperature corrections.
 * @param initial_state   Per-monomer initial state assembled from aggregate files.
 */
struct SimulationConfig {
    std::vector<AggregateConfig>    aggregates;
    std::vector<MaterialEntry>      materials;
    InitialState                    initial_state;
    OutputConfig                    output;
    int                             N_iter    = 0;
    double                          timestep  = 0.0;
    double3                         B_ext     = {0.0, 0.0, 0.0};
    double                          T_dust    = 15.0;
};
