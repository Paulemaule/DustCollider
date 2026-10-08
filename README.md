# DustCollider

A CUDA N-body simulator for collisions between interstellar dust aggregates, composed of spherical 'monomers'.
Monomer-monomer contacts are modelled with JKR contact mechanics (Johnson et al. 1971) with extensions for sliding, rolling and twisting dofs (Dominik & Tielens 1997; Wada et al. 2007).
The state is advanced using a leapfrog integration algorithm on the GPU.
The code can handle multi-material aggregates.

## Requirements

- NVIDIA CUDA toolkit, C++20+ to compile using make (./Makefile).
- An NVIDIA GPU with compute capability of 8.9.
    The code potentially supports other compute capabilities if the appropriate flags are set in ./Makefile and ./tests/Makefile.

## Build

```bash
make                  # Release build (default)
make BUILD=Debug      # Debug build (-g -G -DDEBUG)
make test             # To run the test suite
make clean            # Remove build artifacts
```

The executable is written to `build/dust_collider`.

## Run

A run can be started by providing the executable with a command file.

```bash
./build/dust_collider <path/to/cmd_file>
```

The code writes snapshots of the systems state to disk on completion.
Logging to the command line contains vital information about the run and should be redirected to a log file for later analysis.

## Command files

Command files contain the full configuration for a single dust collider run.
They are plain text and specify run parameters in `<tag>` value pairs.
Comments in the command file are possible by prepending `#` or `!`.
Relative paths are resolved relative to the command files location.

```text
# Aggregate specification
<path_A>            "./aggregate_A"   # input aggregate file
<path_B>            "./aggregate_B"
<pos_A>             -50.1e-9 0 0      # centre-of-mass position [m]
<pos_B>             +50.1e-9 0 0
<vel_A>             +4 0 0            # bulk velocity [m/s]
<vel_B>             -4 0 0
<ang_A>             0 0 0             # angular velocity [rad/s]
<ang_B>             0 0 0

# Run specification (alternatives in physical time: see "Run length and snapshots" below)
<N_iter>            1000000           # number of integration steps
<N_save>            10000             # store a snapshot every N_save steps
<time_step>         0                 # [s]; Recommended to be omitted, the code determines the optimal timestep based on the material properties.

# Output specification
<save_ovito>        1                 # write OVITO .dump files
<save_pos>          1
<save_vel>          1
<save_omega>        1
<save_force>        1
<save_torque>       1
<path_results>      "./out/"          # output directory

# Material properties
# <material id="N">  "name"  gamma  E  nu  rho  xi  tvis  [tss tsl Msat chi Tc]
<material id="1">   "forsterite"  0.07  204e9  0.24  3210  2e-10  1e-12
```

The material parameters are: The surface energy `gamma` [J/m²] ; Youngs modulus `E` [Pa] ; Poisson ratio `nu` ; the density `rho` [kg/m³] ; the critical rolling displacement `xi` [m] (see eg Wada et al. 2007) ; and the viscous damping timescale `tvis` [s].
A material may also define five magnetic parameters (`tss`, `tsl`, `Msat`, `chi`, `Tc`), however they are currently unused.

### Run length and snapshots

The auto-calculated timestep depends on the materials and monomer sizes. The effects of several tags in the command file depend on the exact result of the timestep calculation.

| Group | Tag | Meaning |
|---|---|---|
| run length | `<N_iter>` | number of integration steps |
| | `<t_end>` | simulated time [s], converted to `ceil(t_end / dt)` steps |
| snapshots | `<N_save>` | a snapshot every `N_save` steps (constant cadence in units of the contact timescale when dt is auto-calculated) |
| | `<t_save>` | a snapshot every `t_save` seconds, converted to `ceil(t_save / dt)` steps (comparable across setups) |
| | `<N_snap>` | the total number of snapshots, including the initial and the final state (constant output size) |

Exactly one tag of each group is required.
Snapshots are stored on a uniform grid: snapshot `k` is the state after `k * N_save` steps, i.e. at `t = k * N_save * dt`, and snapshot 0 is the initial state.
The run is extended to a whole number of snapshot intervals, so the final state is always stored; the setup log warns when this happens.
Integer values accept scientific notation (`<N_iter> 1e6`).
The setup log reports the resolved timestep, `N_iter` and `N_save`, and the run time and snapshot cadence in seconds and in units of the smallest contact timescale `tau_N,min`.

Aggregates are identified by the key after the underscore in `<path_X>`, `<pos_X>`, `<vel_X>` and `<ang_X>` (`A`, `B`, …), any number of aggregates should be supported.

## Aggregate file format

Aggregate files specify the initial configuration of a single aggregate.
They are plain text and determine the position, radius and material of each monomer.
Quantities are whitespace separated.
Line 0 is a header containing (The number of monomers, external aggregate radius [nm], effective radius [nm]).
Lines 1-4 contain additional metadata.
Lines 5+ contain information on a single monomer `x y z _ radius _ mat_id` (positions [nm], radius [nm] and material ID 1-indexed)

## Output

Written to `<path_results>/`:

- `binary/` — raw `double`/`int` arrays of the stored snapshots (positions,
  velocities, forces, torques, angular velocities) plus energies and cluster IDs.
  There are `N_iter / N_save + 1` snapshots, snapshot `k` at `t = k * N_save * dt`.
  Each energy file holds one `double` per snapshot, summed over all monomers:
  - `sim_{normal,sliding,rolling,twisting}_pot.bin` — the potential energy stored in
    the contacts, averaged over the `N_save` steps before the snapshot.
  - `sim_normal_damp.bin` — viscous damping of the normal motion.
  - `sim_{sliding,rolling,twisting}_slip.bin` — energy dissipated through inelastic motion beyond critical tangential displacements.
  - `sim_{normal,sliding,rolling,twisting}_break.bin` — when a contact breaks all stored potential energy is dissipated.
  - `sim_normal_form.bin` — when contacts are formed, energy is dissipated.

- `ovito/` — `t_*.dump` files for visualization in OVITO (https://www.ovito.org/)
  (vectors are rescaled for display).