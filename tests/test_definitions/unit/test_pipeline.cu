/**
 * @file test_pipeline.cu
 * @brief Integration tests for Pipeline::run() in simulationSetup/simulationSetup.cuh
 *
 * These tests exercise the full CPU setup path end-to-end:
 *   command file → aggregate loading → initial state → config validation
 *
 * Each test writes real temp files to /tmp, calls Pipeline::run(), and
 * inspects the resulting SimulationConfig via getConfig().
 *
 * NOTE: Pipeline::run() prints its normal console output (headline, log
 * messages, error text) during each test case. This is expected and does
 * not indicate a failure.
 */

#include "utils/config.cuh"                      // VERBOSITY — must come before printing.cuh
#include "simulationSetup/simulationSetup.cuh"
#include <cstdio>
#include <cmath>

static const char* PIPE_CMD  = "/tmp/test_pipe_cmd.tmp";
static const char* PIPE_AGG  = "/tmp/test_pipe_agg.tmp";
static const char* PIPE_AGG2 = "/tmp/test_pipe_agg2.tmp";

static void write_file(const char* path, const char* content) {
    FILE* f = fopen(path, "w");
    if (f) { fputs(content, f); fclose(f); }
}

// Writes a single-monomer aggregate at (x,y,z) nm with given radius (nm) and mat_id (1-indexed).
static void write_monomer_agg(const char* path,
                               double x_nm, double y_nm, double z_nm,
                               double r_nm, int mat_id) {
    char buf[256];
    snprintf(buf, sizeof(buf),
        "1 %.2f %.2f\n# m\n# m\n# m\n# m\n"
        "%.6f %.6f %.6f 0.0 %.6f 0.0 %d\n",
        r_nm * 2.0, r_nm,
        x_nm, y_nm, z_nm, r_nm, mat_id);
    write_file(path, buf);
}

void test_pipeline() {

    // ------------------------------------------------------------------ //
    // Wrong number of arguments → Status::error
    // ------------------------------------------------------------------ //

    {
        Pipeline p;
        const char* argv[] = {"test"};
        CHECK(p.run(1, argv) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Non-existent command file → Status::error
    // ------------------------------------------------------------------ //

    {
        Pipeline p;
        const char* argv[] = {"test", "/tmp/does_not_exist_pipeline_xyz.tmp"};
        CHECK(p.run(2, argv) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Valid setup: single monomer, explicit timestep.
    //
    // Checks: Status::ok, N_iter, N_save, timestep, initial_state sizes,
    //         position (aggregate offset applied), velocity, radius,
    //         mass = (4/3) π ρ r³, moment = (2/5) m r², mat_id 1→0.
    // ------------------------------------------------------------------ //

    {
        write_monomer_agg(PIPE_AGG, 0.0, 0.0, 0.0, 100.0, 1);

        char cmd[512];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 200\n"
            "<N_save> 20\n"
            "<path_A> \"%s\"\n"
            "<pos_A> 1e-6 2e-6 3e-6\n"
            "<vel_A> 0.0 0.0 -1.5\n"
            "<time_step> 1e-10\n"
            "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n",
            PIPE_AGG);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::ok);

        const SimulationConfig& cfg = p.getConfig();

        // Config-level fields
        CHECK(cfg.N_iter == 200);
        CHECK(cfg.output.N_save == 20);
        CHECK_APPROX(cfg.timestep, 1e-10, 1e-14);

        // initial_state arrays all have exactly 1 entry
        CHECK(cfg.initial_state.positions.size()    == 1);
        CHECK(cfg.initial_state.velocities.size()   == 1);
        CHECK(cfg.initial_state.radii.size()        == 1);
        CHECK(cfg.initial_state.masses.size()       == 1);
        CHECK(cfg.initial_state.moments.size()      == 1);
        CHECK(cfg.initial_state.material_ids.size() == 1);

        // Position: monomer at local (0,0,0) + aggregate offset (1μm,2μm,3μm)
        CHECK_APPROX(cfg.initial_state.positions[0].x, 1e-6, 1e-12);
        CHECK_APPROX(cfg.initial_state.positions[0].y, 2e-6, 1e-12);
        CHECK_APPROX(cfg.initial_state.positions[0].z, 3e-6, 1e-12);

        // Velocity: only aggregate velocity (no rotation)
        CHECK_APPROX(cfg.initial_state.velocities[0].x,  0.0, 1e-14);
        CHECK_APPROX(cfg.initial_state.velocities[0].y,  0.0, 1e-14);
        CHECK_APPROX(cfg.initial_state.velocities[0].z, -1.5, 1e-12);

        // Radius: 100 nm → 100e-9 m
        CHECK_APPROX(cfg.initial_state.radii[0], 100e-9, 1e-12);

        // Mass: (4/3) π ρ r³
        double r   = 100e-9;
        double rho = 2200.0;
        double m   = (4.0 / 3.0) * PI * rho * r * r * r;
        CHECK_APPROX(cfg.initial_state.masses[0], m, 1e-10);

        // Moment of inertia: (2/5) m r²
        CHECK_APPROX(cfg.initial_state.moments[0], (2.0 / 5.0) * m * r * r, 1e-10);

        // Material ID: 1-indexed in file → 0-indexed in state
        CHECK(cfg.initial_state.material_ids[0] == 0);
    }

    // ------------------------------------------------------------------ //
    // Tangential velocity from aggregate angular velocity.
    //
    // Monomer at local (d, 0, 0), aggregate angular ω = (0, 0, ω_z).
    //   vel_tang = ω × r_local = (0*0 - ω_z*0, ω_z*d - 0*0, 0) = (0, ω_z*d, 0)
    //   total velocity = agg_vel + vel_tang = (v_x, ω_z*d, 0)
    // ------------------------------------------------------------------ //

    {
        const double d_nm  = 50.0;
        const double d     = d_nm * 1e-9;
        const double omega = 1e6;

        write_monomer_agg(PIPE_AGG, d_nm, 0.0, 0.0, 50.0, 1);

        char cmd[512];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 100\n"
            "<N_save> 10\n"
            "<path_A> \"%s\"\n"
            "<pos_A> 0.0 0.0 1e-6\n"
            "<vel_A> 1.0 0.0 0.0\n"
            "<ang_A> 0.0 0.0 %.6e\n"
            "<time_step> 1e-10\n"
            "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n",
            PIPE_AGG, omega);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::ok);

        const SimulationConfig& cfg = p.getConfig();
        CHECK_APPROX(cfg.initial_state.velocities[0].x, 1.0,       1e-10);
        CHECK_APPROX(cfg.initial_state.velocities[0].y, omega * d,  1e-10);
        CHECK_APPROX(cfg.initial_state.velocities[0].z, 0.0,        1e-10);
    }

    // ------------------------------------------------------------------ //
    // Two aggregates: both monomers appear in initial_state in order.
    // ------------------------------------------------------------------ //

    {
        write_monomer_agg(PIPE_AGG,  0.0, 0.0, 0.0, 100.0, 1);
        write_monomer_agg(PIPE_AGG2, 0.0, 0.0, 0.0, 100.0, 1);

        char cmd[512];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 100\n"
            "<N_save> 10\n"
            "<path_A> \"%s\"\n"
            "<pos_A> 0.0 0.0  1e-6\n"
            "<vel_A> 0.0 0.0 -1.0\n"
            "<path_B> \"%s\"\n"
            "<pos_B> 0.0 0.0 -1e-6\n"
            "<vel_B> 0.0 0.0  1.0\n"
            "<time_step> 1e-10\n"
            "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n",
            PIPE_AGG, PIPE_AGG2);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::ok);

        const SimulationConfig& cfg = p.getConfig();
        CHECK(cfg.initial_state.positions.size() == 2);
        CHECK_APPROX(cfg.initial_state.positions[0].z,  1e-6, 1e-12);
        CHECK_APPROX(cfg.initial_state.positions[1].z, -1e-6, 1e-12);
        CHECK_APPROX(cfg.initial_state.velocities[0].z, -1.0, 1e-12);
        CHECK_APPROX(cfg.initial_state.velocities[1].z,  1.0, 1e-12);
    }

    // ------------------------------------------------------------------ //
    // Auto-timestep: when <time_step> is absent, it is calculated from the
    // minimum JKR contact timescale across all monomer pairs.
    // With 2 monomers the loop executes once; result must be > 0 and finite.
    // ------------------------------------------------------------------ //

    {
        write_monomer_agg(PIPE_AGG,  0.0, 0.0, 0.0, 100.0, 1);
        write_monomer_agg(PIPE_AGG2, 0.0, 0.0, 0.0, 100.0, 1);

        char cmd[512];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 100\n"
            "<N_save> 10\n"
            "<path_A> \"%s\"\n"
            "<pos_A> 0.0 0.0  1e-6\n"
            "<vel_A> 0.0 0.0 -1.0\n"
            "<path_B> \"%s\"\n"
            "<pos_B> 0.0 0.0 -1e-6\n"
            "<vel_B> 0.0 0.0  1.0\n"
            "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n",
            PIPE_AGG, PIPE_AGG2);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::ok);

        const SimulationConfig& cfg = p.getConfig();
        CHECK(cfg.timestep > 0.0);
        CHECK(std::isfinite(cfg.timestep));
    }

    // ------------------------------------------------------------------ //
    // Validation: N_iter = 0 → Status::error
    // ------------------------------------------------------------------ //

    {
        write_monomer_agg(PIPE_AGG, 0.0, 0.0, 0.0, 100.0, 1);

        char cmd[512];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 0\n"
            "<N_save> 10\n"
            "<path_A> \"%s\"\n"
            "<pos_A> 0.0 0.0 1e-6\n"
            "<vel_A> 0.0 0.0 -1.0\n"
            "<time_step> 1e-10\n"
            "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n",
            PIPE_AGG);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Validation: magnetic material with Msat non-zero but Tc = 0 → error
    // (Msat and Tc must both be zero or both non-zero)
    // ------------------------------------------------------------------ //

    {
        write_monomer_agg(PIPE_AGG, 0.0, 0.0, 0.0, 100.0, 1);

        char cmd[512];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 100\n"
            "<N_save> 10\n"
            "<path_A> \"%s\"\n"
            "<pos_A> 0.0 0.0 1e-6\n"
            "<vel_A> 0.0 0.0 -1.0\n"
            "<time_step> 1e-10\n"
            "<material id=\"1\"> \"iron\" 0.03 5e10 0.17 7800.0 2e-10 1e-9 1e-9 1e-8 1e5 0.5 0.0\n",
            PIPE_AGG);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Cleanup
    // ------------------------------------------------------------------ //

    remove(PIPE_CMD);
    remove(PIPE_AGG);
    remove(PIPE_AGG2);
}
