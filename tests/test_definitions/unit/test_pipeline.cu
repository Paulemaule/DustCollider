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
static const char* PIPE_AGG3 = "/tmp/test_pipe_agg3.tmp";

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

// Runs the pipeline on a valid two-monomer setup, completed by the given run length and snapshot tags.
static Status run_schedule_setup(Pipeline& p, const char* schedule_tags) {
    write_monomer_agg(PIPE_AGG,  0.0, 0.0, 0.0, 100.0, 1);
    write_monomer_agg(PIPE_AGG2, 0.0, 0.0, 0.0, 100.0, 1);

    char cmd[1024];
    snprintf(cmd, sizeof(cmd),
        "%s"
        "<path_A> \"%s\"\n"
        "<pos_A> 0.0 0.0  1e-6\n"
        "<vel_A> 0.0 0.0 -1.0\n"
        "<path_B> \"%s\"\n"
        "<pos_B> 0.0 0.0 -1e-6\n"
        "<vel_B> 0.0 0.0  1.0\n"
        "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n",
        schedule_tags, PIPE_AGG, PIPE_AGG2);
    write_file(PIPE_CMD, cmd);

    const char* argv[] = {"test", PIPE_CMD};
    return p.run(2, argv);
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
    // The spin acts about the center of mass, not the file origin.
    //
    // A single monomer at local (d, 0, 0) is its own center of mass, so the aggregate angular velocity ω = (0, 0, ω_z)
    // gives it no velocity: it moves with the aggregate velocity alone and spins with ω.
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
        CHECK_APPROX(cfg.aggregates[0].center_of_mass.x / d, 1.0, 1e-15);
        CHECK_APPROX(cfg.initial_state.velocities[0].x, 1.0, 1e-15);
        CHECK_APPROX(cfg.initial_state.velocities[0].y, 0.0, 1e-15);
        CHECK_APPROX(cfg.initial_state.velocities[0].z, 0.0, 1e-15);
        CHECK(cfg.initial_state.omegas[0].z == omega);
    }

    // ------------------------------------------------------------------ //
    // Rigid motion of an aggregate whose center of mass is neither the file origin nor the centroid.
    //
    // Three monomers of two densities, off the file origin, with an angular velocity ω along all three axes.
    // Each monomer moves with v + ω × (r_k - c), with c the center of mass in the file frame, so the center of mass
    // moves with v. The positions are the file positions offset by pos.
    // ------------------------------------------------------------------ //

    {
        const double x_nm[3][3] = { { -20.0,  5.0,  10.0 },
                                    {  15.0, -8.0,   3.0 },
                                    {   4.0, 25.0, -12.0 } };
        const double r_nm[3]    = { 10.0, 20.0, 15.0 };
        const double rho[3]     = { 2200.0, 5000.0, 2200.0 };
        const double pos[3]     = { 1e-6, -2e-6, 5e-7 };
        const double vel[3]     = { 1.5, -0.5, 2.0 };
        const double w[3]       = { 3e7, -2e7, 5e7 };

        write_file(PIPE_AGG, "3 80.0 30.0\n# m\n# m\n# m\n# m\n"
                             "-20.0  5.0  10.0 0.0 10.0 0.0 1\n"
                             " 15.0 -8.0   3.0 0.0 20.0 0.0 2\n"
                             "  4.0 25.0 -12.0 0.0 15.0 0.0 1\n");

        char cmd[1024];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 100\n"
            "<N_save> 10\n"
            "<path_A> \"%s\"\n"
            "<pos_A> 1e-6 -2e-6 5e-7\n"
            "<vel_A> 1.5 -0.5 2.0\n"
            "<ang_A> 3e7 -2e7 5e7\n"
            "<time_step> 1e-10\n"
            "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n"
            "<material id=\"2\"> \"dense\"  0.03 5e10 0.17 5000.0 2e-10 1e-9\n",
            PIPE_AGG);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::ok);

        const SimulationConfig& cfg = p.getConfig();
        const InitialState&     is  = cfg.initial_state;

        // The center of mass in the file frame.
        double mass = 0.0;
        double c[3] = { 0.0, 0.0, 0.0 };
        for (int k = 0; k < 3; k++) {
            const double r = r_nm[k] * 1e-9;
            const double m = (4.0 / 3.0) * PI * rho[k] * r * r * r;
            mass += m;
            for (int j = 0; j < 3; j++) c[j] += m * x_nm[k][j] * 1e-9;
        }
        for (int j = 0; j < 3; j++) c[j] /= mass;

        CHECK_APPROX(cfg.aggregates[0].center_of_mass.x / c[0], 1.0, 1e-13);
        CHECK_APPROX(cfg.aggregates[0].center_of_mass.y / c[1], 1.0, 1e-13);
        CHECK_APPROX(cfg.aggregates[0].center_of_mass.z / c[2], 1.0, 1e-13);

        double momentum[3] = { 0.0, 0.0, 0.0 };
        for (int k = 0; k < 3; k++) {
            const double d[3] = { x_nm[k][0] * 1e-9 - c[0], x_nm[k][1] * 1e-9 - c[1], x_nm[k][2] * 1e-9 - c[2] };

            CHECK_APPROX(is.velocities[k].x, vel[0] + w[1] * d[2] - w[2] * d[1], 1e-14);
            CHECK_APPROX(is.velocities[k].y, vel[1] + w[2] * d[0] - w[0] * d[2], 1e-14);
            CHECK_APPROX(is.velocities[k].z, vel[2] + w[0] * d[1] - w[1] * d[0], 1e-14);

            CHECK_APPROX(is.positions[k].x / (pos[0] + x_nm[k][0] * 1e-9), 1.0, 1e-15);
            CHECK_APPROX(is.positions[k].y / (pos[1] + x_nm[k][1] * 1e-9), 1.0, 1e-15);
            CHECK_APPROX(is.positions[k].z / (pos[2] + x_nm[k][2] * 1e-9), 1.0, 1e-15);

            CHECK(is.omegas[k].x == w[0] && is.omegas[k].y == w[1] && is.omegas[k].z == w[2]);

            momentum[0] += is.masses[k] * is.velocities[k].x;
            momentum[1] += is.masses[k] * is.velocities[k].y;
            momentum[2] += is.masses[k] * is.velocities[k].z;
        }

        // The center of mass moves with the aggregate velocity.
        CHECK_APPROX(momentum[0] / cfg.aggregates[0].mass, vel[0], 1e-14);
        CHECK_APPROX(momentum[1] / cfg.aggregates[0].mass, vel[1], 1e-14);
        CHECK_APPROX(momentum[2] / cfg.aggregates[0].mass, vel[2], 1e-14);
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
    // Run schedule (issue Z): the run length and snapshot tags are resolved
    // into N_iter and N_save. Snapshots lie on a uniform grid of N_save
    // iterations and N_iter is rounded up to a whole number of intervals.
    // An explicit <time_step> of 1e-12 s keeps the expected counts exact.
    // ------------------------------------------------------------------ //

    struct ScheduleCase {
        const char* tags;
        int         N_iter;
        int         N_save;
    };

    for (const ScheduleCase& c : {
            // <N_iter> + <N_save> that already fit are left unchanged
            ScheduleCase{ "<N_iter> 1000\n<N_save> 100\n",                      1000,  100 },
            // An explicit <N_iter> is extended to end on a snapshot
            ScheduleCase{ "<N_iter> 1003\n<N_save> 100\n",                      1100,  100 },
            // <N_save> larger than <N_iter> gives a single interval
            ScheduleCase{ "<N_iter> 100\n<N_save> 1000\n",                      1000, 1000 },
            // ceil(1e-9 / 1e-12) is 1001 in floating point, which must not extend the run
            ScheduleCase{ "<time_step> 1e-12\n<t_end> 1e-9\n<N_save> 100\n",    1000,  100 },
            ScheduleCase{ "<time_step> 1e-12\n<t_end> 1.0005e-9\n<N_save> 100\n", 1100, 100 },
            // <t_save> is rounded up to whole iterations, with the same tolerance as <t_end>
            ScheduleCase{ "<time_step> 1e-12\n<N_iter> 1000\n<t_save> 1e-10\n",      1000, 100 },
            ScheduleCase{ "<time_step> 1e-12\n<N_iter> 1000\n<t_save> 1.0004e-10\n", 1010, 101 },
            ScheduleCase{ "<time_step> 1e-12\n<N_iter> 1000\n<t_save> 1e-12\n",      1000,   1 },
            ScheduleCase{ "<time_step> 1e-12\n<t_end> 1e-9\n<t_save> 2.5e-10\n",     1000, 250 },
            // <t_end> = 50 * <t_save> gives exactly 50 intervals. Rounding <t_save> to nearest (100 iterations)
            // made the run 22 iterations short of 50 intervals and added a spurious 51st (5100 iterations).
            ScheduleCase{ "<time_step> 1e-12\n<t_end> 5.002e-9\n<t_save> 1.0004e-10\n", 5050, 101 },
            // <N_snap> counts both the initial and the final state
            ScheduleCase{ "<N_iter> 1000\n<N_snap> 11\n",                       1000,  100 },
            ScheduleCase{ "<N_iter> 1001\n<N_snap> 11\n",                       1010,  101 },
            ScheduleCase{ "<N_iter> 5\n<N_snap> 11\n",                            10,    1 },
            ScheduleCase{ "<time_step> 1e-12\n<t_end> 1e-9\n<N_snap> 2\n",      1000, 1000 } }) {
        Pipeline p;
        CHECK(run_schedule_setup(p, c.tags) == Status::ok);
        CHECK(p.getConfig().N_iter        == c.N_iter);
        CHECK(p.getConfig().output.N_save == c.N_save);
    }

    // ------------------------------------------------------------------ //
    // Run schedule: invalid combinations and values → Status::error
    // ------------------------------------------------------------------ //

    for (const char* tags : {
            // Run length: none or both
            "<N_save> 100\n",
            "<time_step> 1e-12\n<N_iter> 1000\n<t_end> 1e-9\n<N_save> 100\n",
            // Snapshots: none or more than one
            "<N_iter> 1000\n",
            "<N_iter> 1000\n<N_save> 100\n<N_snap> 11\n",
            "<time_step> 1e-12\n<N_iter> 1000\n<N_save> 100\n<t_save> 1e-10\n",
            // Invalid values
            "<N_iter> -1000\n<N_save> 100\n",
            "<N_iter> 1000\n<N_save> -100\n",
            "<time_step> 1e-12\n<t_end> -1e-9\n<N_save> 100\n",
            "<time_step> 1e-12\n<N_iter> 1000\n<t_save> -1e-10\n",
            "<time_step> 1e-12\n<N_iter> 1000\n<t_save> 9.99e-13\n",       // shorter than dt
            "<time_step> 1e-12\n<t_end> 1e-2\n<N_save> 100\n",             // 1e10 iterations
            "<N_iter> 1000\n<N_snap> 1\n",
            // Rounding up to whole intervals exceeds the range of int
            "<N_iter> 2147483001\n<N_save> 1000\n" }) {
        Pipeline p;
        CHECK(run_schedule_setup(p, tags) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Run schedule with the auto-calculated timestep: the run covers at
    // least <t_end>, extended by less than one <N_snap> interval.
    // ------------------------------------------------------------------ //

    {
        Pipeline p;
        CHECK(run_schedule_setup(p, "<t_end> 1e-9\n<N_snap> 101\n") == Status::ok);

        const SimulationConfig& cfg = p.getConfig();
        const double t_run = cfg.N_iter * cfg.timestep;
        CHECK(cfg.tau_min > 0.0);
        CHECK(cfg.N_iter % 100 == 0);
        CHECK(cfg.N_iter / cfg.output.N_save == 100);
        CHECK(t_run >= 1e-9 * (1.0 - 1e-9));
        CHECK(t_run <  1e-9 + 100 * cfg.timestep);
    }

    // ------------------------------------------------------------------ //
    // tau_min is also calculated when <time_step> is set explicitly, and the source of the timestep is recorded
    // ------------------------------------------------------------------ //

    {
        Pipeline pa, pb;
        CHECK(run_schedule_setup(pa, "<N_iter> 1000\n<N_save> 100\n") == Status::ok);
        CHECK(run_schedule_setup(pb, "<time_step> 1e-12\n<N_iter> 1000\n<N_save> 100\n") == Status::ok);
        CHECK(pb.getConfig().tau_min > 0.0);
        CHECK(pa.getConfig().tau_min == pb.getConfig().tau_min);
        CHECK(pa.getConfig().timestep == 0.005 * pa.getConfig().tau_min);
        CHECK(!pa.getConfig().timestep_explicit);
        CHECK(pb.getConfig().timestep_explicit);
    }

    // ------------------------------------------------------------------ //
    // Pair scales of a single pair: both extremes hold its value, assembled from the helpers of integrator_utils.cuh
    // ------------------------------------------------------------------ //

    {
        Pipeline p;
        CHECK(run_schedule_setup(p, "<N_iter> 1000\n<N_save> 100\n") == Status::ok);
        const SimulationConfig& cfg = p.getConfig();
        const PairScales&       s   = cfg.pair_scales;

        // Two silica monomers of 100 nm: gamma 0.03, E 5e10, nu 0.17, xi 2e-10.
        const double R      = get_R(100e-9, 100e-9);
        const double M      = 0.5 * cfg.initial_state.masses[0];
        const double gamma  = get_gamma(0.03, 0.03);
        const double G      = get_G_i(5e10, 0.17);
        const double a_0    = get_a_0(gamma, R, get_E_s(5e10, 5e10, 0.17, 0.17));
        const double F_c    = get_F_c(gamma, R);
        const double d_c    = get_delta_N_crit(a_0, R);
        const double k_s    = get_k_s(get_G_s(G, G, 0.17, 0.17), a_0);
        const double k_r    = get_k_r(F_c, R);
        const double k_t    = get_k_t(get_G(G, G), a_0);
        const double d_S    = get_delta_S_crit(0.17, 0.17, a_0);
        const double d_T    = get_delta_T_crit();
        const double E_roll = 0.5 * PI * R * k_r * 2e-10;

        CHECK(s.F_c.min == s.F_c.max);
        CHECK(s.F_c.min_i == 0 && s.F_c.min_j == 1 && s.F_c.max_i == 0 && s.F_c.max_j == 1);

        CHECK_APPROX(s.a_0.min          / a_0,                           1.0,  1e-14);
        CHECK_APPROX(s.delta_N_0.min    / get_delta_N_0(a_0, R),         1.0,  1e-14);
        CHECK_APPROX(s.delta_N_crit.min / d_c,                           1.0,  1e-14);
        CHECK_APPROX(s.F_c.min          / F_c,                           1.0,  1e-14);
        CHECK_APPROX(s.tau_N.min        / sqrt(M * d_c / F_c),           1.0,  1e-14);
        CHECK(s.tau_N.min == cfg.tau_min);

        // The JKR potential between equilibrium and the turning point, 1.54 F_c delta_N_crit in Wada et al. (2007).
        CHECK_APPROX(s.E_break.min / (F_c * d_c),                        1.5425853631546,  1e-9);
        // E_roll = 6 pi^2 gamma R xi of Dominik & Tielens (1997).
        CHECK_APPROX(s.E_roll.min  / E_roll,                             1.0,  1e-14);
        CHECK_APPROX(s.E_roll.min  / (6. * PI * PI * gamma * R * 2e-10), 1.0,  1e-14);

        CHECK_APPROX(s.U_S_crit.min / (0.5 * k_s * d_S * d_S),           1.0,  1e-14);
        CHECK_APPROX(s.U_R_crit.min / (0.5 * k_r * 2e-10 * 2e-10),       1.0,  1e-14);
        CHECK_APPROX(s.U_T_crit.min / (0.5 * k_t * d_T * d_T),           1.0,  1e-14);

        CHECK_APPROX(s.v_break.min / sqrt(2. * s.E_break.min / M),       1.0,  1e-14);
        CHECK_APPROX(s.v_roll.min  / sqrt(2. * E_roll / M),              1.0,  1e-14);
    }

    // ------------------------------------------------------------------ //
    // Pair scales of three monomers of one material with radii 50, 100 and 200 nm: the scales that grow with the
    // reduced radius are smallest for the pair (0, 1) and largest for the pair (1, 2)
    // ------------------------------------------------------------------ //

    {
        write_monomer_agg(PIPE_AGG,  0.0, 0.0, 0.0,  50.0, 1);
        write_monomer_agg(PIPE_AGG2, 0.0, 0.0, 0.0, 100.0, 1);
        write_monomer_agg(PIPE_AGG3, 0.0, 0.0, 0.0, 200.0, 1);

        char cmd[1024];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 1000\n<N_save> 100\n"
            "<path_A> \"%s\"\n<pos_A> 0.0 0.0  1e-6\n"
            "<path_B> \"%s\"\n<pos_B> 0.0 0.0 -1e-6\n"
            "<path_C> \"%s\"\n<pos_C> 0.0 1e-6 0.0\n"
            "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n",
            PIPE_AGG, PIPE_AGG2, PIPE_AGG3);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::ok);
        const PairScales& s = p.getConfig().pair_scales;

        for (const Extremes* e : { &s.a_0, &s.F_c, &s.E_roll }) {
            CHECK(e->min < e->max);
            CHECK(e->min_i == 0 && e->min_j == 1);
            CHECK(e->max_i == 1 && e->max_j == 2);
        }
        CHECK_APPROX(s.F_c.min / get_F_c(0.03, get_R( 50e-9, 100e-9)),  1.0,  1e-14);
        CHECK_APPROX(s.F_c.max / get_F_c(0.03, get_R(100e-9, 200e-9)),  1.0,  1e-14);
        CHECK(s.tau_N.min == p.getConfig().tau_min);
    }

    // ------------------------------------------------------------------ //
    // Aggregate bookkeeping: the monomer range of each aggregate, the header radii and the center of mass in the file
    // frame
    // ------------------------------------------------------------------ //

    {
        // Aggregate A: monomers of 10 and 20 nm at x = -10 and 20 nm. Aggregate B: a single monomer.
        write_file(PIPE_AGG, "2 60.0 20.8\n# m\n# m\n# m\n# m\n"
                             "-10.0 0.0 0.0 10.0 10.0 0.0 1\n"
                             " 20.0 0.0 0.0 20.0 20.0 0.0 1\n");
        write_monomer_agg(PIPE_AGG2, 0.0, 0.0, 0.0, 100.0, 1);

        char cmd[1024];
        snprintf(cmd, sizeof(cmd),
            "<N_iter> 1000\n<N_save> 100\n"
            "<path_A> \"%s\"\n<pos_A> 0.0 0.0  1e-6\n"
            "<path_B> \"%s\"\n<pos_B> 0.0 0.0 -1e-6\n"
            "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n",
            PIPE_AGG, PIPE_AGG2);
        write_file(PIPE_CMD, cmd);

        Pipeline p;
        const char* argv[] = {"test", PIPE_CMD};
        CHECK(p.run(2, argv) == Status::ok);
        const AggregateConfig& A = p.getConfig().aggregates[0];
        const AggregateConfig& B = p.getConfig().aggregates[1];

        CHECK(A.first_monomer == 0 && A.Nmon == 2);
        CHECK(B.first_monomer == 2 && B.Nmon == 1);
        CHECK_APPROX(A.external_radius  / 60e-9,   1.0,  1e-14);
        CHECK_APPROX(A.effective_radius / 20.8e-9, 1.0,  1e-14);

        // One material, so the masses go with r^3: (1 * (-10) + 8 * 20) / 9 nm.
        CHECK_APPROX(A.mass / (4.0 / 3.0 * PI * 2200.0 * (1e-24 + 8e-24)), 1.0,  1e-14);
        CHECK_APPROX(A.center_of_mass.x / (150.0 / 9.0 * 1e-9), 1.0,  1e-14);
        CHECK(A.center_of_mass.y == 0.0 && A.center_of_mass.z == 0.0);
        CHECK(B.center_of_mass.x == 0.0 && B.center_of_mass.y == 0.0 && B.center_of_mass.z == 0.0);
    }

    // ------------------------------------------------------------------ //
    // Cleanup
    // ------------------------------------------------------------------ //

    remove(PIPE_CMD);
    remove(PIPE_AGG);
    remove(PIPE_AGG2);
}
