/**
 * @file test_parser.cu
 * @brief Tests for CommandFile::parse() in simulationSetup/commandFile.cuh
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * Each test writes a small command file to TMP_FILE, parses it, and checks
 * the resulting SimulationConfig.  Error-path tests will print error messages
 * from the parser to stdout — this is expected.
 */

#include "simulationSetup/commandFile.cuh"
#include <cstdio>

static const char* TMP_FILE = "/tmp/test_dustcollider_parser.tmp";

static void write_tmp(const char* content) {
    FILE* f = fopen(TMP_FILE, "w");
    if (f) { fputs(content, f); fclose(f); }
}

void test_parser() {

    // ------------------------------------------------------------------ //
    // File not found → Status::error (parser prints an error message here)
    // ------------------------------------------------------------------ //

    {
        CommandFile cf("/tmp/does_not_exist_dustcollider_xyz.tmp");
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Empty file → Status::ok, struct stays at defaults
    // ------------------------------------------------------------------ //

    {
        write_tmp("");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.N_iter == 0);
    }

    // ------------------------------------------------------------------ //
    // Comment-only file → Status::ok
    // ------------------------------------------------------------------ //

    {
        write_tmp("# this is a comment\n"
                  "# another comment\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
    }

    // ------------------------------------------------------------------ //
    // N_iter and N_save
    // ------------------------------------------------------------------ //

    {
        write_tmp("<N_iter> 500\n"
                  "<N_save> 50\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.N_iter == 500);
        CHECK(cfg.output.N_save == 50);
    }

    // ------------------------------------------------------------------ //
    // path_results: quoted path is stored correctly
    // ------------------------------------------------------------------ //

    {
        write_tmp("<path_results> \"/some/output/dir\"\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.output.path == "/some/output/dir");
    }

    // ------------------------------------------------------------------ //
    // B_ext: 3-component vector
    // ------------------------------------------------------------------ //

    {
        write_tmp("<B_ext> 0.0 0.0 1.5e-4\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK_APPROX(cfg.B_ext.x, 0.0,    1e-14);
        CHECK_APPROX(cfg.B_ext.y, 0.0,    1e-14);
        CHECK_APPROX(cfg.B_ext.z, 1.5e-4, 1e-12);
    }

    // ------------------------------------------------------------------ //
    // time_step and T_dust scalars
    // ------------------------------------------------------------------ //

    {
        write_tmp("<time_step> 1.0e-10\n"
                  "<T_dust> 300.0\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK_APPROX(cfg.timestep, 1.0e-10, 1e-23);
        CHECK_APPROX(cfg.T_dust,   300.0,   1e-12);
    }

    // ------------------------------------------------------------------ //
    // Boolean save flags: "1", "0", and "true"
    // ------------------------------------------------------------------ //

    {
        write_tmp("<save_pos> 1\n"
                  "<save_vel> 0\n"
                  "<save_energy> true\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK( cfg.output.position);
        CHECK(!cfg.output.velocity);
        CHECK( cfg.output.energy);
    }

    // ------------------------------------------------------------------ //
    // Aggregate position and velocity: two tags share one AggregateConfig
    // ------------------------------------------------------------------ //

    {
        write_tmp("<pos_A> 1.0 2.0 3.0\n"
                  "<vel_A> 0.0 0.0 -1.5\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.aggregates.size() == 1);
        CHECK(cfg.aggregates[0].name == "A");
        CHECK_APPROX(cfg.aggregates[0].position.x,  1.0,  1e-14);
        CHECK_APPROX(cfg.aggregates[0].position.y,  2.0,  1e-14);
        CHECK_APPROX(cfg.aggregates[0].position.z,  3.0,  1e-14);
        CHECK_APPROX(cfg.aggregates[0].velocity.z, -1.5,  1e-14);
    }

    // ------------------------------------------------------------------ //
    // Two aggregates (A and B) parsed independently
    // ------------------------------------------------------------------ //

    {
        write_tmp("<pos_A> 1.0 0.0 0.0\n"
                  "<pos_B> -1.0 0.0 0.0\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.aggregates.size() == 2);
        CHECK_APPROX(cfg.aggregates[0].position.x,  1.0, 1e-14);
        CHECK_APPROX(cfg.aggregates[1].position.x, -1.0, 1e-14);
    }

    // ------------------------------------------------------------------ //
    // Material tag (6 params, non-magnetic)
    // ------------------------------------------------------------------ //

    {
        write_tmp("<material id=\"1\"> \"steel\" 0.05 1e11 0.3 8000.0 1e-9 1e-8\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.materials.size() == 1);
        CHECK(cfg.materials[0].name == "steel");
        CHECK_APPROX(cfg.materials[0].gamma, 0.05,   1e-14);
        CHECK_APPROX(cfg.materials[0].E,     1e11,   1e-14);
        CHECK_APPROX(cfg.materials[0].nu,    0.3,    1e-14);
    }

    // ------------------------------------------------------------------ //
    // Unknown tag → Status::error (parser prints an error message here)
    // ------------------------------------------------------------------ //

    {
        write_tmp("<unknown_tag> 42\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Aggregate path and angular velocity tags
    // ------------------------------------------------------------------ //

    {
        write_tmp("<path_A> \"/some/agg/path\"\n"
                  "<ang_A> 0.1 0.2 0.3\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.aggregates.size() == 1);
        CHECK(cfg.aggregates[0].path == "/some/agg/path");
        CHECK_APPROX(cfg.aggregates[0].angular.x, 0.1, 1e-14);
        CHECK_APPROX(cfg.aggregates[0].angular.y, 0.2, 1e-14);
        CHECK_APPROX(cfg.aggregates[0].angular.z, 0.3, 1e-14);
    }

    // ------------------------------------------------------------------ //
    // Aggregate vector tags require exactly 3 components → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("<pos_A> 1.0 2.0\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // B_ext also requires exactly 3 components → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("<B_ext> 1.0 2.0\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Remaining boolean save flags: save_ovito, save_force, save_torque,
    // save_omega (save_pos/save_vel/save_energy are covered above)
    // ------------------------------------------------------------------ //

    {
        write_tmp("<save_ovito> 1\n"
                  "<save_force> 0\n"
                  "<save_torque> true\n"
                  "<save_omega> 0\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK( cfg.output.ovito);
        CHECK(!cfg.output.force);
        CHECK( cfg.output.torque);
        CHECK(!cfg.output.angular);
    }

    // ------------------------------------------------------------------ //
    // Invalid boolean literal → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("<save_pos> maybe\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Non-numeric scalar values → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("<N_iter> abc\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }
    {
        write_tmp("<N_save> abc\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }
    {
        write_tmp("<T_dust> abc\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }
    {
        write_tmp("<time_step> abc\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Legacy tags are accepted (with a warning) but have no effect
    // ------------------------------------------------------------------ //

    {
        write_tmp("<time_start> 0.0\n"
                  "<time_stop> 1.0e-4\n"
                  "<save_cluster> 1\n"
                  "<save_mag> 0\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
    }

    // ------------------------------------------------------------------ //
    // path_results / aggregate path with no quoted value → empty path
    // → Status::error.  A bare unquoted token is used so the tag is not
    // the last character on the line (avoids an unrelated out-of-range
    // substr() in extract_values when a tag has no trailing value at all).
    // ------------------------------------------------------------------ //

    {
        write_tmp("<path_results> not_a_path\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }
    {
        write_tmp("<path_A> not_a_path\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Regression: a tag with no value at all (the '>' is the very last
    // character on the sanitized line) must not crash. extract_values()
    // used to do line.substr(r_pos + 2), which throws std::out_of_range
    // here since r_pos + 2 > line.size() — uncaught, that would abort the
    // whole process instead of returning Status::error.
    // ------------------------------------------------------------------ //

    {
        write_tmp("<path_results>\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }
    {
        write_tmp("<N_iter>\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Material: magnetic variant (11 params)
    // ------------------------------------------------------------------ //

    {
        write_tmp("<material id=\"1\"> \"iron\" "
                  "0.03 5e10 0.17 7800.0 2e-10 1e-9 1e-9 1e-8 1e5 0.5 10.0\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.materials.size() == 1);
        CHECK_APPROX(cfg.materials[0].tss,  1e-9, 1e-14);
        CHECK_APPROX(cfg.materials[0].tsl,  1e-8, 1e-14);
        CHECK_APPROX(cfg.materials[0].Msat, 1e5,  1e-9);
        CHECK_APPROX(cfg.materials[0].chi,  0.5,  1e-14);
        CHECK_APPROX(cfg.materials[0].Tc,   10.0, 1e-12);
    }

    // ------------------------------------------------------------------ //
    // Material: no quoted name present at all → the (only) quoted token
    // is consumed as the name, leaving the ID empty → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("<material> \"steel\" 0.05 1e11 0.3 8000.0 1e-9 1e-8\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Material: non-numeric ID → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("<material id=\"abc\"> \"steel\" 0.05 1e11 0.3 8000.0 1e-9 1e-8\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Material: ID must be >= 1 → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("<material id=\"0\"> \"steel\" 0.05 1e11 0.3 8000.0 1e-9 1e-8\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Material: wrong parameter count (neither 6 nor 11) → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("<material id=\"1\"> \"steel\" 0.05 1e11 0.3 8000.0 1e-9\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Material: a sparse (out-of-order) ID resizes the materials vector,
    // leaving the skipped slots default-constructed
    // ------------------------------------------------------------------ //

    {
        write_tmp("<material id=\"3\"> \"thirdmat\" 0.05 1e11 0.3 8000.0 1e-9 1e-8\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::ok);
        CHECK(cfg.materials.size() == 3);
        CHECK(cfg.materials[0].name == "");
        CHECK(cfg.materials[1].name == "");
        CHECK(cfg.materials[2].name == "thirdmat");
    }

    // ------------------------------------------------------------------ //
    // Malformed line with no '<tag>' at all → Status::error
    // ------------------------------------------------------------------ //

    {
        write_tmp("no_tag_here 123\n");
        CommandFile cf(TMP_FILE);
        SimulationConfig cfg;
        CHECK(cf.parse(cfg) == Status::error);
    }

    // ------------------------------------------------------------------ //
    // Cleanup
    // ------------------------------------------------------------------ //

    remove(TMP_FILE);
}
