/**
 * @file test_output.cu
 * @brief Tests for Simulator::write_output() in simulator.cuh: the files it writes, the formats header.txt gives
 * for them, and the status it returns when the output directory or a file can not be written.
 *
 * Each case sets up a two-monomer run with Pipeline::run() and constructs the Simulator, which allocates device
 * memory, but never runs it. write_output() then writes the zero-initialised snapshot arrays.
 *
 * The failures are provoked by files and directories that block the output paths, which fails for root as well.
 * /dev/full stands in for a full disk: opening and writing succeed, and only the flush on closing fails.
 *
 * NOTE: Pipeline::run() and the Simulator print their normal console output, including the error messages of the
 * failing cases. This is expected and does not indicate a failure.
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 */

#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>     // getpid

#include "utils/config.cuh"                      // VERBOSITY — must come before printing.cuh
#include "simulationSetup/simulationSetup.cuh"
#include "simulator.cuh"

// Per process, so that concurrent runs on a shared machine do not collide.
static const std::filesystem::path OUT_ROOT =
    std::filesystem::temp_directory_path() / ("test_output_" + std::to_string(getpid()));

static const size_t OUT_NMON  = 2;      // monomers of the run
static const size_t OUT_NSNAP = 2;      // snapshots of the run: N_iter = N_save

// Writes a text file, creating its directory.
static void out_write_text(const std::filesystem::path& path, const std::string& content) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream(path) << content;
}

// The size of a file, or -1 if it does not exist.
static uintmax_t out_size(const std::filesystem::path& path) {
    std::error_code ec;
    const uintmax_t size = std::filesystem::file_size(path, ec);
    return ec ? uintmax_t(-1) : size;
}

// Reads binary/header.txt: the dtype of each listed file, and the size in bytes that its dtype and shape give.
struct OutFormat { std::string dtype; uintmax_t bytes; };
static std::map<std::string, OutFormat> out_read_header(const std::filesystem::path& path) {
    std::map<std::string, OutFormat> formats;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream tokens(line);
        std::string name, dtype;
        tokens >> name >> dtype;

        // The shape is the comma separated list between the brackets
        const size_t open = line.find('['), close = line.find(']');
        if (open == std::string::npos || close == std::string::npos) continue;
        uintmax_t count = 1;
        std::istringstream dims(line.substr(open + 1, close - open - 1));
        std::string dim;
        while (std::getline(dims, dim, ',')) count *= std::stoull(dim);

        const uintmax_t item = dtype == "float64" ? 8 : dtype == "int32" ? 4 : dtype == "uint8" ? 1 : 0;
        formats[name] = {dtype, count * item};
    }
    return formats;
}

// Sets up a run of two single-monomer aggregates that saves every state array and the Ovito files into out_dir,
// and returns the status of write_output().
static Status out_write(const std::filesystem::path& out_dir) {
    const std::filesystem::path agg = OUT_ROOT / "aggregate";
    out_write_text(agg, "1 200.0 100.0\n# m\n# m\n# m\n# m\n0.0 0.0 0.0 0.0 100.0 0.0 1\n");

    const std::filesystem::path cmd = OUT_ROOT / "cmd_file";
    out_write_text(cmd,
        "<path_A> \"" + agg.string() + "\"\n"
        "<pos_A> 0.0 0.0  1e-6\n"
        "<vel_A> 0.0 0.0 -1.0\n"
        "<path_B> \"" + agg.string() + "\"\n"
        "<pos_B> 0.0 0.0 -1e-6\n"
        "<vel_B> 0.0 0.0  1.0\n"
        "<material id=\"1\"> \"silica\" 0.03 5e10 0.17 2200.0 2e-10 1e-9\n"
        "<N_iter> 10\n"
        "<N_save> 10\n"
        "<save_pos> 1\n"
        "<save_vel> 1\n"
        "<save_omega> 1\n"
        "<save_force> 1\n"
        "<save_torque> 1\n"
        "<save_ovito> 1\n"
        "<path_results> \"" + out_dir.string() + "\"\n");

    Pipeline p;
    const std::string cmd_path = cmd.string();
    const char* argv[] = {"test", cmd_path.c_str()};
    const Status setup = p.run(2, argv);
    CHECK(setup == Status::ok);
    if (setup != Status::ok) return Status::error;

    Simulator sim(p.getConfig());
    return sim.write_output();
}

void test_output() {
    namespace fs = std::filesystem;

    const uintmax_t vec3_bytes   = OUT_NSNAP * OUT_NMON * sizeof(double3);
    const uintmax_t energy_bytes = OUT_NSNAP * sizeof(double);

    // ------------------------------------------------------------------ //
    // Writable directory → Status::ok, every file with its full size
    // ------------------------------------------------------------------ //

    {
        fs::remove_all(OUT_ROOT);
        const fs::path out = OUT_ROOT / "ok";
        const fs::path bin = out / "binary";

        CHECK(out_write(out) == Status::ok);

        CHECK(out_size(bin / "agg_a_mon.bin")     == OUT_NMON * sizeof(double));
        CHECK(out_size(bin / "agg_mass_mon.bin")  == OUT_NMON * sizeof(double));
        CHECK(out_size(bin / "agg_matid_mon.bin") == OUT_NMON * sizeof(int));
        CHECK(out_size(bin / "sim_cluster.bin")   == OUT_NSNAP * OUT_NMON * sizeof(int));
        CHECK(out_size(bin / "sim_contact.bin")   == OUT_NSNAP * contactMatrixBytes(OUT_NMON));

        for (const char* name : {"sim_pos.bin", "sim_vel.bin", "sim_omega.bin", "sim_force.bin", "sim_torque.bin"}) {
            CHECK(out_size(bin / name) == vec3_bytes);
        }

        for (const char* name : {"sim_normal_pot.bin",   "sim_sliding_pot.bin",   "sim_rolling_pot.bin",
                                 "sim_twisting_pot.bin", "sim_normal_damp.bin",   "sim_sliding_slip.bin",
                                 "sim_rolling_slip.bin", "sim_twisting_slip.bin", "sim_normal_break.bin",
                                 "sim_sliding_break.bin", "sim_rolling_break.bin", "sim_twisting_break.bin",
                                 "sim_normal_form.bin"}) {
            CHECK(out_size(bin / name) == energy_bytes);
        }

        CHECK(fs::is_regular_file(out / "ovito" / "t_00000.dump"));
        CHECK(fs::is_regular_file(out / "ovito" / "t_00001.dump"));

        // header.txt lists every other file of the directory, with the dtype and shape of its size
        const auto formats = out_read_header(bin / "header.txt");
        CHECK(formats.size() == 23);
        for (const auto& [name, format] : formats) {
            CHECK(format.bytes > 0);
            CHECK(out_size(bin / name) == format.bytes);
        }
        for (const auto& entry : fs::directory_iterator(bin)) {
            const std::string name = entry.path().filename().string();
            if (name != "header.txt") CHECK(formats.count(name) == 1);
        }
        CHECK(formats.at("sim_pos.bin").dtype     == "float64");
        CHECK(formats.at("sim_cluster.bin").dtype == "int32");
        CHECK(formats.at("sim_contact.bin").dtype == "uint8");
    }

    // ------------------------------------------------------------------ //
    // Output directory below a regular file → Status::error
    // ------------------------------------------------------------------ //

    {
        fs::remove_all(OUT_ROOT);
        out_write_text(OUT_ROOT / "blocker", "x");

        CHECK(out_write(OUT_ROOT / "blocker" / "run") == Status::error);
    }

    // ------------------------------------------------------------------ //
    // A binary file blocked by a directory → Status::error,
    // the other files are still written
    // ------------------------------------------------------------------ //

    {
        fs::remove_all(OUT_ROOT);
        const fs::path out = OUT_ROOT / "blocked_file";
        fs::create_directories(out / "binary" / "sim_pos.bin");

        CHECK(out_write(out) == Status::error);
        CHECK(out_size(out / "binary" / "sim_vel.bin")         == vec3_bytes);
        CHECK(out_size(out / "binary" / "sim_normal_form.bin") == energy_bytes);
        CHECK(fs::is_regular_file(out / "ovito" / "t_00001.dump"));

        // header.txt lists only the files that were written
        const auto formats = out_read_header(out / "binary" / "header.txt");
        CHECK(formats.size() == 22);
        CHECK(formats.count("sim_pos.bin") == 0);
        CHECK(formats.count("sim_vel.bin") == 1);
    }

    // ------------------------------------------------------------------ //
    // A binary file on a full disk → Status::error
    // ------------------------------------------------------------------ //

    if (fs::exists("/dev/full")) {
        fs::remove_all(OUT_ROOT);
        const fs::path out = OUT_ROOT / "full_disk";
        fs::create_directories(out / "binary");
        fs::create_symlink("/dev/full", out / "binary" / "sim_normal_pot.bin");

        CHECK(out_write(out) == Status::error);
        CHECK(out_size(out / "binary" / "sim_sliding_pot.bin") == energy_bytes);

        const auto formats = out_read_header(out / "binary" / "header.txt");
        CHECK(formats.count("sim_normal_pot.bin")  == 0);
        CHECK(formats.count("sim_sliding_pot.bin") == 1);
    }

    // ------------------------------------------------------------------ //
    // header.txt blocked by a directory → Status::error,
    // the binary files are still written
    // ------------------------------------------------------------------ //

    {
        fs::remove_all(OUT_ROOT);
        const fs::path out = OUT_ROOT / "blocked_header";
        fs::create_directories(out / "binary" / "header.txt");

        CHECK(out_write(out) == Status::error);
        CHECK(out_size(out / "binary" / "sim_pos.bin")         == vec3_bytes);
        CHECK(out_size(out / "binary" / "sim_normal_form.bin") == energy_bytes);
    }

    // ------------------------------------------------------------------ //
    // Ovito directory blocked by a regular file → Status::error,
    // the binary files are still written, also those after the Ovito files
    // ------------------------------------------------------------------ //

    {
        fs::remove_all(OUT_ROOT);
        const fs::path out = OUT_ROOT / "blocked_ovito";
        out_write_text(out / "ovito", "x");

        CHECK(out_write(out) == Status::error);
        CHECK(out_size(out / "binary" / "sim_pos.bin")         == vec3_bytes);
        CHECK(out_size(out / "binary" / "sim_normal_form.bin") == energy_bytes);
    }

    // ------------------------------------------------------------------ //
    // One Ovito file blocked by a directory → Status::error,
    // the other snapshots are still written
    // ------------------------------------------------------------------ //

    {
        fs::remove_all(OUT_ROOT);
        const fs::path out = OUT_ROOT / "blocked_dump";
        fs::create_directories(out / "ovito" / "t_00000.dump");

        CHECK(out_write(out) == Status::error);
        CHECK(fs::is_regular_file(out / "ovito" / "t_00001.dump"));
    }

    fs::remove_all(OUT_ROOT);
}
