/**
 * @file test_aggregate.cu
 * @brief Tests for Aggregate::from_file() in simulationSetup/aggregate.cuh
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * Aggregate file format (positions and radius in nm, mat_id 1-indexed):
 *   Line 0:   Nmon  external_radius_nm  effective_radius_nm
 *   Lines 1-4: metadata — skipped
 *   Lines 5+:  x  y  z  _  radius  _  mat_id   (7 values per line)
 *
 * Error paths throw std::runtime_error (not Status).
 */

#include "simulationSetup/aggregate.cuh"
#include <stdexcept>
#include <cstdio>

static const char* AGG_TMP = "/tmp/test_dustcollider_aggregate.tmp";

static void write_agg_tmp(const char* content) {
    FILE* f = fopen(AGG_TMP, "w");
    if (f) { fputs(content, f); fclose(f); }
}

void test_aggregate() {

    // ------------------------------------------------------------------ //
    // File not found → throws std::runtime_error
    // ------------------------------------------------------------------ //

    {
        bool threw = false;
        try {
            Aggregate::from_file("/tmp/does_not_exist_dustcollider_agg.tmp");
        } catch (const std::runtime_error&) { threw = true; }
        CHECK(threw);
    }

    // ------------------------------------------------------------------ //
    // Valid 2-monomer aggregate: header, 4 skipped lines, 2 monomer lines
    // ------------------------------------------------------------------ //

    write_agg_tmp(
        "2 100.0 50.0\n"                         // line 0: header
        "# metadata 1\n"                          // line 1: skipped
        "# metadata 2\n"                          // line 2: skipped
        "# metadata 3\n"                          // line 3: skipped
        "# metadata 4\n"                          // line 4: skipped
        "10.0 20.0 30.0 0.0 50.0 0.0 1\n"        // monomer 0: mat_id 1 → 0
        "0.0 0.0 -100.0 0.0 75.0 0.0 2\n"        // monomer 1: mat_id 2 → 1
    );

    {
        Aggregate agg = Aggregate::from_file(AGG_TMP);

        // Header
        CHECK(agg.header.Nmon == 2);
        CHECK_APPROX(agg.header.external_radius,  100e-9, 1e-12);
        CHECK_APPROX(agg.header.effective_radius,  50e-9, 1e-12);

        // Monomer count
        CHECK(agg.monomers.positions.size()    == 2);
        CHECK(agg.monomers.radii.size()        == 2);
        CHECK(agg.monomers.material_ids.size() == 2);

        // Monomer 0: position (nm → m), radius (nm → m), mat_id (1 → 0)
        CHECK_APPROX(agg.monomers.positions[0].x,  10e-9, 1e-12);
        CHECK_APPROX(agg.monomers.positions[0].y,  20e-9, 1e-12);
        CHECK_APPROX(agg.monomers.positions[0].z,  30e-9, 1e-12);
        CHECK_APPROX(agg.monomers.radii[0],        50e-9, 1e-12);
        CHECK(agg.monomers.material_ids[0] == 0);

        // Monomer 1
        CHECK_APPROX(agg.monomers.positions[1].x,   0.0,    1e-12);
        CHECK_APPROX(agg.monomers.positions[1].z, -100e-9,  1e-12);
        CHECK_APPROX(agg.monomers.radii[1],         75e-9,  1e-12);
        CHECK(agg.monomers.material_ids[1] == 1);
    }

    // ------------------------------------------------------------------ //
    // Wrong number of values on a monomer line → throws
    // ------------------------------------------------------------------ //

    {
        write_agg_tmp(
            "1 100.0 50.0\n"
            "skip\nskip\nskip\nskip\n"
            "10.0 20.0 30.0\n"           // only 3 values instead of 7
        );
        bool threw = false;
        try {
            Aggregate::from_file(AGG_TMP);
        } catch (const std::runtime_error&) { threw = true; }
        CHECK(threw);
    }

    // ------------------------------------------------------------------ //
    // Cleanup
    // ------------------------------------------------------------------ //

    remove(AGG_TMP);
}
