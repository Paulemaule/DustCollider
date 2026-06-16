/**
 * @file test_cluster.cu
 * @brief Tests for findMonomerClusters() in physics/integrator_utils.cuh
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * The contact_pointer matrix is Nmon×Nmon; entry [i*Nmon+j] is non-zero
 * (any component != 0) iff monomers i and j are in contact.
 */

#include "utils/constant.cuh"
#include "utils/vector.cuh"
#include "physics/integrator_utils.cuh"

static void connect(double3* mat, int Nmon, int i, int j) {
    mat[i * Nmon + j] = {1.0, 0.0, 0.0};
    mat[j * Nmon + i] = {1.0, 0.0, 0.0};
}

void test_cluster() {

    // ------------------------------------------------------------------ //
    // Single isolated monomer → assigned cluster id 0
    // ------------------------------------------------------------------ //

    {
        const int N = 1;
        double3 cp[1] = {};
        int cl[1] = {-1};
        findMonomerClusters(N, cp, cl);
        CHECK(cl[0] == 0);
    }

    // ------------------------------------------------------------------ //
    // Two isolated monomers → distinct ids assigned in discovery order
    // ------------------------------------------------------------------ //

    {
        const int N = 2;
        double3 cp[4] = {};
        int cl[2] = {-1, -1};
        findMonomerClusters(N, cp, cl);
        CHECK(cl[0] == 0);
        CHECK(cl[1] == 1);
        CHECK(cl[0] != cl[1]);
    }

    // ------------------------------------------------------------------ //
    // Two monomers in contact → same cluster id
    // ------------------------------------------------------------------ //

    {
        const int N = 2;
        double3 cp[4] = {};
        connect(cp, N, 0, 1);
        int cl[2] = {-1, -1};
        findMonomerClusters(N, cp, cl);
        CHECK(cl[0] == cl[1]);
    }

    // ------------------------------------------------------------------ //
    // Linear chain 0–1–2 → all three in the same cluster
    // ------------------------------------------------------------------ //

    {
        const int N = 3;
        double3 cp[9] = {};
        connect(cp, N, 0, 1);
        connect(cp, N, 1, 2);
        int cl[3] = {-1, -1, -1};
        findMonomerClusters(N, cp, cl);
        CHECK(cl[0] == cl[1]);
        CHECK(cl[1] == cl[2]);
    }

    // ------------------------------------------------------------------ //
    // Two separate pairs (0–1) and (2–3) → exactly two clusters
    // ------------------------------------------------------------------ //

    {
        const int N = 4;
        double3 cp[16] = {};
        connect(cp, N, 0, 1);
        connect(cp, N, 2, 3);
        int cl[4] = {-1, -1, -1, -1};
        findMonomerClusters(N, cp, cl);
        CHECK(cl[0] == cl[1]);
        CHECK(cl[2] == cl[3]);
        CHECK(cl[0] != cl[2]);
        CHECK(cl[0] == 0);
        CHECK(cl[2] == 1);
    }

    // ------------------------------------------------------------------ //
    // Star topology: center 0 connected to leaves 1, 2, 3 → one cluster
    // ------------------------------------------------------------------ //

    {
        const int N = 4;
        double3 cp[16] = {};
        for (int leaf = 1; leaf < N; leaf++)
            connect(cp, N, 0, leaf);
        int cl[4] = {-1, -1, -1, -1};
        findMonomerClusters(N, cp, cl);
        CHECK(cl[0] == cl[1]);
        CHECK(cl[1] == cl[2]);
        CHECK(cl[2] == cl[3]);
    }
}
