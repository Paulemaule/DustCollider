/**
 * @file test_cluster.cu
 * @brief Tests for findMonomerClusters(), contactMatrixBytes() and packContactMatrix() in
 * physics/integrator_utils.cuh
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

// Bit (a, b) of a packed contact matrix, unpacked in the bit order of numpy.unpackbits.
static bool contact_bit(const uint8_t* packed, int Nmon, int a, int b) {
    const size_t k = size_t(a) * Nmon + b;
    return (packed[k / 8] >> (7 - k % 8)) & 1u;
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

    // ------------------------------------------------------------------ //
    // contactMatrixBytes: one bit per pair, rounded up to whole bytes
    // ------------------------------------------------------------------ //

    {
        CHECK(contactMatrixBytes(0) == 0);
        CHECK(contactMatrixBytes(1) == 1);
        CHECK(contactMatrixBytes(2) == 1);      // 4 bits
        CHECK(contactMatrixBytes(3) == 2);      // 9 bits
        CHECK(contactMatrixBytes(4) == 2);      // 16 bits
        CHECK(contactMatrixBytes(2144) == 574592);
    }

    // ------------------------------------------------------------------ //
    // Bit order: entry k = 1 alone sets the second most significant bit
    // of byte 0, as numpy.packbits does
    // ------------------------------------------------------------------ //

    {
        const int N = 4;
        double3 cp[16] = {};
        cp[0 * N + 1] = {0.0, 0.0, -1.0};
        uint8_t packed[2] = {0xFF, 0xFF};
        packContactMatrix(N, cp, packed);
        CHECK(packed[0] == 0x40);
        CHECK(packed[1] == 0x00);
    }

    // ------------------------------------------------------------------ //
    // Chain 0-1-2 of N = 3: bits 1, 3, 5, 7 in byte 0; byte 1 holds only
    // the diagonal entry (2, 2) and the padding, which are cleared
    // ------------------------------------------------------------------ //

    {
        const int N = 3;
        double3 cp[9] = {};
        connect(cp, N, 0, 1);
        connect(cp, N, 1, 2);
        uint8_t packed[2] = {0xFF, 0xFF};
        packContactMatrix(N, cp, packed);
        CHECK(packed[0] == 0x55);
        CHECK(packed[1] == 0x00);
    }

    // ------------------------------------------------------------------ //
    // N = 4 with contacts 0-2, 1-3 and 2-3: every bit matches the pointer
    // matrix, the matrix is symmetric and its diagonal is 0
    // ------------------------------------------------------------------ //

    {
        const int N = 4;
        double3 cp[16] = {};
        connect(cp, N, 0, 2);
        connect(cp, N, 1, 3);
        connect(cp, N, 2, 3);
        uint8_t packed[2] = {};
        packContactMatrix(N, cp, packed);

        int n_set = 0;
        for (int a = 0; a < N; a++) {
            for (int b = 0; b < N; b++) {
                CHECK(contact_bit(packed, N, a, b) == (vec_length_sq(cp[a * N + b]) != 0.0));
                CHECK(contact_bit(packed, N, a, b) == contact_bit(packed, N, b, a));
                if (contact_bit(packed, N, a, b)) n_set++;
            }
            CHECK(!contact_bit(packed, N, a, a));
        }
        CHECK(n_set == 6);
    }
}
