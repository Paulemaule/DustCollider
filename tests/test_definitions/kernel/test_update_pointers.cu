/**
 * @file test_update_pointers.cu
 * @brief Tests for the pointer correction of inelastic sliding and rolling in updatePointers.
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * Every case builds a contact with prescribed sliding and rolling displacements (see two_monomer.cuh), runs
 * updatePointers and measures the displacements of the corrected pointers with tm_displacements. A displacement
 * beyond its critical value must end up on it, in its old direction, and the other displacement must stay unchanged.
 * The rolling displacement is measured in the contact plane, where it is limited.
 *
 * Every case runs for the radii 15/30, 30/15 and 20/20 nm, with the contact pointers stored unrotated and rotated by
 * 0.7 rad. The energy booking of the correction is tested in test_energy.cu.
 */

#include <cmath>
#include <cstdio>
#include "test_definitions/kernel/two_monomer.cuh"

// The tolerance of the displacements, relative to their critical values. Measuring a displacement subtracts terms of
// the size of the radii, which leaves a rounding floor of ~1e-13.
static const double UP_TOL = 1e-10;

/**
 * The displacement d limited to the length delta_crit, in its own direction.
 * Written out independently of the kernel helpers.
 */
static double3 up_limited(const double3 d, const double delta_crit) {
    double length = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    if (length <= delta_crit) return d;
    return { d.x * delta_crit / length, d.y * delta_crit / length, d.z * delta_crit / length };
}

/**
 * A copy of the system s with the contact pointers of the 'next' state of dev, derotated with q.
 */
static TmSystem up_corrected_system(TmDevice& dev, const TmSystem& s, const double4 q) {
    HostState out(TM_NMON);
    out.pull_from(dev.next, TM_NMON);

    TmSystem c = s;
    c.n[0] = quat_apply_inverse(q, out.view().contact_pointer[TM_SLOT_01]);
    c.n[1] = quat_apply_inverse(q, out.view().contact_pointer[TM_SLOT_10]);
    return c;
}

/**
 * Checks that a measured displacement agrees with its expected value within UP_TOL * delta_crit.
 */
static void up_check_displacement(const double3 measured, const double3 expected, const double delta_crit,
                                  const char* label, const char* dof) {
    double miss = vec_length(vec_diff(measured, expected)) / delta_crit;
    bool ok = miss <= UP_TOL;
    if (!ok)
        printf("    %s: the %s displacement misses its expected value by %.3g of its critical value\n", label, dof, miss);
    CHECK(ok);
}

static bool up_identical(const double3 u, const double3 v) {
    return u.x == v.x && u.y == v.y && u.z == v.z;
}

/**
 * A contact with the sliding and rolling displacements f_S and f_R times their critical values.
 */
static void up_test_correction(const double r0, const double r1, const double f_S, const double f_R, const double4 q) {
    char label[96];
    snprintf(label, sizeof(label), "updatePointers, %.0f/%.0f nm, sliding %.1f, rolling %.1f", r0 * 1e9, r1 * 1e9, f_S, f_R);

    TmPair   p = tm_pair(r0, r1);
    TmSystem s = tm_make_system(r0, r1, p.delta_N_0, true, f_S, f_R, 0.5);

    double normal; double3 sliding, rolling;
    tm_displacements(s, normal, sliding, rolling);

    TmDevice dev(s, q);
    dev.run_update_pointers();
    dev.check_launch(label);
    TmSystem c = up_corrected_system(dev, s, q);

    // The displacements of the corrected pointers: each one is limited to its critical value, independently of the
    // other one.
    double normal_c; double3 sliding_c, rolling_c;
    tm_displacements(c, normal_c, sliding_c, rolling_c);
    up_check_displacement(sliding_c, up_limited(sliding, p.delta_S_crit), p.delta_S_crit, label, "sliding");
    up_check_displacement(tm_tangential(c, rolling_c), up_limited(tm_tangential(s, rolling), p.delta_R_crit),
                          p.delta_R_crit, label, "rolling");

    // The pointers are unit vectors, and each points towards the other monomer.
    double3 axis = vec_diff(s.x[1], s.x[0]);
    CHECK(fabs(vec_length(c.n[0]) - 1.) < 1e-14 && fabs(vec_length(c.n[1]) - 1.) < 1e-14);
    CHECK(vec_dot(c.n[0], axis) > 0. && vec_dot(c.n[1], axis) < 0.);

    // Below both critical values the stored pointers are copied unchanged.
    if (f_S < 1. && f_R < 1.) {
        HostState out_curr(TM_NMON), out_next(TM_NMON);
        out_curr.pull_from(dev.curr, TM_NMON);
        out_next.pull_from(dev.next, TM_NMON);
        CHECK(up_identical(out_next.view().contact_pointer[TM_SLOT_01], out_curr.view().contact_pointer[TM_SLOT_01]));
        CHECK(up_identical(out_next.view().contact_pointer[TM_SLOT_10], out_curr.view().contact_pointer[TM_SLOT_10]));
    }

    // A second step from the corrected pointers finds no inelastic motion: it books nothing and keeps the
    // displacements.
    char label_2[112];
    snprintf(label_2, sizeof(label_2), "%s, second step", label);

    TmDevice dev_2(c, q);
    dev_2.run_update_pointers();
    dev_2.check_energy(EnergyRecord{}, p.F_c * p.delta_N_crit, label_2);
    TmSystem c_2 = up_corrected_system(dev_2, c, q);

    double normal_2; double3 sliding_2, rolling_2;
    tm_displacements(c_2, normal_2, sliding_2, rolling_2);
    up_check_displacement(sliding_2, sliding_c, p.delta_S_crit, label_2, "sliding");
    up_check_displacement(tm_tangential(c_2, rolling_2), tm_tangential(c, rolling_c), p.delta_R_crit, label_2, "rolling");
}

static void up_test_radii(const double r0, const double r1, const double4 q) {
    up_test_correction(r0, r1, 0.5,  1.5,  q);    // Rolling only.
    up_test_correction(r0, r1, 1.5,  0.5,  q);    // Sliding only.
    up_test_correction(r0, r1, 2.,   2.,   q);    // Both.
    up_test_correction(r0, r1, 10.,  10.,  q);    // Both, far beyond the critical values.
    up_test_correction(r0, r1, 0.5,  0.5,  q);    // Neither.
}

void test_update_pointers() {
    const double4 rotations[2] = { TM_IDENTITY, TM_ROTATION };
    for (const double4 q : rotations) {
        up_test_radii(15e-9, 30e-9, q);
        up_test_radii(30e-9, 15e-9, q);
        up_test_radii(20e-9, 20e-9, q);
    }
}
