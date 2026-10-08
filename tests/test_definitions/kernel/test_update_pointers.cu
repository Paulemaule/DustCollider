/**
 * @file test_update_pointers.cu
 * @brief Tests for the inelastic sliding and rolling corrections of the contact pointers:
 * correct_contact_pointer() and the updatePointers kernel in physics/integrator.cuh.
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * Every test uses a single contact between two touching monomers (normal displacement 0, so the contact
 * does not break). The contact pointers are built for a prescribed sliding and rolling displacement, in
 * multiples of the critical values. After the correction:
 *  - a displacement above its critical value lies exactly on it, with its direction unchanged,
 *  - a displacement below its critical value is unchanged, also when the other one was corrected
 *    (for r_i != r_j the earlier pointer shifts coupled the two, see CLAUDE.md (AC)),
 *  - the pointers are unit vectors.
 * The displacements are measured with up_displacements, which is written out independently of the code
 * under test and follows the original kernel formulas.
 *
 * The host part calls correct_contact_pointer for both threads of the pair, (0,1) and (1,0), each with the
 * uncorrected partner pointer, as the kernel does. The kernel part runs updatePointers on a two monomer
 * system, which adds the pair matrix layout, the corotation of the pointers and the energy bookkeeping.
 *
 * Tolerance: measuring a displacement subtracts terms of order r from each other, which leaves a rounding
 * floor of ~1e-13 of the critical displacement. UP_TOL keeps a wide margin to that floor. The pointer shifts
 * before (AC) missed the critical value by ~1e-4 for equal radii and by O(1) for unequal radii.
 */

#include <cmath>
#include <cstdio>
#include "utils/constant.cuh"
#include "utils/vector.cuh"
#include "utils/buffer.cuh"
#include "physics/integrator_utils.cuh"
#include "physics/integrator.cuh"
#include "physics/state.cuh"
#include "physics/materials.cuh"

// Relative tolerance of the displacements, in units of the critical displacement.
static const double UP_TOL = 1e-10;

// Forsterite, as in the shipped command files.
static const double UP_GAMMA = 0.07;
static const double UP_E     = 204e9;
static const double UP_NU    = 0.24;
static const double UP_XI    = 2e-10;

/**
 * A contact between two touching monomers.
 */
struct UpContact {
    double  r0, r1;             // The monomer radii.
    double3 x0, x1;             // The monomer positions.
    double3 n0, n1;             // The contact pointers of monomer 0 and 1 in the lab frame.
    double  delta_S_crit;       // The critical sliding displacement of the pair.
    double  delta_R_crit;       // The critical rolling displacement of the pair.
};

/**
 * The prescribed displacements of a test case, in multiples of the critical values.
 */
struct UpCase {
    const char* label;
    double      f_S;            // The sliding displacement in units of delta_S_crit.
    double      f_R;            // The rolling displacement in units of delta_R_crit.
};

static const UpCase UP_CASES[] = {
    { "rolling only",  0.5,  1.5 },
    { "sliding only",  1.5,  0.5 },
    { "both",          2.0,  2.0 },
    { "both, 10x",    10.0, 10.0 },
    { "neither",       0.5,  0.5 },
};

// The radius pairs. Swapping the radii swaps the roles of the two threads, equal radii are the reference case.
static const double UP_RADII[][2] = {
    { 15e-9, 30e-9 },
    { 30e-9, 15e-9 },
    { 20e-9, 20e-9 },
};

// a * u + b * v
static double3 up_lin(const double a, const double3 u, const double b, const double3 v) {
    return { a * u.x + b * v.x, a * u.y + b * v.y, a * u.z + b * v.z };
}

static bool up_equal(const double3 u, const double3 v) {
    return u.x == v.x && u.y == v.y && u.z == v.z;
}

// CHECK with the test case printed above a failure.
static void up_check(const bool ok, const char* what, const UpContact& c, const char* where) {
    if (!ok) printf("    updatePointers [%s, r = %.0f/%.0f nm]: %s\n", where, c.r0 * 1e9, c.r1 * 1e9, what);
    CHECK(ok);
}

/**
 * Builds a touching contact whose sliding displacement is about f_S * delta_S_crit and whose rolling
 * displacement is about f_R * delta_R_crit. The two point in tangential directions 60 degrees apart, and
 * the contact axis is not aligned with a coordinate axis, so no component vanishes trivially.
 * The pointers are shifted from the reference configuration in the decoupled basis
 *      n0 = norm(-n + s / (r0 + r1) * t_S + rho / r0 * t_R)
 *      n1 = norm( n - s / (r0 + r1) * t_S + rho / r1 * t_R)
 * The normalization changes the displacements slightly, so the tests measure them instead of using f_S, f_R.
 */
static UpContact up_make_contact(const double r0, const double r1, const double f_S, const double f_R) {
    UpContact c;
    c.r0 = r0;
    c.r1 = r1;

    // The critical displacements, as the kernel computes them.
    double R   = get_R(r0, r1);
    double a_0 = get_a_0(get_gamma(UP_GAMMA, UP_GAMMA), R, get_E_s(UP_E, UP_E, UP_NU, UP_NU));
    c.delta_S_crit = get_delta_S_crit(UP_NU, UP_NU, a_0);
    c.delta_R_crit = 0.5 * (UP_XI + UP_XI);

    double3 axis = { 1.0, 0.3, -0.2 };      // The direction from monomer 0 to monomer 1.
    vec_normalize(axis);
    double3 t_S = vec_cross(axis, { 0.0, 0.0, 1.0 });
    vec_normalize(t_S);
    double3 t_R = up_lin(0.5, t_S, 0.5 * sqrt(3.0), vec_cross(axis, t_S));

    c.x0 = { 1.0e-7, -2.0e-8, 5.0e-8 };
    c.x1 = up_lin(1.0, c.x0, r0 + r1, axis);

    // The contact normal as seen by the (0,1) thread, pointing from monomer 1 to monomer 0.
    double3 n = vec_get_normal(c.x0, c.x1);

    double s   = f_S * c.delta_S_crit;
    double rho = f_R * c.delta_R_crit;

    c.n0 = up_lin(-1.0, n, 1.0, up_lin( s / (r0 + r1), t_S, rho / r0, t_R));
    c.n1 = up_lin( 1.0, n, 1.0, up_lin(-s / (r0 + r1), t_S, rho / r1, t_R));
    vec_normalize(c.n0);
    vec_normalize(c.n1);

    return c;
}

/**
 * The sliding displacement and the tangential part of the rolling displacement of the contact for the
 * pointers n0, n1 (Wada et al. 2007, as in the kernels): the projections of
 *      zeta_0 = r0 * n0 - r1 * n1 + (r0 + r1) * n      and      R * (n0 + n1)
 * onto the contact plane.
 */
static void up_displacements(const UpContact& c, const double3 n0, const double3 n1, double3& sliding, double3& rolling) {
    double3 n = vec_get_normal(c.x0, c.x1);
    double  R = c.r0 * c.r1 / (c.r0 + c.r1);

    double3 zeta_0 = up_lin(1.0, up_lin(c.r0, n0, -c.r1, n1), c.r0 + c.r1, n);
    double3 xi     = up_lin(R, n0, R, n1);

    sliding = up_lin(1.0, zeta_0, -vec_dot(zeta_0, n), n);
    rolling = up_lin(1.0, xi,     -vec_dot(xi, n),     n);
}

// The displacement after a correction: scaled onto the critical value if it exceeds it, unchanged otherwise.
static double3 up_expected(const double3 d, const double crit) {
    double length = vec_length(d);
    return length > crit ? up_lin(crit / length, d, 0.0, d) : d;
}

// The length by which a displacement exceeds its critical value, 0 if it does not.
static double up_excess(const double3 d, const double crit) {
    return fmax(0.0, vec_length(d) - crit);
}

/**
 * Checks the corrected lab frame pointers n0_new, n1_new against the displacements before the correction.
 */
static void up_check_corrected(const UpContact& c, const double3 n0_new, const double3 n1_new, const char* where) {
    double3 sliding, rolling, sliding_new, rolling_new;
    up_displacements(c, c.n0,   c.n1,   sliding,     rolling);
    up_displacements(c, n0_new, n1_new, sliding_new, rolling_new);

    double3 sliding_exp = up_expected(sliding, c.delta_S_crit);
    double3 rolling_exp = up_expected(rolling, c.delta_R_crit);

    up_check(vec_length(vec_diff(sliding_new, sliding_exp)) / c.delta_S_crit < UP_TOL, "sliding displacement after the correction", c, where);
    up_check(vec_length(vec_diff(rolling_new, rolling_exp)) / c.delta_R_crit < UP_TOL, "rolling displacement after the correction", c, where);
    up_check(fabs(vec_length(n0_new) - 1.0) < 1e-14, "pointer 0 is a unit vector", c, where);
    up_check(fabs(vec_length(n1_new) - 1.0) < 1e-14, "pointer 1 is a unit vector", c, where);
}

/**
 * Runs correct_contact_pointer for both threads of the pair on the host.
 */
static void up_test_host(const UpContact& c, const char* label) {
    char where[96];
    snprintf(where, sizeof(where), "%s, host", label);

    double3 n0_new = c.n0;
    double3 n1_new = c.n1;
    double  excess_S0, excess_R0, excess_S1, excess_R1;

    // Each thread sees the contact normal pointing towards its own monomer and the uncorrected partner pointer.
    bool corrected_0 = correct_contact_pointer(n0_new, c.n1, vec_get_normal(c.x0, c.x1), c.r0, c.r1, c.delta_S_crit, c.delta_R_crit, excess_S0, excess_R0);
    bool corrected_1 = correct_contact_pointer(n1_new, c.n0, vec_get_normal(c.x1, c.x0), c.r1, c.r0, c.delta_S_crit, c.delta_R_crit, excess_S1, excess_R1);

    double3 sliding, rolling;
    up_displacements(c, c.n0, c.n1, sliding, rolling);
    double excess_S = up_excess(sliding, c.delta_S_crit);
    double excess_R = up_excess(rolling, c.delta_R_crit);
    bool   expect_correction = excess_S > 0. || excess_R > 0.;

    up_check(corrected_0 == expect_correction && corrected_1 == expect_correction, "correction flag", c, where);
    up_check(fabs(excess_S0 - excess_S) / c.delta_S_crit < UP_TOL && fabs(excess_S1 - excess_S) / c.delta_S_crit < UP_TOL, "sliding excess", c, where);
    up_check(fabs(excess_R0 - excess_R) / c.delta_R_crit < UP_TOL && fabs(excess_R1 - excess_R) / c.delta_R_crit < UP_TOL, "rolling excess", c, where);

    if (!expect_correction) {
        up_check(up_equal(n0_new, c.n0) && up_equal(n1_new, c.n1), "pointers are untouched below both critical values", c, where);
    }

    up_check_corrected(c, n0_new, n1_new, where);
}

/**
 * Runs the updatePointers kernel on a system of the two monomers of the contact. Both stored pointers are
 * corotated with the quaternion q.
 */
static void up_test_kernel(const UpContact& c, const double4 q, const char* label) {
    char where[96];
    snprintf(where, sizeof(where), "%s, kernel%s", label, q.w == 1. ? "" : ", rotated");

    const int Nmon = 2;

    // Pair matrix layout: entry i + j * Nmon belongs to thread (i, j) and holds the pointer of monomer i.
    const int slot_01 = 0 + 1 * Nmon;
    const int slot_10 = 1 + 0 * Nmon;

    HostMaterials     host_materials(Nmon);
    HostMaterialsView hm = host_materials.view();
    for (int k = 0; k < Nmon; k++) {
        hm.radius[k]            = k == 0 ? c.r0 : c.r1;
        hm.mass[k]              = 1.;
        hm.moment[k]            = 1.;
        hm.matID[k]             = 0;
        hm.density[k]           = 3210.;
        hm.surface_energy[k]    = UP_GAMMA;
        hm.youngs_modulus[k]    = UP_E;
        hm.poisson_number[k]    = UP_NU;
        hm.damping_timescale[k] = 1e-12;
        hm.crit_rolling_disp[k] = UP_XI;
    }

    HostState     host_state(Nmon);
    HostStateView hs = host_state.view();
    for (int k = 0; k < Nmon; k++) {
        hs.velocity[k] = { 0., 0., 0. };
        hs.omega[k]    = { 0., 0., 0. };
        hs.force[k]    = { 0., 0., 0. };
        hs.torque[k]   = { 0., 0., 0. };
    }
    for (int k = 0; k < Nmon * Nmon; k++) {
        hs.contact_compression[k] = 0.;
        hs.contact_twist[k]       = 0.;
        hs.contact_pointer[k]     = { 0., 0., 0. };
        hs.contact_normal[k]      = { 0., 0., 0. };
        hs.contact_rotation[k]    = { 0., 0., 0., 0. };
    }
    hs.position[0] = c.x0;
    hs.position[1] = c.x1;

    // The kernel recovers the lab frame pointers with quat_apply_inverse(rotation, pointer).
    hs.contact_pointer[slot_01]  = quat_apply(q, c.n0);
    hs.contact_pointer[slot_10]  = quat_apply(q, c.n1);
    hs.contact_rotation[slot_01] = q;
    hs.contact_rotation[slot_10] = q;

    DeviceMaterials device_materials(Nmon);
    host_materials.push_to(device_materials, Nmon);

    // updatePointers reads the twisting displacement from the 'next' buffer, so both buffers get the same state.
    DeviceState curr(Nmon);
    DeviceState next(Nmon);
    host_state.push_to(curr, Nmon);
    host_state.push_to(next, Nmon);

    DeviceBuffer<double4> inelastic(1);
    CHECK_CUDA(cudaMemset(inelastic.data(), 0, sizeof(double4)));

    DeviceStateView     cv = curr.view();
    DeviceStateView     nv = next.view();
    DeviceMaterialsView mv = device_materials.view();

    updatePointers<<<1, Nmon * Nmon>>>(
        cv.position, cv.contact_pointer,
        cv.contact_rotation, cv.contact_compression,
        nv.contact_pointer, nv.contact_rotation,
        nv.contact_twist, nv.contact_compression,
        inelastic.data(),
        mv.radius, mv.youngs_modulus, mv.poisson_number,
        mv.surface_energy, mv.crit_rolling_disp,
        Nmon
    );
    cudaError_t launch_error = cudaGetLastError();
    cudaError_t sync_error   = cudaDeviceSynchronize();
    up_check(launch_error == cudaSuccess && sync_error == cudaSuccess, "kernel execution", c, where);

    HostState     out_state(Nmon);
    out_state.pull_from(next, Nmon);
    HostStateView out = out_state.view();

    double4 booked;
    CHECK_CUDA(cudaMemcpy(&booked, inelastic.data(), sizeof(double4), cudaMemcpyDeviceToHost));

    // Pointers.
    up_check(vec_length_sq(out.contact_pointer[0]) == 0. && vec_length_sq(out.contact_pointer[3]) == 0., "the diagonal of the pair matrix is untouched", c, where);

    double3 sliding, rolling;
    up_displacements(c, c.n0, c.n1, sliding, rolling);
    double excess_S = up_excess(sliding, c.delta_S_crit);
    double excess_R = up_excess(rolling, c.delta_R_crit);

    if (excess_S == 0. && excess_R == 0.) {
        up_check(up_equal(out.contact_pointer[slot_01], hs.contact_pointer[slot_01]) && up_equal(out.contact_pointer[slot_10], hs.contact_pointer[slot_10]),
                 "pointers are copied unchanged below both critical values", c, where);
    }

    up_check_corrected(c, quat_apply_inverse(q, out.contact_pointer[slot_01]), quat_apply_inverse(q, out.contact_pointer[slot_10]), where);

    // Energy bookkeeping: k * delta_crit * excess for the pair, half of it from each thread.
    double R     = get_R(c.r0, c.r1);
    double G     = get_G_i(UP_E, UP_NU);
    double gamma = get_gamma(UP_GAMMA, UP_GAMMA);
    double a_0   = get_a_0(gamma, R, get_E_s(UP_E, UP_E, UP_NU, UP_NU));
    double k_s   = 8. * get_G_s(G, G, UP_NU, UP_NU) * a_0;
    double k_r   = 4. * (3. * PI * gamma * R) / R;

    double dissipated_S = k_s * c.delta_S_crit * excess_S;
    double dissipated_R = k_r * c.delta_R_crit * excess_R;

    if (dissipated_S > 0.) up_check(fabs(booked.x / dissipated_S - 1.) < 1e-9, "booked sliding dissipation", c, where);
    else                   up_check(booked.x == 0., "no sliding dissipation below the critical value", c, where);
    if (dissipated_R > 0.) up_check(fabs(booked.y / dissipated_R - 1.) < 1e-9, "booked rolling dissipation", c, where);
    else                   up_check(booked.y == 0., "no rolling dissipation below the critical value", c, where);
    up_check(booked.z == 0. && booked.w == 0., "no twisting or normal dissipation", c, where);
}

void test_update_pointers() {
    const double4 identity   = { 0., 0., 0., 1. };
    // A rotation by 0.7 rad around the axis (1, 2, 2) / 3.
    const double  half_angle = 0.35;
    const double4 rotation   = { sin(half_angle) / 3., 2. * sin(half_angle) / 3., 2. * sin(half_angle) / 3., cos(half_angle) };

    for (const auto& radii : UP_RADII) {
        for (const UpCase& test_case : UP_CASES) {
            UpContact c = up_make_contact(radii[0], radii[1], test_case.f_S, test_case.f_R);

            up_test_host(c, test_case.label);
            up_test_kernel(c, identity, test_case.label);
            up_test_kernel(c, rotation, test_case.label);
        }
    }
}
