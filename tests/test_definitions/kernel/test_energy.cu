/**
 * @file test_energy.cu
 * @brief Tests for the energy trackers (physics/energy.cuh): which kernel books which energy into which tracker.
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * Every test runs a kernel on a system of two monomers with one contact, or none. Each thread of the pair books
 * half of the energy of the pair into the slot of its own monomer. So both slots of a tracker must hold half of
 * the expected pair energy, and all other trackers must stay empty.
 * The expected energies use the potentials of integrator_utils.cuh. The displacements they are evaluated at come
 * from en_displacements, which is written out independently of the kernels.
 *
 * Every case runs with the contact pointers stored unrotated and rotated by 0.7 rad, the kernels have to undo
 * the rotation.
 *
 * The potentials are not booked by evaluate but by contact_potentials, which runs on the stored state at every
 * snapshot.
 */

#include <cmath>
#include <cstdio>
#include "utils/constant.cuh"
#include "utils/vector.cuh"
#include "utils/buffer.cuh"
#include "physics/integrator_utils.cuh"
#include "physics/integrator.cuh"
#include "physics/energy.cuh"
#include "physics/state.cuh"
#include "physics/materials.cuh"

// Relative tolerance of the booked energies.
static const double EN_TOL = 1e-9;

// Forsterite, as in the shipped command files.
static const double EN_GAMMA = 0.07;
static const double EN_E     = 204e9;
static const double EN_NU    = 0.24;
static const double EN_XI    = 2e-10;
static const double EN_TVIS  = 1e-12;

// A timestep of the order the auto-timestep picks for these monomers.
static const double EN_DT    = 2e-14;

// The test systems consist of two monomers. Thread (0,1) works on the pair matrix slot 0 + 1 * Nmon,
// thread (1,0) on the slot 1 + 0 * Nmon.
static const int EN_NMON    = 2;
static const int EN_SLOT_01 = 2;
static const int EN_SLOT_10 = 1;

/**
 * The constants of a pair of forsterite monomers.
 */
struct EnPair {
    double R;               // The reduced radius.
    double E_s;             // The combined Youngs modulus.
    double a_0;             // The equilibrium contact radius.
    double F_c;             // The critical force.
    double k_s, k_r, k_t;   // The sliding, rolling and twisting stiffness.
    double delta_N_0;       // The equilibrium normal displacement.
    double delta_N_crit;    // The critical normal displacement.
    double delta_S_crit;    // The critical sliding displacement.
    double delta_R_crit;    // The critical rolling displacement.
    double delta_T_crit;    // The critical twisting displacement.
};

static EnPair en_pair(const double r0, const double r1) {
    EnPair p;
    double G     = get_G_i(EN_E, EN_NU);
    double gamma = get_gamma(EN_GAMMA, EN_GAMMA);

    p.R            = get_R(r0, r1);
    p.E_s          = get_E_s(EN_E, EN_E, EN_NU, EN_NU);
    p.a_0          = get_a_0(gamma, p.R, p.E_s);
    p.F_c          = 3. * PI * gamma * p.R;
    p.k_s          = 8. * get_G_s(G, G, EN_NU, EN_NU) * p.a_0;
    p.k_r          = 4. * p.F_c / p.R;
    p.k_t          = 16. * (G * G / (G + G)) * p.a_0 * p.a_0 * p.a_0 / 3.;
    p.delta_N_0    = p.a_0 * p.a_0 / (3. * p.R);
    p.delta_N_crit = get_delta_N_crit(p.a_0, p.R);
    p.delta_S_crit = get_delta_S_crit(EN_NU, EN_NU, p.a_0);
    p.delta_R_crit = EN_XI;
    p.delta_T_crit = 1. / (16. * PI);
    return p;
}

/**
 * A system of two monomers, with or without a contact.
 */
struct EnSystem {
    double  r[2];               // The monomer radii.
    double3 x[2];               // The monomer positions.
    double3 n[2];               // The contact pointers of monomer 0 and 1 in the lab frame, zero without a contact.
    double  twist;              // The twisting displacement of the contact.
    double  compression_old;    // The normal displacement stored in the contact (evaluate's previous value).
};

/**
 * Builds two monomers whose normal displacement is delta, along an axis that is not aligned with the coordinate
 * axes. With contact = true the pointers give sliding, rolling and twisting displacements of f_S, f_R and f_T times
 * their critical values; the sliding and rolling displacements are 60 degrees apart.
 */
static EnSystem en_make_system(const double r0, const double r1, const double delta, const bool contact,
                               const double f_S, const double f_R, const double f_T) {
    EnPair p = en_pair(r0, r1);

    // The axis from monomer 0 to monomer 1 and two tangential directions.
    double3 axis = { 1., 0.3, -0.2 };
    vec_normalize(axis);
    double3 e1 = vec_cross(axis, { 0., 0., 1. });
    vec_normalize(e1);
    double3 e2 = vec_cross(axis, e1);

    EnSystem s;
    s.r[0] = r0;
    s.r[1] = r1;

    double d = r0 + r1 - delta;
    s.x[0] = { 1e-7, -2e-8, 5e-8 };
    s.x[1] = { s.x[0].x + d * axis.x, s.x[0].y + d * axis.y, s.x[0].z + d * axis.z };

    s.twist           = contact ? f_T * p.delta_T_crit : 0.;
    s.compression_old = delta;

    if (!contact) {
        s.n[0] = { 0., 0., 0. };
        s.n[1] = { 0., 0., 0. };
        return s;
    }

    // The tangential displacements.
    double zeta = f_S * p.delta_S_crit;
    double xi   = f_R * p.delta_R_crit;
    double3 zeta_v = { zeta * e1.x, zeta * e1.y, zeta * e1.z };
    double3 xi_v   = { xi * (0.5 * e1.x + 0.5 * sqrt(3.) * e2.x),
                       xi * (0.5 * e1.y + 0.5 * sqrt(3.) * e2.y),
                       xi * (0.5 * e1.z + 0.5 * sqrt(3.) * e2.z) };

    // The tangential parts of the pointers that produce them: sliding = r0 t0 - r1 t1, rolling = R (t0 + t1).
    double3 t0 = { xi_v.x / r0 + zeta_v.x / (r0 + r1), xi_v.y / r0 + zeta_v.y / (r0 + r1), xi_v.z / r0 + zeta_v.z / (r0 + r1) };
    double3 t1 = { xi_v.x / r1 - zeta_v.x / (r0 + r1), xi_v.y / r1 - zeta_v.y / (r0 + r1), xi_v.z / r1 - zeta_v.z / (r0 + r1) };

    // Pointer 0 points towards monomer 1, pointer 1 towards monomer 0.
    double c0 = sqrt(1. - vec_length_sq(t0));
    double c1 = sqrt(1. - vec_length_sq(t1));
    s.n[0] = {  c0 * axis.x + t0.x,  c0 * axis.y + t0.y,  c0 * axis.z + t0.z };
    s.n[1] = { -c1 * axis.x + t1.x, -c1 * axis.y + t1.y, -c1 * axis.z + t1.z };
    return s;
}

/**
 * The displacements of the contact as seen by thread (0,1), following Wada et al. (2007).
 * Written out independently of the kernels.
 */
static void en_displacements(const EnSystem& s, double& normal, double3& sliding, double3& rolling) {
    double R = get_R(s.r[0], s.r[1]);

    double3 diff = vec_diff(s.x[0], s.x[1]);
    double  dist = vec_length(diff);
    double3 n    = { diff.x / dist, diff.y / dist, diff.z / dist };

    normal = s.r[0] + s.r[1] - dist;

    double3 zeta_0 = { s.r[0] * s.n[0].x - s.r[1] * s.n[1].x + (s.r[0] + s.r[1]) * n.x,
                       s.r[0] * s.n[0].y - s.r[1] * s.n[1].y + (s.r[0] + s.r[1]) * n.y,
                       s.r[0] * s.n[0].z - s.r[1] * s.n[1].z + (s.r[0] + s.r[1]) * n.z };
    double zeta_n = vec_dot(zeta_0, n);
    sliding = { zeta_0.x - zeta_n * n.x, zeta_0.y - zeta_n * n.y, zeta_0.z - zeta_n * n.z };
    rolling = { R * (s.n[0].x + s.n[1].x), R * (s.n[0].y + s.n[1].y), R * (s.n[0].z + s.n[1].z) };
}

/**
 * Checks that both monomer slots of an energy tracker hold the expected energy per monomer.
 * The comparison has an absolute floor of 1e-12 * energy_scale, for the trackers that are expected to be empty
 * but receive rounding residues (e.g. the sliding potential of a contact without sliding).
 */
static void en_check_tracker(const double* booked, const double expected, const double energy_scale,
                             const char* label, const char* tracker) {
    for (int k = 0; k < EN_NMON; k++) {
        bool ok = fabs(booked[k] - expected) <= EN_TOL * fabs(expected) + 1e-12 * energy_scale;
        if (!ok)
            printf("    %s: %s of monomer %d: booked %.17g, expected %.17g\n", label, tracker, k, booked[k], expected);
        CHECK(ok);
    }
}

/**
 * The system in device memory. Both state buffers hold the same state, because the kernels read the contact
 * rotation and twist from the 'next' buffer. The stored pointers are rotated by q and the contact rotation is q.
 */
struct EnDevice {
    DeviceMaterials      materials;
    DeviceState          curr;
    DeviceState          next;
    DeviceEnergy         energy;

    EnDevice(const EnSystem& s, const double4 q)
        : materials(EN_NMON), curr(EN_NMON), next(EN_NMON), energy(EN_NMON)
    {
        HostMaterials     host_materials(EN_NMON);
        HostMaterialsView hm = host_materials.view();
        for (int k = 0; k < EN_NMON; k++) {
            hm.radius[k]            = s.r[k];
            hm.mass[k]              = 1.;
            hm.moment[k]            = 1.;
            hm.matID[k]             = 0;
            hm.density[k]           = 3210.;
            hm.surface_energy[k]    = EN_GAMMA;
            hm.youngs_modulus[k]    = EN_E;
            hm.poisson_number[k]    = EN_NU;
            hm.damping_timescale[k] = EN_TVIS;
            hm.crit_rolling_disp[k] = EN_XI;
        }
        host_materials.push_to(materials, EN_NMON);

        HostState     host_state(EN_NMON);
        HostStateView hs = host_state.view();
        for (int k = 0; k < EN_NMON; k++) {
            hs.position[k] = s.x[k];
            hs.velocity[k] = { 0., 0., 0. };
            hs.omega[k]    = { 0., 0., 0. };
            hs.force[k]    = { 0., 0., 0. };
            hs.torque[k]   = { 0., 0., 0. };
        }
        for (int k = 0; k < EN_NMON * EN_NMON; k++) {
            hs.contact_compression[k] = -1.;
            hs.contact_twist[k]       = 0.;
            hs.contact_pointer[k]     = { 0., 0., 0. };
            hs.contact_normal[k]      = { 0., 0., 0. };
            hs.contact_rotation[k]    = { 0., 0., 0., 0. };
        }
        if (vec_length_sq(s.n[0]) != 0.) {
            hs.contact_pointer[EN_SLOT_01]     = quat_apply(q, s.n[0]);
            hs.contact_pointer[EN_SLOT_10]     = quat_apply(q, s.n[1]);
            hs.contact_rotation[EN_SLOT_01]    = q;
            hs.contact_rotation[EN_SLOT_10]    = q;
            hs.contact_twist[EN_SLOT_01]       = s.twist;
            hs.contact_twist[EN_SLOT_10]       = s.twist;
            hs.contact_compression[EN_SLOT_01] = s.compression_old;
            hs.contact_compression[EN_SLOT_10] = s.compression_old;
        }
        host_state.push_to(curr, EN_NMON);
        host_state.push_to(next, EN_NMON);

        energy.zero();
    }

    /** Checks the kernel launch, then that every energy tracker holds the expected energy per monomer. */
    void check_energy(const EnergyRecord& expected, const double energy_scale, const char* label) {
        cudaError_t launch_error = cudaGetLastError();
        cudaError_t sync_error   = cudaDeviceSynchronize();
        if (launch_error != cudaSuccess || sync_error != cudaSuccess)
            printf("    %s: kernel failed: %s / %s\n", label, cudaGetErrorString(launch_error), cudaGetErrorString(sync_error));
        CHECK(launch_error == cudaSuccess && sync_error == cudaSuccess);

        HostEnergy host_energy(EN_NMON);
        host_energy.pull_from(energy, EN_NMON);
        HostEnergyView b = host_energy.view();

        en_check_tracker(b.normal_pot,     expected.normal_pot,     energy_scale, label, "normal_pot");
        en_check_tracker(b.sliding_pot,    expected.sliding_pot,    energy_scale, label, "sliding_pot");
        en_check_tracker(b.rolling_pot,    expected.rolling_pot,    energy_scale, label, "rolling_pot");
        en_check_tracker(b.twisting_pot,   expected.twisting_pot,   energy_scale, label, "twisting_pot");
        en_check_tracker(b.normal_damp,    expected.normal_damp,    energy_scale, label, "normal_damp");
        en_check_tracker(b.sliding_slip,   expected.sliding_slip,   energy_scale, label, "sliding_slip");
        en_check_tracker(b.rolling_slip,   expected.rolling_slip,   energy_scale, label, "rolling_slip");
        en_check_tracker(b.twisting_slip,  expected.twisting_slip,  energy_scale, label, "twisting_slip");
        en_check_tracker(b.normal_break,   expected.normal_break,   energy_scale, label, "normal_break");
        en_check_tracker(b.sliding_break,  expected.sliding_break,  energy_scale, label, "sliding_break");
        en_check_tracker(b.rolling_break,  expected.rolling_break,  energy_scale, label, "rolling_break");
        en_check_tracker(b.twisting_break, expected.twisting_break, energy_scale, label, "twisting_break");
        en_check_tracker(b.normal_form,    expected.normal_form,    energy_scale, label, "normal_form");
    }

    void run_update_pointers() {
        DeviceStateView     cv = curr.view();
        DeviceStateView     nv = next.view();
        DeviceMaterialsView mv = materials.view();
        DeviceEnergyView    ev = energy.view();
        updatePointers<<<1, EN_NMON * EN_NMON>>>(
            cv.position, cv.contact_pointer,
            cv.contact_rotation, cv.contact_compression,
            nv.contact_pointer, nv.contact_rotation,
            nv.contact_twist, nv.contact_compression,
            ev.sliding_slip, ev.rolling_slip, ev.twisting_slip,
            ev.normal_break, ev.sliding_break, ev.rolling_break, ev.twisting_break,
            ev.normal_form,
            mv.radius, mv.youngs_modulus, mv.poisson_number,
            mv.surface_energy, mv.crit_rolling_disp,
            EN_NMON
        );
    }

    void run_evaluate() {
        DeviceStateView     cv = curr.view();
        DeviceStateView     nv = next.view();
        DeviceMaterialsView mv = materials.view();
        DeviceEnergyView    ev = energy.view();
        evaluate<<<1, EN_NMON * EN_NMON>>>(
            nv.position, cv.contact_pointer,
            nv.contact_rotation, nv.contact_twist, cv.contact_compression,
            nv.force, nv.torque,
            ev.normal_damp,
            mv.mass, mv.radius, mv.youngs_modulus, mv.poisson_number,
            mv.surface_energy, mv.crit_rolling_disp, mv.damping_timescale,
            EN_DT, EN_NMON
        );
    }

    void run_contact_potentials() {
        DeviceStateView     cv = curr.view();
        DeviceMaterialsView mv = materials.view();
        DeviceEnergyView    ev = energy.view();
        contact_potentials<<<1, EN_NMON * EN_NMON>>>(
            cv.position, cv.contact_pointer, cv.contact_rotation, cv.contact_twist,
            ev.normal_pot, ev.sliding_pot, ev.rolling_pot, ev.twisting_pot,
            mv.radius, mv.youngs_modulus, mv.poisson_number, mv.surface_energy,
            EN_NMON
        );
    }
};

/**
 * updatePointers: breaking, formation, inelastic motion and the cases in which nothing is booked.
 */
static void en_test_update_pointers(const double4 q) {
    // Breaking: everything stored in the contact is lost.
    {
        const char* label = "updatePointers, contact breaks";
        EnPair   p = en_pair(15e-9, 30e-9);
        EnSystem s = en_make_system(15e-9, 30e-9, -1.05 * p.delta_N_crit, true, 0.5, 0.5, 0.5);

        double normal; double3 sliding, rolling;
        en_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.normal_break   = 0.5 * get_U_N(p.F_c, p.delta_N_crit, get_contact_radius(normal, p.a_0, p.R), p.a_0);
        expected.sliding_break  = 0.5 * get_U_S(p.k_s, sliding);
        expected.rolling_break  = 0.5 * get_U_R(p.k_r, rolling);
        expected.twisting_break = 0.5 * get_U_T(p.k_t, s.twist);

        EnDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);

        HostState out(EN_NMON);
        out.pull_from(dev.next, EN_NMON);
        CHECK(vec_length_sq(out.view().contact_pointer[EN_SLOT_01]) == 0. && vec_length_sq(out.view().contact_pointer[EN_SLOT_10]) == 0.);
    }

    // Formation: the normal potential of the new contact is compensated.
    {
        const char* label = "updatePointers, contact forms";
        EnPair   p = en_pair(15e-9, 30e-9);
        EnSystem s = en_make_system(15e-9, 30e-9, p.delta_N_0, false, 0., 0., 0.);

        double normal; double3 sliding, rolling;
        en_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.normal_form = - 0.5 * get_U_N(p.F_c, p.delta_N_crit, get_contact_radius(normal, p.a_0, p.R), p.a_0);

        EnDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);

        HostState out(EN_NMON);
        out.pull_from(dev.next, EN_NMON);
        CHECK(vec_length_sq(out.view().contact_pointer[EN_SLOT_01]) != 0. && vec_length_sq(out.view().contact_pointer[EN_SLOT_10]) != 0.);
    }

    // Inelastic motion in all three degrees of freedom. Equal radii, for which the old pointer correction is
    // accurate (see CLAUDE.md (AC)). The booking does not depend on it.
    {
        const char* label = "updatePointers, inelastic motion";
        EnPair   p = en_pair(20e-9, 20e-9);
        EnSystem s = en_make_system(20e-9, 20e-9, p.delta_N_0, true, 2., 2., 2.);

        double normal; double3 sliding, rolling;
        en_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.sliding_slip  = 0.5 * p.k_s * p.delta_S_crit * (vec_length(sliding) - p.delta_S_crit);
        expected.rolling_slip  = 0.5 * p.k_r * p.delta_R_crit * (vec_length(rolling) - p.delta_R_crit);
        expected.twisting_slip = 0.5 * p.k_t * p.delta_T_crit * (fabs(s.twist)     - p.delta_T_crit);

        EnDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }

    // Below all critical displacements nothing is booked.
    {
        const char* label = "updatePointers, elastic contact";
        EnPair   p = en_pair(15e-9, 30e-9);
        EnSystem s = en_make_system(15e-9, 30e-9, p.delta_N_0, true, 0.5, 0.5, 0.5);

        EnergyRecord expected = {};

        EnDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }

    // Monomers apart without a contact: nothing is booked.
    {
        const char* label = "updatePointers, no contact";
        EnPair   p = en_pair(15e-9, 30e-9);
        EnSystem s = en_make_system(15e-9, 30e-9, -2. * p.delta_N_crit, false, 0., 0., 0.);

        EnergyRecord expected = {};

        EnDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }
}

/**
 * evaluate: the damping work and the force of the contact. It books no potentials.
 */
static void en_test_evaluate(const double4 q) {
    // A compressed contact, compressed further since the last step, with the pointers on the contact axis.
    {
        const char* label = "evaluate, damped normal motion";
        EnPair   p = en_pair(15e-9, 30e-9);
        EnSystem s = en_make_system(15e-9, 30e-9, 2. * p.delta_N_0, true, 0., 0., 0.);
        s.compression_old = 2. * p.delta_N_0 - 0.01 * p.delta_N_0;

        double normal; double3 sliding, rolling;
        en_displacements(s, normal, sliding, rolling);

        double a       = get_contact_radius(normal, p.a_0, p.R);
        double step    = normal - s.compression_old;
        double damping = 2. * EN_TVIS / (EN_NU * EN_NU) * p.E_s * a * step / EN_DT;

        EnergyRecord expected = {};
        expected.normal_damp = 0.5 * damping * step;

        EnDevice dev(s, q);
        dev.run_evaluate();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);

        // The force on monomer 0 points along (x_0 - x_1), the force on monomer 1 is opposite.
        HostState out(EN_NMON);
        out.pull_from(dev.next, EN_NMON);
        double3 F_0 = out.view().force[0];
        double3 F_1 = out.view().force[1];

        double3 diff = vec_diff(s.x[0], s.x[1]);
        vec_normalize(diff);
        double F_N = 4. * p.F_c * (pow(a / p.a_0, 3.) - pow(a / p.a_0, 1.5)) + damping;
        CHECK(fabs(vec_dot(F_0, diff) / F_N - 1.) < EN_TOL);
        CHECK(vec_length(vec_cross(F_0, diff)) < 1e-12 * fabs(F_N));
        CHECK(vec_length({ F_0.x + F_1.x, F_0.y + F_1.y, F_0.z + F_1.z }) < 1e-12 * fabs(F_N));
    }

    // A contact at equilibrium compression with sliding, rolling and twisting displacements and no normal motion:
    // nothing is booked.
    {
        const char* label = "evaluate, displaced contact";
        EnPair   p = en_pair(15e-9, 30e-9);
        EnSystem s = en_make_system(15e-9, 30e-9, p.delta_N_0, true, 0.5, 0.7, 0.3);

        double normal; double3 sliding, rolling;
        en_displacements(s, normal, sliding, rolling);
        s.compression_old = normal;

        EnergyRecord expected = {};

        EnDevice dev(s, q);
        dev.run_evaluate();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);

        HostState out(EN_NMON);
        out.pull_from(dev.next, EN_NMON);
        double3 F_0 = out.view().force[0];
        double3 F_1 = out.view().force[1];
        CHECK(vec_length({ F_0.x + F_1.x, F_0.y + F_1.y, F_0.z + F_1.z }) < 1e-10 * vec_length(F_0));
    }
}

/**
 * contact_potentials: the potentials of the stored state.
 */
static void en_test_contact_potentials(const double4 q) {
    // A contact with displacements in all four degrees of freedom.
    {
        const char* label = "contact_potentials, displaced contact";
        EnPair   p = en_pair(15e-9, 30e-9);
        EnSystem s = en_make_system(15e-9, 30e-9, 2. * p.delta_N_0, true, 0.5, 0.7, 0.3);

        double normal; double3 sliding, rolling;
        en_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.normal_pot   = 0.5 * get_U_N(p.F_c, p.delta_N_crit, get_contact_radius(normal, p.a_0, p.R), p.a_0);
        expected.sliding_pot  = 0.5 * get_U_S(p.k_s, sliding);
        expected.rolling_pot  = 0.5 * get_U_R(p.k_r, rolling);
        expected.twisting_pot = 0.5 * get_U_T(p.k_t, s.twist);

        EnDevice dev(s, q);
        dev.run_contact_potentials();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }

    // A contact stretched close to breaking, where the normal potential is large and positive.
    {
        const char* label = "contact_potentials, stretched contact";
        EnPair   p = en_pair(20e-9, 20e-9);
        EnSystem s = en_make_system(20e-9, 20e-9, -0.9 * p.delta_N_crit, true, 0.2, 0.2, 0.);

        double normal; double3 sliding, rolling;
        en_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.normal_pot  = 0.5 * get_U_N(p.F_c, p.delta_N_crit, get_contact_radius(normal, p.a_0, p.R), p.a_0);
        expected.sliding_pot = 0.5 * get_U_S(p.k_s, sliding);
        expected.rolling_pot = 0.5 * get_U_R(p.k_r, rolling);

        EnDevice dev(s, q);
        dev.run_contact_potentials();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }

    // Overlapping monomers without a registered contact store no potential energy.
    {
        const char* label = "contact_potentials, no contact";
        EnPair   p = en_pair(15e-9, 30e-9);
        EnSystem s = en_make_system(15e-9, 30e-9, p.delta_N_0, false, 0., 0., 0.);

        EnergyRecord expected = {};

        EnDevice dev(s, q);
        dev.run_contact_potentials();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }
}

void test_energy() {
    const double4 identity   = { 0., 0., 0., 1. };
    // A rotation by 0.7 rad around the axis (1, 2, 2) / 3.
    const double  half_angle = 0.35;
    const double4 rotation   = { sin(half_angle) / 3., 2. * sin(half_angle) / 3., 2. * sin(half_angle) / 3., cos(half_angle) };

    en_test_update_pointers(identity);
    en_test_update_pointers(rotation);

    en_test_evaluate(identity);
    en_test_evaluate(rotation);

    en_test_contact_potentials(identity);
    en_test_contact_potentials(rotation);
}
