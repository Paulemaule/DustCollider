/**
 * @file two_monomer.cuh
 * @brief Shared fixture for the kernel tests: a system of two forsterite monomers with one contact, or none.
 *
 * It provides:
 *  - TmPair / tm_pair:            the constants of the monomer pair, from the same helpers the kernels use.
 *  - TmSystem / tm_make_system:   two monomers with a prescribed normal, sliding, rolling and twisting displacement.
 *  - tm_displacements:            the displacements of thread (0,1), written out independently of the kernels.
 *  - TmDevice:                    the system in device memory, the kernel launches and the energy tracker checks.
 *  - TM_IDENTITY / TM_ROTATION:   the contact rotations every case runs with. The stored pointers are rotated by
 *                                 them, so the kernels have to undo the rotation.
 *
 * Thread (0,1) works on the pair matrix slot 0 + 1 * Nmon, thread (1,0) on the slot 1 + 0 * Nmon.
 *
 * This file is #included by the kernel test files; testkit.h is expected to be already in scope.
 */

#pragma once

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

// Relative tolerance of the booked energies and forces.
static const double TM_TOL = 1e-9;

// Forsterite, as in the shipped command files.
static const double TM_GAMMA = 0.07;
static const double TM_E     = 204e9;
static const double TM_NU    = 0.24;
static const double TM_XI    = 2e-10;
static const double TM_TVIS  = 1e-12;

// A timestep of the order the auto-timestep picks for these monomers.
static const double TM_DT    = 2e-14;

// The test systems consist of two monomers. Thread (0,1) works on the pair matrix slot 0 + 1 * Nmon,
// thread (1,0) on the slot 1 + 0 * Nmon.
static const int TM_NMON    = 2;
static const int TM_SLOT_01 = 2;
static const int TM_SLOT_10 = 1;

// The contact rotations: the identity, and a rotation by 0.7 rad around the axis (1, 2, 2) / 3.
static const double4 TM_IDENTITY = { 0., 0., 0., 1. };
static const double4 TM_ROTATION = { sin(0.35) / 3., 2. * sin(0.35) / 3., 2. * sin(0.35) / 3., cos(0.35) };

/**
 * The constants of a pair of forsterite monomers.
 */
struct TmPair {
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

static TmPair tm_pair(const double r0, const double r1) {
    TmPair p;
    double G     = get_G_i(TM_E, TM_NU);
    double gamma = get_gamma(TM_GAMMA, TM_GAMMA);

    p.R            = get_R(r0, r1);
    p.E_s          = get_E_s(TM_E, TM_E, TM_NU, TM_NU);
    p.a_0          = get_a_0(gamma, p.R, p.E_s);
    p.F_c          = 3. * PI * gamma * p.R;
    p.k_s          = 8. * get_G_s(G, G, TM_NU, TM_NU) * p.a_0;
    p.k_r          = 4. * p.F_c / p.R;
    p.k_t          = 16. * (G * G / (G + G)) * p.a_0 * p.a_0 * p.a_0 / 3.;
    p.delta_N_0    = p.a_0 * p.a_0 / (3. * p.R);
    p.delta_N_crit = get_delta_N_crit(p.a_0, p.R);
    p.delta_S_crit = get_delta_S_crit(TM_NU, TM_NU, p.a_0);
    p.delta_R_crit = TM_XI;
    p.delta_T_crit = 1. / (16. * PI);
    return p;
}

/**
 * A system of two monomers, with or without a contact.
 */
struct TmSystem {
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
static TmSystem tm_make_system(const double r0, const double r1, const double delta, const bool contact,
                               const double f_S, const double f_R, const double f_T) {
    TmPair p = tm_pair(r0, r1);

    // The axis from monomer 0 to monomer 1 and two tangential directions.
    double3 axis = { 1., 0.3, -0.2 };
    vec_normalize(axis);
    double3 e1 = vec_cross(axis, { 0., 0., 1. });
    vec_normalize(e1);
    double3 e2 = vec_cross(axis, e1);

    TmSystem s;
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
static void tm_displacements(const TmSystem& s, double& normal, double3& sliding, double3& rolling) {
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
static void tm_check_tracker(const double* booked, const double expected, const double energy_scale,
                             const char* label, const char* tracker) {
    for (int k = 0; k < TM_NMON; k++) {
        bool ok = fabs(booked[k] - expected) <= TM_TOL * fabs(expected) + 1e-12 * energy_scale;
        if (!ok)
            printf("    %s: %s of monomer %d: booked %.17g, expected %.17g\n", label, tracker, k, booked[k], expected);
        CHECK(ok);
    }
}

/**
 * The system in device memory. Both state buffers hold the same state, because the kernels read the contact
 * rotation and twist from the 'next' buffer. The stored pointers are rotated by q and the contact rotation is q.
 */
struct TmDevice {
    DeviceMaterials      materials;
    DeviceState          curr;
    DeviceState          next;
    DeviceEnergy         energy;

    TmDevice(const TmSystem& s, const double4 q)
        : materials(TM_NMON), curr(TM_NMON), next(TM_NMON), energy(TM_NMON)
    {
        HostMaterials     host_materials(TM_NMON);
        HostMaterialsView hm = host_materials.view();
        for (int k = 0; k < TM_NMON; k++) {
            hm.radius[k]            = s.r[k];
            hm.mass[k]              = 1.;
            hm.moment[k]            = 1.;
            hm.matID[k]             = 0;
            hm.density[k]           = 3210.;
            hm.surface_energy[k]    = TM_GAMMA;
            hm.youngs_modulus[k]    = TM_E;
            hm.poisson_number[k]    = TM_NU;
            hm.damping_timescale[k] = TM_TVIS;
            hm.crit_rolling_disp[k] = TM_XI;
        }
        host_materials.push_to(materials, TM_NMON);

        HostState     host_state(TM_NMON);
        HostStateView hs = host_state.view();
        for (int k = 0; k < TM_NMON; k++) {
            hs.position[k] = s.x[k];
            hs.velocity[k] = { 0., 0., 0. };
            hs.omega[k]    = { 0., 0., 0. };
            hs.force[k]    = { 0., 0., 0. };
            hs.torque[k]   = { 0., 0., 0. };
        }
        for (int k = 0; k < TM_NMON * TM_NMON; k++) {
            hs.contact_compression[k] = -1.;
            hs.contact_twist[k]       = 0.;
            hs.contact_pointer[k]     = { 0., 0., 0. };
            hs.contact_normal[k]      = { 0., 0., 0. };
            hs.contact_rotation[k]    = { 0., 0., 0., 0. };
        }
        if (vec_length_sq(s.n[0]) != 0.) {
            hs.contact_pointer[TM_SLOT_01]     = quat_apply(q, s.n[0]);
            hs.contact_pointer[TM_SLOT_10]     = quat_apply(q, s.n[1]);
            hs.contact_rotation[TM_SLOT_01]    = q;
            hs.contact_rotation[TM_SLOT_10]    = q;
            hs.contact_twist[TM_SLOT_01]       = s.twist;
            hs.contact_twist[TM_SLOT_10]       = s.twist;
            hs.contact_compression[TM_SLOT_01] = s.compression_old;
            hs.contact_compression[TM_SLOT_10] = s.compression_old;
        }
        host_state.push_to(curr, TM_NMON);
        host_state.push_to(next, TM_NMON);

        energy.zero();
    }

    /** Checks the kernel launch, then that every energy tracker holds the expected energy per monomer. */
    void check_energy(const EnergyRecord& expected, const double energy_scale, const char* label) {
        cudaError_t launch_error = cudaGetLastError();
        cudaError_t sync_error   = cudaDeviceSynchronize();
        if (launch_error != cudaSuccess || sync_error != cudaSuccess)
            printf("    %s: kernel failed: %s / %s\n", label, cudaGetErrorString(launch_error), cudaGetErrorString(sync_error));
        CHECK(launch_error == cudaSuccess && sync_error == cudaSuccess);

        HostEnergy host_energy(TM_NMON);
        host_energy.pull_from(energy, TM_NMON);
        HostEnergyView b = host_energy.view();

        tm_check_tracker(b.normal_pot,     expected.normal_pot,     energy_scale, label, "normal_pot");
        tm_check_tracker(b.sliding_pot,    expected.sliding_pot,    energy_scale, label, "sliding_pot");
        tm_check_tracker(b.rolling_pot,    expected.rolling_pot,    energy_scale, label, "rolling_pot");
        tm_check_tracker(b.twisting_pot,   expected.twisting_pot,   energy_scale, label, "twisting_pot");
        tm_check_tracker(b.normal_damp,    expected.normal_damp,    energy_scale, label, "normal_damp");
        tm_check_tracker(b.sliding_slip,   expected.sliding_slip,   energy_scale, label, "sliding_slip");
        tm_check_tracker(b.rolling_slip,   expected.rolling_slip,   energy_scale, label, "rolling_slip");
        tm_check_tracker(b.twisting_slip,  expected.twisting_slip,  energy_scale, label, "twisting_slip");
        tm_check_tracker(b.normal_break,   expected.normal_break,   energy_scale, label, "normal_break");
        tm_check_tracker(b.sliding_break,  expected.sliding_break,  energy_scale, label, "sliding_break");
        tm_check_tracker(b.rolling_break,  expected.rolling_break,  energy_scale, label, "rolling_break");
        tm_check_tracker(b.twisting_break, expected.twisting_break, energy_scale, label, "twisting_break");
        tm_check_tracker(b.normal_form,    expected.normal_form,    energy_scale, label, "normal_form");
    }

    void run_update_pointers() {
        DeviceStateView     cv = curr.view();
        DeviceStateView     nv = next.view();
        DeviceMaterialsView mv = materials.view();
        DeviceEnergyView    ev = energy.view();
        updatePointers<<<1, TM_NMON * TM_NMON>>>(
            cv.position, cv.contact_pointer,
            cv.contact_rotation, cv.contact_compression,
            nv.contact_pointer, nv.contact_rotation,
            nv.contact_twist, nv.contact_compression,
            ev.sliding_slip, ev.rolling_slip, ev.twisting_slip,
            ev.normal_break, ev.sliding_break, ev.rolling_break, ev.twisting_break,
            ev.normal_form,
            mv.radius, mv.youngs_modulus, mv.poisson_number,
            mv.surface_energy, mv.crit_rolling_disp,
            TM_NMON
        );
    }

    void run_evaluate() {
        DeviceStateView     cv = curr.view();
        DeviceStateView     nv = next.view();
        DeviceMaterialsView mv = materials.view();
        DeviceEnergyView    ev = energy.view();
        evaluate<<<1, TM_NMON * TM_NMON>>>(
            nv.position, cv.contact_pointer,
            nv.contact_rotation, nv.contact_twist, cv.contact_compression,
            nv.force, nv.torque,
            ev.normal_damp,
            mv.mass, mv.radius, mv.youngs_modulus, mv.poisson_number,
            mv.surface_energy, mv.crit_rolling_disp, mv.damping_timescale,
            TM_DT, TM_NMON
        );
    }

    void run_contact_potentials() {
        DeviceStateView     cv = curr.view();
        DeviceMaterialsView mv = materials.view();
        DeviceEnergyView    ev = energy.view();
        contact_potentials<<<1, TM_NMON * TM_NMON>>>(
            cv.position, cv.contact_pointer, cv.contact_rotation, cv.contact_twist,
            ev.normal_pot, ev.sliding_pot, ev.rolling_pot, ev.twisting_pot,
            mv.radius, mv.youngs_modulus, mv.poisson_number, mv.surface_energy,
            TM_NMON
        );
    }
};
