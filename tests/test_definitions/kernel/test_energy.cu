/**
 * @file test_energy.cu
 * @brief Tests for the energy trackers (physics/energy.cuh): which kernel books which energy into which tracker.
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * Every test runs a kernel on a system of two monomers with one contact, or none (see two_monomer.cuh). Each
 * thread of the pair books half of the energy of the pair into the slot of its own monomer. So both slots of a
 * tracker must hold half of the expected pair energy, and all other trackers must stay empty.
 * The expected energies use the potentials of integrator_utils.cuh. The displacements they are evaluated at come
 * from tm_displacements, which is written out independently of the kernels.
 *
 * Every case runs with the contact pointers stored unrotated and rotated by 0.7 rad, the kernels have to undo
 * the rotation.
 *
 * The potentials are not booked by evaluate but by contact_potentials, which runs on the stored state at every
 * snapshot.
 */

#include <cmath>
#include "test_definitions/kernel/two_monomer.cuh"

/**
 * updatePointers: breaking, formation, inelastic motion and the cases in which nothing is booked.
 */
static void en_test_update_pointers(const double4 q) {
    // Breaking: everything stored in the contact is lost.
    {
        const char* label = "updatePointers, contact breaks";
        TmPair   p = tm_pair(15e-9, 30e-9);
        TmSystem s = tm_make_system(15e-9, 30e-9, -1.05 * p.delta_N_crit, true, 0.5, 0.5, 0.5);

        double normal; double3 sliding, rolling;
        tm_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.normal_break   = 0.5 * get_U_N(p.F_c, p.delta_N_crit, get_contact_radius(normal, p.a_0, p.R), p.a_0);
        expected.sliding_break  = 0.5 * get_U_S(p.k_s, sliding);
        expected.rolling_break  = 0.5 * get_U_R(p.k_r, rolling);
        expected.twisting_break = 0.5 * get_U_T(p.k_t, s.twist);

        TmDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);

        HostState out(TM_NMON);
        out.pull_from(dev.next, TM_NMON);
        CHECK(vec_length_sq(out.view().contact_pointer[TM_SLOT_01]) == 0. && vec_length_sq(out.view().contact_pointer[TM_SLOT_10]) == 0.);
    }

    // Formation: the normal potential of the new contact is compensated.
    {
        const char* label = "updatePointers, contact forms";
        TmPair   p = tm_pair(15e-9, 30e-9);
        TmSystem s = tm_make_system(15e-9, 30e-9, p.delta_N_0, false, 0., 0., 0.);

        double normal; double3 sliding, rolling;
        tm_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.normal_form = - 0.5 * get_U_N(p.F_c, p.delta_N_crit, get_contact_radius(normal, p.a_0, p.R), p.a_0);

        TmDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);

        HostState out(TM_NMON);
        out.pull_from(dev.next, TM_NMON);
        CHECK(vec_length_sq(out.view().contact_pointer[TM_SLOT_01]) != 0. && vec_length_sq(out.view().contact_pointer[TM_SLOT_10]) != 0.);
    }

    // Inelastic motion in all three degrees of freedom. Equal radii, for which the old pointer correction is
    // accurate (see CLAUDE.md (AC)). The booking does not depend on it.
    {
        const char* label = "updatePointers, inelastic motion";
        TmPair   p = tm_pair(20e-9, 20e-9);
        TmSystem s = tm_make_system(20e-9, 20e-9, p.delta_N_0, true, 2., 2., 2.);

        double normal; double3 sliding, rolling;
        tm_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.sliding_slip  = 0.5 * p.k_s * p.delta_S_crit * (vec_length(sliding) - p.delta_S_crit);
        expected.rolling_slip  = 0.5 * p.k_r * p.delta_R_crit * (vec_length(rolling) - p.delta_R_crit);
        expected.twisting_slip = 0.5 * p.k_t * p.delta_T_crit * (fabs(s.twist)     - p.delta_T_crit);

        TmDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }

    // Below all critical displacements nothing is booked.
    {
        const char* label = "updatePointers, elastic contact";
        TmPair   p = tm_pair(15e-9, 30e-9);
        TmSystem s = tm_make_system(15e-9, 30e-9, p.delta_N_0, true, 0.5, 0.5, 0.5);

        EnergyRecord expected = {};

        TmDevice dev(s, q);
        dev.run_update_pointers();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }

    // Monomers apart without a contact: nothing is booked.
    {
        const char* label = "updatePointers, no contact";
        TmPair   p = tm_pair(15e-9, 30e-9);
        TmSystem s = tm_make_system(15e-9, 30e-9, -2. * p.delta_N_crit, false, 0., 0., 0.);

        EnergyRecord expected = {};

        TmDevice dev(s, q);
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
        TmPair   p = tm_pair(15e-9, 30e-9);
        TmSystem s = tm_make_system(15e-9, 30e-9, 2. * p.delta_N_0, true, 0., 0., 0.);
        s.compression_old = 2. * p.delta_N_0 - 0.01 * p.delta_N_0;

        double normal; double3 sliding, rolling;
        tm_displacements(s, normal, sliding, rolling);

        double a       = get_contact_radius(normal, p.a_0, p.R);
        double step    = normal - s.compression_old;
        double damping = 2. * TM_TVIS / (TM_NU * TM_NU) * p.E_s * a * step / TM_DT;

        EnergyRecord expected = {};
        expected.normal_damp = 0.5 * damping * step;

        TmDevice dev(s, q);
        dev.run_evaluate();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);

        // The force on monomer 0 points along (x_0 - x_1), the force on monomer 1 is opposite.
        HostState out(TM_NMON);
        out.pull_from(dev.next, TM_NMON);
        double3 F_0 = out.view().force[0];
        double3 F_1 = out.view().force[1];

        double3 diff = vec_diff(s.x[0], s.x[1]);
        vec_normalize(diff);
        double F_N = 4. * p.F_c * (pow(a / p.a_0, 3.) - pow(a / p.a_0, 1.5)) + damping;
        CHECK(fabs(vec_dot(F_0, diff) / F_N - 1.) < TM_TOL);
        CHECK(vec_length(vec_cross(F_0, diff)) < 1e-12 * fabs(F_N));
        CHECK(vec_length({ F_0.x + F_1.x, F_0.y + F_1.y, F_0.z + F_1.z }) < 1e-12 * fabs(F_N));
    }

    // A contact at equilibrium compression with sliding, rolling and twisting displacements and no normal motion:
    // nothing is booked.
    {
        const char* label = "evaluate, displaced contact";
        TmPair   p = tm_pair(15e-9, 30e-9);
        TmSystem s = tm_make_system(15e-9, 30e-9, p.delta_N_0, true, 0.5, 0.7, 0.3);

        double normal; double3 sliding, rolling;
        tm_displacements(s, normal, sliding, rolling);
        s.compression_old = normal;

        EnergyRecord expected = {};

        TmDevice dev(s, q);
        dev.run_evaluate();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);

        HostState out(TM_NMON);
        out.pull_from(dev.next, TM_NMON);
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
        TmPair   p = tm_pair(15e-9, 30e-9);
        TmSystem s = tm_make_system(15e-9, 30e-9, 2. * p.delta_N_0, true, 0.5, 0.7, 0.3);

        double normal; double3 sliding, rolling;
        tm_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.normal_pot   = 0.5 * get_U_N(p.F_c, p.delta_N_crit, get_contact_radius(normal, p.a_0, p.R), p.a_0);
        expected.sliding_pot  = 0.5 * get_U_S(p.k_s, sliding);
        expected.rolling_pot  = 0.5 * get_U_R(p.k_r, rolling);
        expected.twisting_pot = 0.5 * get_U_T(p.k_t, s.twist);

        TmDevice dev(s, q);
        dev.run_contact_potentials();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }

    // A contact stretched close to breaking, where the normal potential is large and positive.
    {
        const char* label = "contact_potentials, stretched contact";
        TmPair   p = tm_pair(20e-9, 20e-9);
        TmSystem s = tm_make_system(20e-9, 20e-9, -0.9 * p.delta_N_crit, true, 0.2, 0.2, 0.);

        double normal; double3 sliding, rolling;
        tm_displacements(s, normal, sliding, rolling);

        EnergyRecord expected = {};
        expected.normal_pot  = 0.5 * get_U_N(p.F_c, p.delta_N_crit, get_contact_radius(normal, p.a_0, p.R), p.a_0);
        expected.sliding_pot = 0.5 * get_U_S(p.k_s, sliding);
        expected.rolling_pot = 0.5 * get_U_R(p.k_r, rolling);

        TmDevice dev(s, q);
        dev.run_contact_potentials();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }

    // Overlapping monomers without a registered contact store no potential energy.
    {
        const char* label = "contact_potentials, no contact";
        TmPair   p = tm_pair(15e-9, 30e-9);
        TmSystem s = tm_make_system(15e-9, 30e-9, p.delta_N_0, false, 0., 0., 0.);

        EnergyRecord expected = {};

        TmDevice dev(s, q);
        dev.run_contact_potentials();
        dev.check_energy(expected, p.F_c * p.delta_N_crit, label);
    }
}

void test_energy() {
    en_test_update_pointers(TM_IDENTITY);
    en_test_update_pointers(TM_ROTATION);

    en_test_evaluate(TM_IDENTITY);
    en_test_evaluate(TM_ROTATION);

    en_test_contact_potentials(TM_IDENTITY);
    en_test_contact_potentials(TM_ROTATION);
}
