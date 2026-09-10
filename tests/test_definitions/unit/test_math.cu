/**
 * @file test_math.cu
 * @brief Tests for the __host__ __device__ math helpers in physics/integrator_utils.cuh
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * Every helper here is __host__ __device__, so it can be exercised directly on
 * the CPU. The cluster helpers (dfs_iterative / findMonomerClusters) and the
 * CALC_MONOMER_INDICES macro from the same header are covered elsewhere
 * (test_cluster.cu) and are not retested here.
 */

#include "utils/constant.cuh"
#include "utils/vector.cuh"
#include "physics/integrator_utils.cuh"
#include <cmath>

void test_math() {

    // ------------------------------------------------------------------ //
    // get_R: reduced radius  r_i * r_j / (r_i + r_j)
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_R(1.0, 1.0),   0.5,  1e-14);
    CHECK_APPROX(get_R(2.0, 2.0),   1.0,  1e-14);
    // 2*6 / (2+6) = 12/8 = 1.5
    CHECK_APPROX(get_R(2.0, 6.0),   1.5,  1e-14);
    // Symmetric in its arguments.
    CHECK_APPROX(get_R(6.0, 2.0),   get_R(2.0, 6.0),  1e-14);
    // Large-disparity limit: r_i >> r_j  ->  result approaches r_j.
    CHECK_APPROX(get_R(1e6, 1.0),   1.0,  1e-3);
    // Realistic nm-scale radii must not under/overflow: equal radii -> r/2.
    CHECK_APPROX(get_R(1e-9, 1e-9), 0.5e-9,  1e-24);

    // ------------------------------------------------------------------ //
    // get_G_i: shear modulus  E / (2 * (1 + nu))
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_G_i(2.0, 0.0),  1.0,  1e-14);  // 2/(2*1) = 1
    CHECK_APPROX(get_G_i(3.0, 0.5),  1.0,  1e-14);  // 3/(2*1.5) = 1   (incompressible limit -> E/3)
    CHECK_APPROX(get_G_i(5.0, 0.25), 2.0,  1e-14);  // 5/(2*1.25) = 2  (typical Poisson ratio)
    // Realistic Young's modulus must not lose precision.
    CHECK_APPROX(get_G_i(1e11, 0.25), 4e10, 1e-4);

    // ------------------------------------------------------------------ //
    // get_E_s: combined Young's modulus
    //   1 / ((1-nu_i^2)/E_i + (1-nu_j^2)/E_j)
    // ------------------------------------------------------------------ //

    // Equal materials, nu=0 -> E/2.
    CHECK_APPROX(get_E_s(4.0, 4.0, 0.0, 0.0),  2.0,  1e-14);
    // Equal materials, non-zero nu -> E / (2*(1-nu^2)).  E=4, nu=0.5 -> 4/(2*0.75).
    CHECK_APPROX(get_E_s(4.0, 4.0, 0.5, 0.5),  4.0 / 1.5,  1e-14);
    // Unequal moduli, nu=0 -> harmonic-style combination E_i*E_j/(E_i+E_j).
    CHECK_APPROX(get_E_s(2.0, 6.0, 0.0, 0.0),  1.5,  1e-14);
    // Unequal Poisson ratios: 1 / (0.99/2 + 0.84/8) = 1 / 0.6.
    CHECK_APPROX(get_E_s(2.0, 8.0, 0.1, 0.4),  1.0 / 0.6,  1e-14);
    // Symmetric under simultaneous (E,nu) swap of i and j.
    CHECK_APPROX(get_E_s(2.0, 8.0, 0.1, 0.4),  get_E_s(8.0, 2.0, 0.4, 0.1),  1e-14);

    // ------------------------------------------------------------------ //
    // get_G_s: combined shear modulus
    //   1 / G* = (2 - nu_i) / G_i + (2 - nu_j) / G_j
    // ------------------------------------------------------------------ //

    // Equal materials, nu=0 -> 1/((2/G) + (2/G)) = G/4.  G=4 -> 1.
    CHECK_APPROX(get_G_s(4.0, 4.0, 0.0, 0.0),  1.0,  1e-14);
    // Equal materials, nu=0.5 -> 1/(2 * (1.5/4)) = 4/3.
    CHECK_APPROX(get_G_s(4.0, 4.0, 0.5, 0.5),  4.0 / 3.0,  1e-14);
    // Unequal: 1/((2-0.1)/2 + (2-0.4)/8) = 1/(0.95 + 0.2) = 1/1.15.
    CHECK_APPROX(get_G_s(2.0, 8.0, 0.1, 0.4),  1.0 / 1.15,  1e-14);
    // Symmetric under simultaneous (G, nu) swap of i and j.
    CHECK_APPROX(get_G_s(2.0, 8.0, 0.1, 0.4),  get_G_s(8.0, 2.0, 0.4, 0.1),  1e-14);
    // Sanity check
    CHECK(get_G_s(3.0, 7.0, 0.2, 0.35) != get_E_s(3.0, 7.0, 0.2, 0.35));

    // ------------------------------------------------------------------ //
    // get_gamma: combined surface energy
    //   gamma_i + gamma_j - 2 / (1/gamma_i + 1/gamma_j)
    //   closed form: (gamma_i^2 + gamma_j^2) / (gamma_i + gamma_j)
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_gamma(1.0, 1.0),  1.0,  1e-14);
    CHECK_APPROX(get_gamma(2.0, 2.0),  2.0,  1e-14);
    // 1 + 3 - 2*1*3/(1+3) = 4 - 1.5 = 2.5
    CHECK_APPROX(get_gamma(1.0, 3.0),  2.5,  1e-14);
    // Symmetric in its arguments.
    CHECK_APPROX(get_gamma(3.0, 1.0),  get_gamma(1.0, 3.0),  1e-14);
    // Strongly asymmetric, checked against the closed form: (4+64)/(2+8) = 6.8.
    CHECK_APPROX(get_gamma(2.0, 8.0),  6.8,  1e-14);
    // Realistic surface energies (~0.01-0.1 J/m^2): (0.01^2+0.03^2)/(0.01+0.03).
    CHECK_APPROX(get_gamma(0.01, 0.03),
                 (0.01 * 0.01 + 0.03 * 0.03) / (0.01 + 0.03),  1e-16);

    // ------------------------------------------------------------------ //
    // get_a_0: equilibrium contact radius  (9*PI*gamma*R^2 / E_s)^(1/3)
    // ------------------------------------------------------------------ //

    // gamma=1, R=1, E_s=9*PI  ->  a_0 = 1
    CHECK_APPROX(get_a_0(1.0, 1.0, 9.0 * PI),  1.0,  1e-12);
    // Cube-root scaling in gamma: doubling gamma -> a_0 * 2^(1/3).
    CHECK_APPROX(get_a_0(2.0, 1.0, 9.0 * PI),  std::cbrt(2.0),  1e-12);
    // R enters as R^2 under the cube root: R=2 -> (4)^(1/3).
    CHECK_APPROX(get_a_0(1.0, 2.0, 9.0 * PI),  std::cbrt(4.0),  1e-12);
    // Inverse cube-root scaling in E_s: E_s/8 -> a_0 * 2.
    CHECK_APPROX(get_a_0(1.0, 1.0, 9.0 * PI / 8.0),  2.0,  1e-12);
    // Realistic micro-scale inputs must not overflow.
    CHECK(std::isfinite(get_a_0(0.05, 1e-7, 1e11)));

    // ------------------------------------------------------------------ //
    // get_delta_N_crit: critical normal displacement
    //   0.5 * a_0^2 / (R * cbrt(6))
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_delta_N_crit(1.0, 1.0),  0.5 / std::cbrt(6.0),  1e-14);
    // a_0=2, R=1 -> 0.5*4 / cbrt(6) = 2/cbrt(6)
    CHECK_APPROX(get_delta_N_crit(2.0, 1.0),  2.0 / std::cbrt(6.0),  1e-14);
    // Inverse scaling in R: a_0=1, R=2 -> 0.25/cbrt(6).
    CHECK_APPROX(get_delta_N_crit(1.0, 2.0),  0.25 / std::cbrt(6.0),  1e-14);
    // a_0^2 scaling: a_0=3, R=1 -> 4.5/cbrt(6).
    CHECK_APPROX(get_delta_N_crit(3.0, 1.0),  4.5 / std::cbrt(6.0),  1e-14);

    // ------------------------------------------------------------------ //
    // get_delta_S_crit: critical sliding displacement
    //   (2 - 0.5*(nu_i+nu_j)) * a_0 / (16*PI)
    // ------------------------------------------------------------------ //

    // nu=0, a_0=16*PI -> result = 2
    CHECK_APPROX(get_delta_S_crit(0.0, 0.0, 16.0 * PI),  2.0,  1e-12);
    // nu_i=nu_j=0.5, a_0=16*PI -> (2 - 0.5) = 1.5
    CHECK_APPROX(get_delta_S_crit(0.5, 0.5, 16.0 * PI),  1.5,  1e-12);
    // Asymmetric nu: (2 - 0.5*0.6) = 1.7
    CHECK_APPROX(get_delta_S_crit(0.2, 0.4, 16.0 * PI),  1.7,  1e-12);
    // Symmetric in the two Poisson ratios.
    CHECK_APPROX(get_delta_S_crit(0.4, 0.2, 16.0 * PI),
                 get_delta_S_crit(0.2, 0.4, 16.0 * PI),  1e-14);
    // Linear in a_0: doubling a_0 doubles the result.
    CHECK_APPROX(get_delta_S_crit(0.0, 0.0, 32.0 * PI),  4.0,  1e-12);

    // ------------------------------------------------------------------ //
    // get_U_N: normal (JKR) potential
    //   F_c * delta_N_crit * (C0 + 4*cbrt(6)*( (4/5)(a/a_0)^5
    //                         - (4/3)(a/a_0)^(7/2) + (1/3)(a/a_0)^2 ))
    //   with C0 = 0.84661389438303971.
    // ------------------------------------------------------------------ //

    // At a=0 the polynomial vanishes, leaving the constant prefactor C0.
    CHECK_APPROX(get_U_N(1.0, 1.0, 0.0, 1.0),  0.84661389438303971,  1e-14);
    // Linear in F_c (delta_N_crit, a/a_0 fixed).
    CHECK_APPROX(get_U_N(2.0, 1.0, 0.0, 1.0),  2.0 * 0.84661389438303971,  1e-14);
    // Linear in delta_N_crit (F_c, a/a_0 fixed).
    CHECK_APPROX(get_U_N(1.0, 3.0, 0.0, 1.0),  3.0 * 0.84661389438303971,  1e-14);
    // Linear in the F_c*delta_N_crit product jointly.
    CHECK_APPROX(get_U_N(2.0, 3.0, 0.0, 1.0),  6.0 * 0.84661389438303971,  1e-14);
    // At a/a_0 = 1 the bracketed sum is 4/5 - 4/3 + 1/3 = -1/5, so the whole
    // expression reduces to C0 - 0.8*cbrt(6).
    CHECK_APPROX(get_U_N(1.0, 1.0, 1.0, 1.0),
                 0.84661389438303971 - 0.8 * std::pow(6.0, 1.0 / 3.0),  1e-14);
    // The potential depends on a and a_0 only through their ratio.
    CHECK_APPROX(get_U_N(1.0, 1.0, 2.0, 2.0),  get_U_N(1.0, 1.0, 1.0, 1.0),  1e-14);
    CHECK_APPROX(get_U_N(1.0, 1.0, 1.0, 2.0),  get_U_N(1.0, 1.0, 2.0, 4.0),  1e-14);

    // ------------------------------------------------------------------ //
    // get_U_S and get_U_R: quadratic sliding / rolling potential
    //   0.5 * k * |disp|^2
    // ------------------------------------------------------------------ //

    // k=2, {3,4,0} -> 0.5*2*25 = 25
    CHECK_APPROX(get_U_S(2.0, {3.0, 4.0, 0.0}),  25.0,  1e-14);
    CHECK_APPROX(get_U_R(2.0, {3.0, 4.0, 0.0}),  25.0,  1e-14);
    // Genuine 3-D displacement: |{1,2,2}|^2 = 9 -> 0.5*2*9 = 9.
    CHECK_APPROX(get_U_S(2.0, {1.0, 2.0, 2.0}),  9.0,  1e-14);
    CHECK_APPROX(get_U_R(2.0, {1.0, 2.0, 2.0}),  9.0,  1e-14);
    // Zero displacement and zero stiffness both give zero energy.
    CHECK_APPROX(get_U_S(5.0, {0.0, 0.0, 0.0}),  0.0,  1e-14);
    CHECK_APPROX(get_U_S(0.0, {1.0, 2.0, 3.0}),  0.0,  1e-14);
    // Sign of the displacement components is irrelevant (it is squared).
    CHECK_APPROX(get_U_S(2.0, {-3.0, 4.0, 0.0}),  25.0,  1e-14);
    // Linear in the stiffness: doubling k doubles the energy.
    CHECK_APPROX(get_U_S(4.0, {3.0, 4.0, 0.0}),  50.0,  1e-14);

    // ------------------------------------------------------------------ //
    // get_U_T: quadratic twisting potential  0.5 * k * twist^2
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_U_T(2.0, 3.0),   9.0,  1e-14);
    // Zero twist -> zero energy.
    CHECK_APPROX(get_U_T(5.0, 0.0),   0.0,  1e-14);
    // Sign of the twist is irrelevant (it is squared).
    CHECK_APPROX(get_U_T(2.0, -3.0),  9.0,  1e-14);
    // Linear in the stiffness.
    CHECK_APPROX(get_U_T(4.0, 3.0),   18.0,  1e-14);

    // ------------------------------------------------------------------ //
    // get_contact_radius: Newton solver for the JKR contact radius.
    //
    // Solves  3x^2 - 2*sqrt(x) - y = 0  with x = a/a_0 and
    // y = delta_N / delta_N_0,  delta_N_0 = a_0^2 / (3R).
    // The initial guess is x = 1.
    //
    // Below the critical displacement the result is clamped to
    // a_crit = (1/6)^(2/3) * a_0 (the location of the polynomial's minimum).
    // ------------------------------------------------------------------ //

    // Equilibrium: delta_N = delta_N_0 (y=1) is an exact root at x=1 -> a = a_0.
    {
        const double a_0 = 1.0, R = 1.0;
        const double delta_N_0 = a_0 * a_0 / (3.0 * R);
        CHECK_APPROX(get_contact_radius(delta_N_0, a_0, R),  a_0,  1e-14);
    }

    // Zero displacement (y=0): a = (2/3)^(2/3) * a_0.
    {
        CHECK_APPROX(get_contact_radius(0.0, 1.0, 1.0),
                     std::pow(2.0 / 3.0, 2.0 / 3.0),  1e-10);
    }

    // Compressive branch (x > 1): pick x=2, feed the delta_N that maps to it,
    // and confirm the solver recovers a = 2*a_0. This actually exercises the
    // Newton iteration (unlike the y=1 case, where the guess is already a root).
    {
        const double a_0 = 1.0, R = 1.0, x_target = 2.0;
        const double delta_N_0 = a_0 * a_0 / (3.0 * R);
        const double delta_N = (3.0 * x_target * x_target - 2.0 * std::sqrt(x_target)) * delta_N_0;
        CHECK_APPROX(get_contact_radius(delta_N, a_0, R),  x_target * a_0,  1e-9);
    }

    // Lower valid branch (x < 1): x=0.5 maps to y = -0.664..., which is above
    // the displacement-controlled pull-off threshold y_crit = -(9/16)^(1/3) =
    // -0.8255..., so the solver still runs and must converge to a = 0.5*a_0.
    {
        const double a_0 = 1.0, R = 1.0, x_target = 0.5;
        const double delta_N_0 = a_0 * a_0 / (3.0 * R);
        const double delta_N = (3.0 * x_target * x_target - 2.0 * std::sqrt(x_target)) * delta_N_0;
        CHECK_APPROX(get_contact_radius(delta_N, a_0, R),  x_target * a_0,  1e-9);
    }

    // Displacement-controlled JKR pull-off (correct physical behaviour).
    //
    // The contact breaks at the MINIMUM of y(x) = 3x^2 - 2*sqrt(x), i.e. at
    //   delta_N_crit = -(9/16)^(1/3) * delta_N_0  ~ -0.8255 * delta_N_0,
    // where the contact radius is exactly a_crit = (1/6)^(2/3) * a_0. The radius
    // therefore decreases CONTINUOUSLY down to a_crit and only then clamps —
    // there is no jump.
    //
    // NOTE: the current source uses the exponent (9/16)^(2/3) ~ 0.6814 for the
    // threshold, which clamps prematurely and introduces a spurious
    // discontinuity. The two in-window checks below encode the intended physics
    // and WILL FAIL until the exponent in integrator_utils.cuh is corrected to
    // (9/16)^(1/3).
    {
        const double a_0 = 1.0, R = 1.0;
        const double delta_N_0 = a_0 * a_0 / (3.0 * R);
        const double a_crit    = std::pow(1.0 / 6.0, 2.0 / 3.0) * a_0;
        const double delta_N_crit = -std::pow(9.0 / 16.0, 1.0 / 3.0) * delta_N_0;

        // At and below the pull-off displacement the radius is clamped to a_crit.
        CHECK_APPROX(get_contact_radius(delta_N_crit,          a_0, R),  a_crit,  1e-14);
        CHECK_APPROX(get_contact_radius(1.0001 * delta_N_crit, a_0, R),  a_crit,  1e-14);

        // Just ABOVE pull-off the contact still exists: x=0.4 maps to
        // y = -0.7849, between the true threshold (-0.8255) and the buggy one
        // (-0.6814). The solver must return a = 0.4*a_0, NOT the clamp value.
        {
            const double x_target = 0.4;
            const double delta_N = (3.0 * x_target * x_target - 2.0 * std::sqrt(x_target)) * delta_N_0;
            CHECK_APPROX(get_contact_radius(delta_N, a_0, R),  x_target * a_0,  1e-9);
        }

        // Continuity at pull-off: approaching the threshold from above the solver
        // value tends to a_crit with no jump. x=0.31 is just above x_min=0.3029.
        {
            const double x_near = 0.31;
            const double delta_N = (3.0 * x_near * x_near - 2.0 * std::sqrt(x_near)) * delta_N_0;
            CHECK_APPROX(get_contact_radius(delta_N, a_0, R),  x_near * a_0,  1e-6);
        }
    }

    // Far below critical: still clamped to a_crit (no solution exists).
    {
        CHECK_APPROX(get_contact_radius(-1e10, 1.0, 1.0),
                     std::pow(1.0 / 6.0, 2.0 / 3.0),  1e-14);
    }

    // R-invariance: for a fixed y the result is independent of R (R only enters
    // through delta_N_0). Hold y=2 and vary R, adjusting delta_N accordingly.
    {
        const double a_0 = 1.0, y = 2.0;
        const double a_R1 = get_contact_radius(y * a_0 * a_0 / (3.0 * 1.0), a_0, 1.0);
        const double a_R3 = get_contact_radius(y * a_0 * a_0 / (3.0 * 3.0), a_0, 3.0);
        CHECK_APPROX(a_R1, a_R3,  1e-12);
    }

    // a_0 linearity: for a fixed y the result scales linearly in a_0.
    {
        const double y = 2.0;
        const double a1 = get_contact_radius(y * 1.0 * 1.0 / (3.0 * 1.0), 1.0, 1.0);
        const double a2 = get_contact_radius(y * 2.0 * 2.0 / (3.0 * 1.0), 2.0, 1.0);
        CHECK_APPROX(a2, 2.0 * a1,  1e-12);
    }

    // Generic correctness: for several valid displacements (y above the clamp
    // threshold), the returned a must satisfy the original equation to within
    // solver tolerance. This validates convergence without precomputed roots.
    {
        const double a_0 = 1.0, R = 1.0;
        const double delta_N_0 = a_0 * a_0 / (3.0 * R);
        const double ys[] = {5.0, 2.0, 1.0, 0.0, -0.5};
        for (double y : ys) {
            const double a = get_contact_radius(y * delta_N_0, a_0, R);
            const double x = a / a_0;
            const double residual = 3.0 * x * x - 2.0 * std::sqrt(x) - y;
            CHECK_APPROX(residual, 0.0,  1e-9);
        }
    }
}
