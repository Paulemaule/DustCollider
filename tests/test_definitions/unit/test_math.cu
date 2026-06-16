/**
 * @file test_math.cu
 * @brief Tests for the __host__ __device__ math helpers in physics/integrator_utils.cuh
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 *
 * NOTE: get_contact_radius is __device__ only and is therefore not tested here.
 * To test it on CPU, add __host__ to its declaration in integrator_utils.cuh.
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

    // ------------------------------------------------------------------ //
    // get_G_i: shear modulus  E / (2 * (1 + nu))
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_G_i(2.0, 0.0),  1.0,  1e-14);  // 2/(2*1) = 1
    CHECK_APPROX(get_G_i(3.0, 0.5),  1.0,  1e-14);  // 3/(2*1.5) = 1

    // ------------------------------------------------------------------ //
    // get_E_s: combined Young's modulus
    //   1 / ((1-nu_i^2)/E_i + (1-nu_j^2)/E_j)
    //   equal materials, nu=0 → E/2
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_E_s(4.0, 4.0, 0.0, 0.0),  2.0,  1e-14);

    // ------------------------------------------------------------------ //
    // get_G_s: combined shear modulus (same formula as get_E_s)
    //   equal materials, nu=0 → G/2
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_G_s(4.0, 4.0, 0.0, 0.0),  2.0,  1e-14);

    // ------------------------------------------------------------------ //
    // get_gamma: combined surface energy
    //   gamma_i + gamma_j - 2 / (1/gamma_i + 1/gamma_j)
    //   equal materials: g + g - 2/(2/g) = g
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_gamma(1.0, 1.0),  1.0,  1e-14);
    CHECK_APPROX(get_gamma(2.0, 2.0),  2.0,  1e-14);
    // 1 + 3 - 2*1*3/(1+3) = 4 - 1.5 = 2.5
    CHECK_APPROX(get_gamma(1.0, 3.0),  2.5,  1e-14);

    // ------------------------------------------------------------------ //
    // get_a_0: equilibrium contact radius  (9*PI*gamma*R^2 / E_s)^(1/3)
    //   gamma=1, R=1, E_s=9*PI  →  a_0 = 1
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_a_0(1.0, 1.0, 9.0 * PI),  1.0,  1e-12);

    // ------------------------------------------------------------------ //
    // get_delta_N_crit: critical normal displacement
    //   0.5 * a_0^2 / (R * cbrt(6))
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_delta_N_crit(1.0, 1.0),  0.5 / std::cbrt(6.0),  1e-14);
    // a_0=2, R=1 → 0.5*4 / cbrt(6) = 2/cbrt(6)
    CHECK_APPROX(get_delta_N_crit(2.0, 1.0),  2.0 / std::cbrt(6.0),  1e-14);

    // ------------------------------------------------------------------ //
    // get_delta_S_crit: critical sliding displacement
    //   (2 - 0.5*(nu_i+nu_j)) * a_0 / (16*PI)
    //   nu=0, a_0=16*PI → result = 2
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_delta_S_crit(0.0, 0.0, 16.0 * PI),  2.0,  1e-12);

    // ------------------------------------------------------------------ //
    // get_U_S and get_U_R: quadratic sliding / rolling potential
    //   0.5 * k * |disp|^2
    //   k=2, {3,4,0} → 0.5*2*25 = 25
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_U_S(2.0, {3.0, 4.0, 0.0}),  25.0,  1e-14);
    CHECK_APPROX(get_U_R(2.0, {3.0, 4.0, 0.0}),  25.0,  1e-14);

    // ------------------------------------------------------------------ //
    // get_U_T: quadratic twisting potential  0.5 * k * twist^2
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_U_T(2.0, 3.0),  9.0,  1e-14);

    // ------------------------------------------------------------------ //
    // get_U_N: normal (JKR) potential
    //   At a=0 the polynomial terms (a/a_0)^n all vanish, leaving only the
    //   constant prefactor: F_c * delta_N_crit * 0.84661389438303971
    // ------------------------------------------------------------------ //

    CHECK_APPROX(get_U_N(1.0, 1.0, 0.0, 1.0),  0.84661389438303971,  1e-14);
    // Scales linearly with F_c and delta_N_crit
    CHECK_APPROX(get_U_N(2.0, 3.0, 0.0, 1.0),  6.0 * 0.84661389438303971,  1e-14);

    // ------------------------------------------------------------------ //
    // get_contact_radius: Newton solver for JKR contact radius
    //
    // The equation solved is  3x^2 - 2*sqrt(x) - y = 0
    // where x = a/a_0 and y = delta_N / delta_N_0.
    //
    // At y=1 (delta_N = delta_N_0), x=1 is an exact root:
    //   3*1 - 2*1 - 1 = 0.  Newton starts there and stays → a = a_0.
    //
    // At y=0 (delta_N = 0), the root is x = (2/3)^(2/3) (JKR adhesion
    //   contact at zero external force).
    //
    // Below the critical displacement, the function clamps to
    //   a_crit = (1/6)^(2/3) * a_0.
    // ------------------------------------------------------------------ //

    // Equilibrium: delta_N = delta_N_0 = a_0^2/(3R) → a = a_0
    {
        const double a_0 = 1.0, R = 1.0;
        const double delta_N_0 = a_0 * a_0 / (3.0 * R);
        CHECK_APPROX(get_contact_radius(delta_N_0, a_0, R),  a_0,  1e-14);
    }

    // Zero displacement: a = (2/3)^(2/3) * a_0
    {
        CHECK_APPROX(get_contact_radius(0.0, 1.0, 1.0),
                     std::pow(2.0 / 3.0, 2.0 / 3.0),  1e-10);
    }

    // Below critical displacement: clamped to (1/6)^(2/3) * a_0
    {
        CHECK_APPROX(get_contact_radius(-1e10, 1.0, 1.0),
                     std::pow(1.0 / 6.0, 2.0 / 3.0),  1e-14);
    }
}
