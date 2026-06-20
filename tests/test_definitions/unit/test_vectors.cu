/**
 * @file test_vectors.cu
 * @brief Tests for utils/vector.cuh
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 */

#include "utils/vector.cuh"
#include <cmath>

void test_vectors() {

    // ------------------------------------------------------------------ //
    // Vector length
    // ------------------------------------------------------------------ //

    // 3-4-5 Pythagorean triple
    CHECK_APPROX(vec_length({3.0, 4.0, 0.0}), 5.0, 1e-14);
    // Unit vector
    CHECK_APPROX(vec_length({1.0, 0.0, 0.0}), 1.0, 1e-14);

    // ------------------------------------------------------------------ //
    // Squared vector length
    // ------------------------------------------------------------------ //

    CHECK_APPROX(vec_length_sq({3.0, 4.0, 0.0}), 25.0, 1e-14);
    CHECK_APPROX(vec_length_sq({0.0, 0.0, 5.0}), 25.0, 1e-14);

    // ------------------------------------------------------------------ //
    // Dot product
    // ------------------------------------------------------------------ //

    // Orthogonal basis vectors
    CHECK_APPROX(vec_dot({3.0, 0.0, 0.0}, {0.0, 2.0, 0.0}),  0.0, 1e-14);
    // Parallel
    CHECK_APPROX(vec_dot({2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}),  6.0, 1e-14);
    // Anti-parallel
    CHECK_APPROX(vec_dot({3.0, 0.0, 0.0}, {-2.0, 0.0, 0.0}), -6.0, 1e-14);
    // General: {2,3,4} . {1,1,1} = 9
    CHECK_APPROX(vec_dot({2.0, 3.0, 4.0}, {1.0, 1.0, 1.0}),  9.0, 1e-14);

    // ------------------------------------------------------------------ //
    // Cross product
    // ------------------------------------------------------------------ //

    // x × y = z  (right-hand rule)
    double3 c1 = vec_cross({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    CHECK_APPROX(c1.x, 0.0, 1e-14);
    CHECK_APPROX(c1.y, 0.0, 1e-14);
    CHECK_APPROX(c1.z, 1.0, 1e-14);

    // y × x = -z  (antisymmetry)
    double3 c2 = vec_cross({0.0, 1.0, 0.0}, {1.0, 0.0, 0.0});
    CHECK_APPROX(c2.z, -1.0, 1e-14);

    // u × u = 0
    double3 c3 = vec_cross({1.0, 2.0, 3.0}, {1.0, 2.0, 3.0});
    CHECK_APPROX(vec_length(c3), 0.0, 1e-14);

    // ------------------------------------------------------------------ //
    // Vector difference
    // ------------------------------------------------------------------ //

    double3 d = vec_diff({3.0, 2.0, 1.0}, {1.0, 1.0, 1.0});
    CHECK_APPROX(d.x, 2.0, 1e-14);
    CHECK_APPROX(d.y, 1.0, 1e-14);
    CHECK_APPROX(d.z, 0.0, 1e-14);

    // ------------------------------------------------------------------ //
    // Distance between two points
    // ------------------------------------------------------------------ //

    CHECK_APPROX(vec_dist_len({4.0, 0.0, 0.0}, {1.0, 0.0, 0.0}), 3.0, 1e-14);
    // 3-4-5 triangle from origin
    CHECK_APPROX(vec_dist_len({3.0, 4.0, 0.0}, {0.0, 0.0, 0.0}), 5.0, 1e-14);

    // ------------------------------------------------------------------ //
    // Squared distance
    // ------------------------------------------------------------------ //

    CHECK_APPROX(vec_dist_len_sq({4.0, 0.0, 0.0}, {1.0, 0.0, 0.0}), 9.0,  1e-14);
    CHECK_APPROX(vec_dist_len_sq({3.0, 4.0, 0.0}, {0.0, 0.0, 0.0}), 25.0, 1e-14);

    // ------------------------------------------------------------------ //
    // vec_normalize (in-place)
    // ------------------------------------------------------------------ //

    double3 n = {3.0, 4.0, 0.0};
    vec_normalize(n);
    CHECK_APPROX(vec_length(n),  1.0,       1e-14);
    CHECK_APPROX(n.x,            3.0 / 5.0, 1e-14);
    CHECK_APPROX(n.y,            4.0 / 5.0, 1e-14);

    // Zero vector: guard in the implementation leaves it unchanged
    double3 nz = {0.0, 0.0, 0.0};
    vec_normalize(nz);
    CHECK_APPROX(vec_length(nz), 0.0, 1e-14);

    // ------------------------------------------------------------------ //
    // vec_get_normal: unit vector pointing from v to u, i.e. (u-v)/|u-v|
    // ------------------------------------------------------------------ //

    double3 gn = vec_get_normal({4.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
    CHECK_APPROX(gn.x,            1.0, 1e-14);
    CHECK_APPROX(gn.y,            0.0, 1e-14);
    CHECK_APPROX(gn.z,            0.0, 1e-14);
    CHECK_APPROX(vec_length(gn),  1.0, 1e-14);

    // ------------------------------------------------------------------ //
    // vec_get_normalized: normalises the argument in-place AND returns it
    // ------------------------------------------------------------------ //

    double3 vgn = {0.0, 3.0, 4.0};
    double3 ret = vec_get_normalized(vgn);
    // The input is modified
    CHECK_APPROX(vec_length(vgn), 1.0, 1e-14);
    // The return value equals the modified input
    CHECK_APPROX(ret.x, vgn.x, 1e-15);
    CHECK_APPROX(ret.y, vgn.y, 1e-15);
    CHECK_APPROX(ret.z, vgn.z, 1e-15);

    // ------------------------------------------------------------------ //
    // vec_set
    // ------------------------------------------------------------------ //

    double3 vs;
    vec_set(vs, 7.0);
    CHECK_APPROX(vs.x, 7.0, 1e-15);
    CHECK_APPROX(vs.y, 7.0, 1e-15);
    CHECK_APPROX(vs.z, 7.0, 1e-15);

    // ------------------------------------------------------------------ //
    // vec_is_zero
    // NOTE: the implementation checks x+y+z == 0, not component-wise.
    //       Vectors like {1,-1,0} also return true (documented warning).
    //       Tests here only cover the unambiguous cases.
    // ------------------------------------------------------------------ //

    CHECK( vec_is_zero({0.0, 0.0,  0.0}));
    CHECK(!vec_is_zero({1.0, 0.0,  0.0}));
    CHECK(!vec_is_zero({1.0, 1.0, -1.0}));

    // ------------------------------------------------------------------ //
    // Quaternion length
    // ------------------------------------------------------------------ //

    // Identity quaternion {x,y,z,w} = {0,0,0,1}
    CHECK_APPROX(quat_length({0.0, 0.0, 0.0, 1.0}), 1.0, 1e-14);
    // 3-4-5 triple in the w and x components
    CHECK_APPROX(quat_length({3.0, 0.0, 0.0, 4.0}), 5.0, 1e-14);

    // ------------------------------------------------------------------ //
    // quat_normalize (in-place)
    // ------------------------------------------------------------------ //

    double4 q = {3.0, 0.0, 0.0, 4.0};
    quat_normalize(q);
    CHECK_APPROX(quat_length(q), 1.0,       1e-14);
    CHECK_APPROX(q.x,            3.0 / 5.0, 1e-14);
    CHECK_APPROX(q.w,            4.0 / 5.0, 1e-14);

    // ------------------------------------------------------------------ //
    // quat_apply: identity quaternion leaves vector unchanged
    // ------------------------------------------------------------------ //

    double4 q_id = {0.0, 0.0, 0.0, 1.0};
    double3 v    = {1.0, 2.0, 3.0};
    double3 r    = quat_apply(q_id, v);
    CHECK_APPROX(r.x, 1.0, 1e-13);
    CHECK_APPROX(r.y, 2.0, 1e-13);
    CHECK_APPROX(r.z, 3.0, 1e-13);

    // ------------------------------------------------------------------ //
    // quat_apply: 180-degree rotation around z
    // q.w = cos(π/2) = 0, q.z = sin(π/2) = 1  =>  {1,0,0} -> {-1,0,0}
    // ------------------------------------------------------------------ //

    double4 q_180z = {0.0, 0.0, 1.0, 0.0};
    double3 r2     = quat_apply(q_180z, {1.0, 0.0, 0.0});
    CHECK_APPROX(r2.x, -1.0, 1e-13);
    CHECK_APPROX(r2.y,  0.0, 1e-13);
    CHECK_APPROX(r2.z,  0.0, 1e-13);

    // ------------------------------------------------------------------ //
    // quat_apply_inverse: round-trip restores original vector
    // ------------------------------------------------------------------ //

    double4 q_rot = {0.0, 0.0, sin(0.4), cos(0.4)};  // rotation around z
    double3 v2    = {1.5, -2.0, 0.7};
    double3 rt    = quat_apply_inverse(q_rot, quat_apply(q_rot, v2));
    CHECK_APPROX(rt.x, v2.x, 1e-12);
    CHECK_APPROX(rt.y, v2.y, 1e-12);
    CHECK_APPROX(rt.z, v2.z, 1e-12);
}
