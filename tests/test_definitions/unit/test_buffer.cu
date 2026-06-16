/**
 * @file test_buffer.cu
 * @brief Tests for the RAII GPU memory wrappers: DeviceBuffer/HostBuffer
 * (utils/buffer.cuh), DeviceState/HostState (physics/state.cuh) and
 * DeviceMaterials/HostMaterials (physics/materials.cuh).
 *
 * These tests perform real CUDA allocations and host<->device copies on
 * whatever GPU is available in the test environment.
 *
 * This file is #included into test_main.cu; testkit.h is expected to be
 * already in scope.
 */

#include <utility>
#include "utils/buffer.cuh"
#include "physics/state.cuh"
#include "physics/materials.cuh"

void test_buffer() {

    // ------------------------------------------------------------------ //
    // HostBuffer: allocation, size accessors, and direct read/write
    // (pinned host memory is directly addressable by host code).
    // ------------------------------------------------------------------ //

    {
        HostBuffer<double> hb(5);
        CHECK(hb.size() == 5);
        CHECK(hb.size_bytes() == 5 * sizeof(double));
        CHECK(hb.data() != nullptr);

        for (int i = 0; i < 5; i++) hb.data()[i] = i * 2.0;
        for (int i = 0; i < 5; i++) CHECK_APPROX(hb.data()[i], i * 2.0, 1e-14);
    }

    // ------------------------------------------------------------------ //
    // HostBuffer: move constructor transfers ownership and data
    // ------------------------------------------------------------------ //

    {
        HostBuffer<double> hb(3);
        for (int i = 0; i < 3; i++) hb.data()[i] = i + 1.0;

        HostBuffer<double> hb2(std::move(hb));
        CHECK(hb2.size() == 3);
        for (int i = 0; i < 3; i++) CHECK_APPROX(hb2.data()[i], i + 1.0, 1e-14);
    }

    // ------------------------------------------------------------------ //
    // HostBuffer: move assignment frees the old allocation and transfers
    // ownership and data from the source.
    // ------------------------------------------------------------------ //

    {
        HostBuffer<double> hb(4);
        for (int i = 0; i < 4; i++) hb.data()[i] = i + 10.0;

        HostBuffer<double> hb2(2);
        hb2 = std::move(hb);
        CHECK(hb2.size() == 4);
        for (int i = 0; i < 4; i++) CHECK_APPROX(hb2.data()[i], i + 10.0, 1e-14);
    }

    // ------------------------------------------------------------------ //
    // HostBuffer: zero-element allocation must not crash on construction
    // or destruction.
    // ------------------------------------------------------------------ //

    {
        HostBuffer<double> hb(0);
        CHECK(hb.size() == 0);
        CHECK(hb.size_bytes() == 0);
    }

    // ------------------------------------------------------------------ //
    // DeviceBuffer: allocation and size accessors. Device memory is not
    // host-addressable, so correctness of the allocation is verified via
    // a host -> device -> host round trip through cudaMemcpy.
    // ------------------------------------------------------------------ //

    {
        DeviceBuffer<double> db(5);
        CHECK(db.size() == 5);
        CHECK(db.size_bytes() == 5 * sizeof(double));
        CHECK(db.data() != nullptr);

        double src[5];
        for (int i = 0; i < 5; i++) src[i] = i * 3.0;
        cudaMemcpy(db.data(), src, sizeof(src), cudaMemcpyHostToDevice);

        double dst[5] = {};
        cudaMemcpy(dst, db.data(), sizeof(dst), cudaMemcpyDeviceToHost);
        for (int i = 0; i < 5; i++) CHECK_APPROX(dst[i], i * 3.0, 1e-14);
    }

    // ------------------------------------------------------------------ //
    // DeviceBuffer: move constructor transfers ownership; the data
    // previously written into the source buffer remains accessible
    // through the moved-to buffer's pointer.
    // ------------------------------------------------------------------ //

    {
        DeviceBuffer<double> db(3);
        double src[3] = {1.0, 2.0, 3.0};
        cudaMemcpy(db.data(), src, sizeof(src), cudaMemcpyHostToDevice);

        DeviceBuffer<double> db2(std::move(db));
        CHECK(db2.size() == 3);

        double dst[3] = {};
        cudaMemcpy(dst, db2.data(), sizeof(dst), cudaMemcpyDeviceToHost);
        for (int i = 0; i < 3; i++) CHECK_APPROX(dst[i], src[i], 1e-14);
    }

    // ------------------------------------------------------------------ //
    // DeviceBuffer: move assignment frees the old allocation and
    // transfers ownership and data from the source.
    // ------------------------------------------------------------------ //

    {
        DeviceBuffer<double> db(4);
        double src[4] = {1.0, 2.0, 3.0, 4.0};
        cudaMemcpy(db.data(), src, sizeof(src), cudaMemcpyHostToDevice);

        DeviceBuffer<double> db2(2);
        db2 = std::move(db);
        CHECK(db2.size() == 4);

        double dst[4] = {};
        cudaMemcpy(dst, db2.data(), sizeof(dst), cudaMemcpyDeviceToHost);
        for (int i = 0; i < 4; i++) CHECK_APPROX(dst[i], src[i], 1e-14);
    }

    // ------------------------------------------------------------------ //
    // DeviceBuffer: zero-element allocation must not crash on
    // construction or destruction.
    // ------------------------------------------------------------------ //

    {
        DeviceBuffer<double> db(0);
        CHECK(db.size() == 0);
        CHECK(db.size_bytes() == 0);
    }

    // ------------------------------------------------------------------ //
    // HostState / DeviceState: push_to() followed by pull_from() into a
    // second HostState round-trips every array, including the Nmon x Nmon
    // contact matrices.
    // ------------------------------------------------------------------ //

    {
        const size_t Nmon = 3;

        HostState hs(Nmon);
        HostStateView v = hs.view();
        for (size_t i = 0; i < Nmon; i++) {
            v.position[i] = {double(i), double(i) + 0.1, double(i) + 0.2};
            v.velocity[i] = {double(i) * 2.0, 0.0, 0.0};
            v.omega[i]    = {0.0, double(i) * 3.0, 0.0};
            v.force[i]    = {0.0, 0.0, double(i) * 4.0};
            v.torque[i]   = {double(i) + 5.0, 0.0, 0.0};
        }
        for (size_t i = 0; i < Nmon * Nmon; i++) {
            v.contact_compression[i] = double(i) * 0.5;
            v.contact_twist[i]       = double(i) * 0.25;
            v.contact_pointer[i]     = {double(i), 0.0, 0.0};
            v.contact_normal[i]      = {0.0, double(i), 0.0};
            v.contact_rotation[i]    = {double(i), 0.0, 0.0, 1.0};
        }

        DeviceState ds(Nmon);
        hs.push_to(ds, Nmon);

        HostState hs2(Nmon);
        hs2.pull_from(ds, Nmon);
        HostStateView v2 = hs2.view();

        for (size_t i = 0; i < Nmon; i++) {
            CHECK_APPROX(v2.position[i].x, v.position[i].x, 1e-14);
            CHECK_APPROX(v2.position[i].y, v.position[i].y, 1e-14);
            CHECK_APPROX(v2.position[i].z, v.position[i].z, 1e-14);
            CHECK_APPROX(v2.velocity[i].x, v.velocity[i].x, 1e-14);
            CHECK_APPROX(v2.omega[i].y,    v.omega[i].y,    1e-14);
            CHECK_APPROX(v2.force[i].z,    v.force[i].z,    1e-14);
            CHECK_APPROX(v2.torque[i].x,   v.torque[i].x,   1e-14);
        }
        for (size_t i = 0; i < Nmon * Nmon; i++) {
            CHECK_APPROX(v2.contact_compression[i], v.contact_compression[i], 1e-14);
            CHECK_APPROX(v2.contact_twist[i],       v.contact_twist[i],       1e-14);
            CHECK_APPROX(v2.contact_pointer[i].x,   v.contact_pointer[i].x,   1e-14);
            CHECK_APPROX(v2.contact_normal[i].y,    v.contact_normal[i].y,    1e-14);
            CHECK_APPROX(v2.contact_rotation[i].w,  v.contact_rotation[i].w,  1e-14);
        }
    }

    // ------------------------------------------------------------------ //
    // HostMaterials / DeviceMaterials: push_to() followed by pull_from()
    // into a second HostMaterials round-trips every array.
    // ------------------------------------------------------------------ //

    {
        const size_t Nmon = 4;

        HostMaterials hm(Nmon);
        HostMaterialsView v = hm.view();
        for (size_t i = 0; i < Nmon; i++) {
            v.radius[i]            = 1e-9 * (i + 1);
            v.mass[i]               = 1e-18 * (i + 1);
            v.moment[i]              = 1e-27 * (i + 1);
            v.matID[i]               = int(i);
            v.density[i]             = 3000.0 + i;
            v.surface_energy[i]      = 0.01 * (i + 1);
            v.youngs_modulus[i]      = 1e11 * (i + 1);
            v.poisson_number[i]      = 0.2 + i * 0.01;
            v.damping_timescale[i]   = 1e-10 * (i + 1);
            v.crit_rolling_disp[i]   = 1e-11 * (i + 1);
        }

        DeviceMaterials dm(Nmon);
        hm.push_to(dm, Nmon);

        HostMaterials hm2(Nmon);
        hm2.pull_from(dm, Nmon);
        HostMaterialsView v2 = hm2.view();

        for (size_t i = 0; i < Nmon; i++) {
            CHECK_APPROX(v2.radius[i],            v.radius[i],            1e-20);
            CHECK_APPROX(v2.mass[i],               v.mass[i],               1e-25);
            CHECK_APPROX(v2.moment[i],              v.moment[i],              1e-34);
            CHECK(v2.matID[i] == v.matID[i]);
            CHECK_APPROX(v2.density[i],             v.density[i],             1e-10);
            CHECK_APPROX(v2.surface_energy[i],      v.surface_energy[i],      1e-14);
            CHECK_APPROX(v2.youngs_modulus[i],      v.youngs_modulus[i],      1e-6);
            CHECK_APPROX(v2.poisson_number[i],      v.poisson_number[i],      1e-14);
            CHECK_APPROX(v2.damping_timescale[i],   v.damping_timescale[i],   1e-20);
            CHECK_APPROX(v2.crit_rolling_disp[i],   v.crit_rolling_disp[i],   1e-20);
        }
    }
}
