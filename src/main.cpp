#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <string>
//#include <studio.h>
#include <iostream>
#include<fstream>
#include <algorithm>

using namespace std;
using namespace std::chrono;

#include "typedefs.h"
#include "vector.h"
#include "CPipeline.h"
#include "physics.h"

int main(const int argc, const char** argv)
{
    CPipeline pipeline;
    auto start = high_resolution_clock::now();

    if (!pipeline.init(argc, argv))
        return -1;

    if (!pipeline.parse())
        return -1;

    if (!pipeline.checkParameters())
        return -1;
    
    material* mat=0;
    int Nmat = 0;

    pipeline.prepareMaterial(mat, Nmat);

    vec3D* vel = 0;
    vec3D* omega = 0;     //angular velocity due to monomer spin
    vec3D* omega_tot = 0; //angular velocity due to monomer spin + curved trajectory
    vec3D* mag = 0;

    vec3D* pos_old = 0;
    vec3D* pos_new = 0;
    vec3D* force_old = 0;
    vec3D* force_new = 0;
    vec3D* torque_old = 0;
    vec3D* torque_new = 0;

    vec3D* dMdt_old = 0;
    vec3D* dMdt_new = 0;


    vec3D* storage_pos = 0;
    vec3D* storage_vel = 0;
    vec3D* storage_force = 0;
    vec3D* storage_torque = 0;
    vec3D* storage_omega = 0;
    vec3D* storage_mag = 0;
    int * storage_cluster = 0;

    vec3D* matrix_con = 0; //contact pointer between monomers
    vec3D* matrix_norm = 0; //normal vectors between monomers
    quat* matrix_rot = 0; //contact pointer rotation direction
    double* matrix_comp = 0; //old compression lengths, also used to track connection
    double* matrix_twist = 0; //old twisting displacement

    int* matIDs = 0;
    double* amon = 0;
    double* moment = 0;
    double* mass = 0;
    int Nmon = 0;
    int* clusterIDs = 0;

    pipeline.prepareData(pos_old, vel, omega_tot, mag, amon, mass, moment, matIDs, Nmon);
    pipeline.printParameters();

    omega = new vec3D[Nmon];
    
    pos_new = new vec3D[Nmon];

    force_old = new vec3D[Nmon];
    force_new = new vec3D[Nmon];

    torque_old = new vec3D[Nmon];
    torque_new = new vec3D[Nmon];

    dMdt_old = new vec3D[Nmon];
    dMdt_new = new vec3D[Nmon];

    memset(omega, 0, Nmon * sizeof(vec3D));
  
    memset(pos_new, 0, Nmon * sizeof(vec3D));

    memset(force_old, 0, Nmon * sizeof(vec3D));
    memset(force_new, 0, Nmon * sizeof(vec3D));

    memset(torque_old, 0, Nmon * sizeof(vec3D));
    memset(torque_new, 0, Nmon * sizeof(vec3D));

    memset(dMdt_old, 0, Nmon * sizeof(vec3D));
    memset(dMdt_new, 0, Nmon * sizeof(vec3D));

    int Ncon = Nmon * Nmon;

    matrix_con = new vec3D[2 * Ncon];
    matrix_rot = new quat[2 * Ncon];

    matrix_norm = new vec3D[Ncon];

    matrix_comp = new double[Ncon];
    matrix_twist = new double[Ncon];

    memset(matrix_con, 0, 2 * Ncon * sizeof(vec3D));
    memset(matrix_rot, 0, 2 * Ncon * sizeof(quat));

    memset(matrix_norm, 0, Ncon * sizeof(vec3D));

    fill(matrix_comp, matrix_comp + Ncon, -1.0);
    //memset(matrix_comp, -1, Ncon * sizeof(double));
    memset(matrix_twist, 0, Ncon * sizeof(double));

    ullong N_iter = pipeline.getNIter();
    ullong N_save = pipeline.getNSave();
    double time_step = pipeline.getTimeStep();
    bool save_ovito = pipeline.saveOvito();
    vec3D B_ext = pipeline.getBext();

    int N_store = 0;

    if (N_save > 0)
    {
        N_store = Nmon * int((double(N_iter) / double(N_save)+0.5));

        clusterIDs = new int[Nmon];
        fill(clusterIDs, clusterIDs + Nmon, -1);

        if (pipeline.savePos())
        {
            storage_pos = new vec3D[N_store];
            memset(storage_pos, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveVel())
        {
            storage_vel = new vec3D[N_store];
            memset(storage_vel, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveForce())
        {
            storage_force = new vec3D[N_store];
            memset(storage_force, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveTorque())
        {
            storage_torque = new vec3D[N_store];
            memset(storage_torque, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveOmega())
        {
            storage_omega = new vec3D[N_store];
            memset(storage_omega, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveMag())
        {
            storage_mag = new vec3D[N_store];
            memset(storage_mag, 0, N_store * sizeof(vec3D));
        }

        if (pipeline.saveCluster())
        {
            storage_cluster = new int[N_store];
            memset(storage_cluster, 0, N_store * sizeof(int));
        }
    }

    int nBlocks = (Nmon + BLOCK_SIZE + 1) / BLOCK_SIZE;

    ullong counter_save = 0;

    for (ullong iter = 0; iter < N_iter; iter++)
    {
        predictor(pos_old, pos_new, force_old, vel, mass, time_step, Nmon);

        updateNeighbourhoodRelations(pos_new, matrix_con, matrix_norm, matrix_rot, matrix_comp, matrix_twist, amon, mat, matIDs, Nmon);
        
        // todo: rotate magnetization here
        updateContacts(omega, omega_tot, torque_old, mag, matrix_rot, matrix_comp, moment, Nmon, time_step);

        updateParticleInteraction(pos_new, force_new, torque_old, torque_new, dMdt_new, matrix_con, matrix_norm, omega, omega_tot,mag, matrix_rot, matrix_comp, matrix_twist, amon, moment, mat, matIDs, B_ext, Nmon, time_step);
        corrector(force_old, force_new, torque_old, torque_new, dMdt_old, dMdt_new, vel, omega, omega_tot, mag, mass, moment, mat, matIDs, time_step, Nmon);

        switch_pointer(pos_old, pos_new, force_old, force_new, torque_old, torque_new, dMdt_old, dMdt_new);

        if (N_save > 0)
        {
            if (iter % N_save == 0)
            {
                ullong start_index = counter_save * Nmon;

                if (long(N_store) - long(start_index) < Nmon)
                {
                    cout << "ERROR: Stirage overrun!  \n";
                    cout << long(N_store) - long(start_index) << "\t" << Nmon << endl;
                    return -1;
                }

                if (start_index < N_store) //just a failesave if there was a rounding error in N_store
                {
                    if (storage_pos != 0)
                        copy(pos_old, pos_old + Nmon, storage_pos + start_index);

                    if (storage_vel != 0)
                        copy(vel, vel + Nmon, storage_vel + start_index);

                    if (storage_force != 0)
                        copy(force_old, force_old + Nmon, storage_force + start_index);

                    if (storage_torque != 0)
                        copy(torque_old, torque_old + Nmon, storage_torque + start_index);

                    if (storage_omega != 0)
                        copy(omega_tot, omega_tot + Nmon, storage_omega + start_index);

                    if (storage_mag != 0)
                        copy(mag, mag + Nmon, storage_mag + start_index);

                    if (storage_cluster != 0)
                    {
                        fill(clusterIDs, clusterIDs + Nmon, -1);
                        findConnectedComponents(Nmon, matrix_comp, clusterIDs);
                        copy(clusterIDs, clusterIDs + Nmon, storage_cluster + start_index);
                    }
                }
                counter_save++;
            }
        }

        if (iter % 2000 == 0)
        {
            printf("progress: %.4f %%      \r", 100.0 * float(iter) / N_iter);
        }
    }

    cout << CLR_LINE;

    if (N_save > 0)
    {
        cout << SEP_LINE;

        cout << "Writing simulation data:    \n" << flush;

        // write ovito files
        if (save_ovito && storage_pos != 0)
        {
            if (!pipeline.writeAllOVITO(storage_pos, storage_vel, storage_force, storage_torque, storage_omega, storage_mag, storage_cluster, amon, matIDs, Nmon, N_store))
                return -1;
        }

        //write header
        if (!pipeline.writeHeader())
            return -1;

        //write aggregate parameters
        if (!pipeline.writeBinaryDouble("agg_a_mon.bin", amon, Nmon))
            return -1;

        if (!pipeline.writeBinaryDouble("agg_mass_mon.bin", mass, Nmon))
            return -1;

        for (int i = 0; i < Nmon; i++)
            matIDs[i]= matIDs[i]+1;

        if (!pipeline.writeBinaryInt("agg_matid_mon.bin", matIDs, Nmon))
            return -1;

        //write material parameters
        //if (!pipeline.writeBinaryDouble("mat_gamma.bin", mat_gamma, Nmat))
        //    return -1;

        //if (!pipeline.writeBinaryDouble("mat_E.bin", mat_E, Nmat))
        //    return -1;

        //if (!pipeline.writeBinaryDouble("mat_nu.bin", mat_nu, Nmat))
        //    return -1;

        //if (!pipeline.writeBinaryDouble("mat_rho.bin", mat_rho, Nmat))
        //    return -1;

        //if (!pipeline.writeBinaryDouble("mat_xi.bin", mat_xi, Nmat))
         //   return -1;

        //if (!pipeline.writeBinaryDouble("mat_Tvis.bin", mat_Tvis, Nmat))
        //    return -1;

        // write simulation data
        if (storage_pos != 0)
        {
            if (!pipeline.writeBinaryVec("sim_pos.bin", storage_pos, N_store))
                return -1;
        }

        if (storage_vel != 0)
            if (!pipeline.writeBinaryVec("sim_vel.bin", storage_vel, N_store))
                return -1;

        if (storage_force != 0)
            if (!pipeline.writeBinaryVec("sim_force.bin", storage_force, N_store))
                return -1;

        if (storage_torque != 0)
            if (!pipeline.writeBinaryVec("sim_torque.bin", storage_torque, N_store))
                return -1;

        if(storage_cluster != 0)
            if (!pipeline.writeBinaryInt("sim_cluster.bin", storage_cluster, N_store))
                return -1;
    }

    cout << "-> Final clanup ...              \r" << flush;

    if (omega != 0)
        delete[] omega;

    if (omega_tot != 0)
        delete[] omega_tot;

    if (vel != 0)
        delete[] vel;

    if (mag != 0)
        delete[] mag;

    if (pos_old != 0)
        delete[] pos_old;

    if (pos_new != 0)
        delete[] pos_new;

    if (matIDs != 0)
        delete[] matIDs;

    if (amon != 0)
        delete[] amon;

    if (moment != 0)
        delete[] moment;

    if (matrix_con != 0)
        delete[] matrix_con;

    if (matrix_norm != 0)
        delete[] matrix_norm;

    if (matrix_rot != 0)
        delete[] matrix_rot;

    if (matrix_comp != 0)
        delete[] matrix_comp;

    if (matrix_twist != 0)
        delete[] matrix_twist;

    if (mass != 0)
        delete[] mass;

    if (force_old != 0)
        delete[] force_old;

    if (force_new != 0)
        delete[] force_new;

    if (torque_old != 0)
        delete[] torque_old;

    if (torque_new != 0)
        delete[] torque_new;

    if (dMdt_old != 0)
        delete[] dMdt_old;

    if (dMdt_new != 0)
        delete[] dMdt_new;

    if (mat != 0)
        delete[] mat;



    if (storage_pos != 0)
        delete[] storage_pos;

    if (storage_vel != 0)
        delete[] storage_vel;

    if (storage_force != 0)
        delete[] storage_force;

    if (storage_torque != 0)
        delete[] storage_torque;

    if (storage_omega != 0)
        delete[] storage_omega;

    if (storage_mag != 0)
        delete[] storage_mag;

    if (storage_cluster != 0)
        delete[] storage_cluster;


    cout << SEP_LINE;
    cout << "  - Final clanup: done              \n" << flush;
    cout << SEP_LINE;

    auto end = high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("Run time for %llu iterations : %.3f seconds.\n", N_iter, elapsed.count() * 1e-9);
    cout << SEP_LINE;

    //cout << "All done\n";
    return 0;
}