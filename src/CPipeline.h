#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include<fstream>
using namespace std;

#include "typedefs.h"
//

#ifndef CPIPELINE
#define CPIPELINE

class CPipeline
{
public:
    CPipeline();
    ~CPipeline() {}

    bool init(int argc, const char** argv);

    void formatLine(string& line);
    string seperateString(string& str);
    dlist parseValues(string& str);
    bool parse();
    bool checkParameters();
    bool parseLine(string cmd, string data);
    bool prepareData(vec3D*& pos, vec3D*& vel, vec3D*& omega_tot, vec3D*& mag, double*& amon, double*& mass, double*& moment, int*& matIDs, int& Nmon);
    void prepareMaterial(material*& mat, int& Nmat);
    bool isMatID(int id);
    

    bool writeOVITO(ullong time_id, double b_size, const vec3D* pos, const double* amon, const int* matID, int Nmon);
    bool writeAllOVITO(const vec3D* pos, const vec3D* vel, const vec3D* force, const vec3D* torque, const vec3D* omega, const vec3D* mag, const int* cluserIDs, const double* amon, const int* matID, int Nmon, int Nsto);


    void printParameters();

    bool savePos() { return save_pos; }
    bool saveVel() { return save_vel; }
    bool saveForce() { return save_force; }
    bool saveTorque() { return save_torque; }
    bool saveCluster() { return save_cluster; }
    bool saveOvito() { return save_ovito; }
    
    bool saveOmega() { return save_omega; }
    bool saveMag() { return save_mag; }
    
    vec3D getBext() { return Bext; }

    bool writeHeader();
    bool writeMonomers();

    bool writeBinaryVec(string name_file, const vec3D* data, ullong N);
    bool writeBinaryDouble(string name_file, const double* data, ullong N);
    bool writeBinaryInt(string name_file, const int* data, ullong N);


    //getter and setter
    ullong getNIter() { return N_iter; }
    ullong getNSave() { return N_save; }

    double getTimeStep() { return time_step; };
private:
    bool createPath(string path);
    double getMin(dlist& lst);
    double getMax(dlist & lst);

    string cmd_filename;
    string path_results;

    string path_binary;
    string path_ovito;
    string path_plots;

    string path_A;
    string path_B;

    vec3D pos_A;
    vec3D pos_B;

    vec3D vel_A;
    vec3D vel_B;

    vec3D ang_A;
    vec3D ang_B;

    ullong N_iter;
    ullong N_save;

    double time_start;
    double time_stop;
    double time_step;

    strlist lst_matName;
    ilist lst_matID;  //material id
    dlist lst_gamma;  //surface energy
    dlist lst_E;      //Young's modulus
    dlist lst_nu;     //Poisson number 
    dlist lst_rho;    //density
    dlist lst_xi;     //critical rolling length
    dlist lst_Tvis;   //viscous dumping time 

    dlist lst_tss;     //spin-spin time scale
    dlist lst_tsl;     //spin-lattice time scale
    dlist lst_Msat;    //spont. magnetization
    dlist lst_chi;     //mag. susceptibillity
    dlist lst_Tc;      //Curie temperature

    int mat_type;

    double min_gamma;
    double min_E;
    double min_nu;
    double min_rho;
    double min_xi;
    double min_Tvis;

    double min_tss;
    double min_tsl;
    double min_Msat;
    double min_chi;
    double min_Tc;

    double max_gamma;
    double max_E;
    double max_nu;
    double max_rho;
    double max_xi;
    double max_Tvis;

    double max_tss;
    double max_tsl;
    double max_Msat;
    double max_chi;
    double max_Tc;

    double a_eff_A;
    double a_eff_B;

    double a_out_A;
    double a_out_B;

    double a_mon_min;
    double a_mon_max;

    int Nmon_A;
    int Nmon_B;

    string agg_filename_A;
    string agg_filename_B;

    bool save_pos;
    bool save_vel;
    bool save_force;
    bool save_torque;
    bool save_cluster;
    bool save_omega;
    bool save_mag;
    bool save_ovito;

    vec3D Bext;
    double T_dust;
};

#endif


