#include <sys/stat.h>
#include <sstream>
#include <math.h>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include "CPipeline.h"
#include "vector.h"

CPipeline::CPipeline()
{
    pos_A.x = 0;
    pos_A.y = 0;
    pos_A.z = 0;

    pos_B.x = 0;
    pos_B.y = 0;
    pos_B.z = 0;


    vel_A.x = 0;
    vel_A.y = 0;
    vel_A.z = 0;

    vel_B.x = 0;
    vel_B.y = 0;
    vel_B.z = 0;


    ang_A.x = 0;
    ang_A.y = 0;
    ang_A.z = 0;

    ang_B.x = 0;
    ang_B.y = 0;
    ang_B.z = 0;

    N_iter = 0;
    N_save = 0;

    time_start = 0;
    time_stop = 0;
    time_step = 0;

    a_eff_A = 0;
    a_eff_B = 0;

    Nmon_A = 0;
    Nmon_B = 0;

    a_out_A = 0;
    a_out_B = 0;

    a_mon_min = 0;
    a_mon_max = 0;

    save_pos = false;
    save_vel = false;
    save_force = false;
    save_torque = false;
    save_cluster = false;
    save_omega = false;
    save_mag = false;
    save_ovito = false;
    

    Bext.x = 0;  //T
    Bext.y = 0;
    Bext.z = 0;
    T_dust = 15; //K

    mat_type = MAT_TYPE_NONE;

    min_gamma = MIN_DEF;
    min_E = MIN_DEF;
    min_nu = MIN_DEF;
    min_rho = MIN_DEF;
    min_xi = MIN_DEF;
    min_Tvis = MIN_DEF;

    min_tss = MIN_DEF;
    min_tsl = MIN_DEF;
    min_Msat = MIN_DEF;
    min_chi = MIN_DEF;
    min_Tc = MIN_DEF;

    max_gamma = MAX_DEF;
    max_E = MAX_DEF;
    max_nu = MAX_DEF;
    max_rho = MAX_DEF;
    max_xi = MAX_DEF;
    max_Tvis = MAX_DEF;

    max_tss = MAX_DEF;
    max_tsl = MAX_DEF;
    max_Msat = MAX_DEF;
    max_chi = MAX_DEF;
    max_Tc = MAX_DEF;
}

void CPipeline::formatLine(string& line)
{
    string::size_type pos = 0;
    string tmp_str = seperateString(line);

    if (line.size() == 0)
        return;

    if (line.find(">") != string::npos)
    {
        pos = line.find(">");
        line.replace(pos, 1, "> ");
    }

    if (line.find("=") != string::npos)
    {
        pos = line.find("=");
        line.replace(pos, 1, " = ");
    }

    while (line.find(";") != string::npos)
    {
        pos = line.find(";");
        line.replace(pos, 1, " ");
    }

    while (line.find("?") != string::npos)
    {
        pos = line.find("?");
        line.replace(pos, 1, " ");
    }

    while (line.find("*") != string::npos)
    {
        pos = line.find("*");
        line.replace(pos, 1, " ");
    }

    while (line.find('\t') != string::npos)
    {
        pos = line.find('\t');
        line.replace(pos, 1, " ");
    }

    while (line.find(" \r\n") != string::npos)
    {
        pos = line.find(" \r\n");
        line.replace(pos, 3, " ");
    }

    while (line.find(" \r") != string::npos)
    {
        pos = line.find(" \r");
        line.replace(pos, 2, " ");
    }

    while (line.find(" \n") != string::npos)
    {
        pos = line.find(" \n");
        line.replace(pos, 2, " ");
    }

    while (line.find("\r\n") != string::npos)
    {
        pos = line.find("\r\n");
        line.replace(pos, 2, " ");
    }

    while (line.find("\r") != string::npos)
    {
        pos = line.find("\r");
        line.replace(pos, 1, " ");
    }

    while (line.find("\n") != string::npos)
    {
        pos = line.find("\n");
        line.replace(pos, 1, " ");
    }
    //***

    while (line.find("  ") != string::npos)
    {
        pos = line.find("  ");
        line.replace(pos, 2, " ");
    }

    while (line.find(",") != string::npos)
    {
        pos = line.find(",");
        line.replace(pos, 1, ".");
    }

    if (line == " ")
        line = "";

    if (line.size() > 0)
    {
        while (line.c_str()[line.size() - 1] == ' ')
        {
            pos = line.find_last_of(' ');
            line.erase(pos, 1);
        }
    }

    while (line.c_str()[0] == ' ')
    {
        // pos=line.find_first_of(' ');
        line.erase(0, 1);
    }

    if (line.find_first_of("#") != string::npos)
    {
        pos = line.find("#");
        if (pos == 0)
            tmp_str = "";

        line.erase(pos, line.length() - pos);
    }

    if (line.find_first_of("!") != string::npos)
    {
        pos = line.find("!");
        if (pos == 0)
            tmp_str = "";

        line.erase(pos, line.length() - pos);
    }

    if (tmp_str.size() != 0)
        line += " \"" + tmp_str + "\"";
}

bool CPipeline::writeOVITO(ullong time_id, double b_size, const vec3D* pos, const double* amon, const int* matID, int Nmon) 
{
    char str_tmp[2048];
    char str_end[2048];

    strcpy(str_tmp, "t_%05lu.dump");
    sprintf(str_end, str_tmp, time_id);

    string str_file = path_ovito + str_end;

    ofstream writer(str_file.c_str());

    if (writer.fail())
    {
        cout << "\nERROR: Cannot open file:\n" << str_file << endl;
        return false;
    }

    writer << "ITEM: TIMESTEP\n";
    writer << time_id << "\n";

    writer << "ITEM: NUMBER OF ATOMS\n";
    //writer << int(3 * Nmol)  << "\n";
    writer << Nmon << "\n";

    writer << "ITEM: BOX BOUNDS pp pp pp\n";

    strcpy(str_tmp, "%.4f %.4f\n");
    sprintf(str_end, str_tmp, -1.0e9 * b_size, 1.0e9 * b_size);

    writer << str_end;
    writer << str_end;
    writer << str_end;

    writer << "ITEM: ATOMS id mol type x y z radius\n"; //add velocity later

    for (int i = 0; i < Nmon; i++)
    {
        double x = 1.0e9 * pos[i].x;
        double y = 1.0e9 * pos[i].y;
        double z = 1.0e9 * pos[i].z;

        int id = matID[i]+1;

        double r = 1.0e9 * amon[i];

        strcpy(str_tmp, "%d %d %d %.5f %.5f %.5f %.5f\n");
        sprintf(str_end, str_tmp, i, 0, id, x, y, z, r);
        writer << str_end;
    }

    writer.close();

    return true;
}

bool CPipeline::writeAllOVITO(const vec3D* pos, const vec3D* vel, const vec3D* force, const vec3D* torque, const vec3D* omega, const vec3D* mag, const int* cluserIDs, const double* amon, const int* matID, int Nmon, int Nsto)
{
    char str_tmp[1024];
    char str_end[1024];

    double b_size = 0;

    for (int i = 0; i < Nsto; i++)
    {
        int index = i % Nmon;
        double a = amon[index];
        double distance = vec3D_length(pos[i]) + a;
        
        if (b_size < distance)
            b_size = distance;
    }

    b_size *= 1.15e9;

    if (b_size > 10000)
        b_size = 10000;

    int steps = Nsto / Nmon;

    double max_vel = 0;
    double max_force = 0;
    double max_torque = 0;
    double max_omega = 0;
    double max_mag = 0;

    for (int i = 0; i < Nmon; i++)
    {
        if (vel != 0)
        {
            double len_vel = vec3D_length(vel[i]);
            if (max_vel < len_vel)
                max_vel = sqrt(len_vel);
        }

        if (force != 0)
        {
            double len_force = vec3D_length(force[i]);
            if (max_force < len_force)
                max_force = sqrt(len_force);
        }

        if (torque != 0)
        {
            double len_torque = vec3D_length(torque[i]);
            if (max_torque < len_torque)
                max_torque = sqrt(len_torque);
        }

        if (omega != 0)
        {
            double len_omega = vec3D_length(omega[i]);
            if (max_omega < len_omega)
                max_omega = sqrt(len_omega);
        }

        if (mag != 0)
        {
            double len_mag = vec3D_length(mag[i]);
            if (max_mag < len_mag)
                max_mag = sqrt(len_mag);
        }
    }


    for (int i = 0; i < steps; i++)
    {
        strcpy(str_tmp, "t_%05lu.dump");
        sprintf(str_end, str_tmp, i);

        string str_file = path_ovito + str_end;

        ofstream writer(str_file.c_str());

        if (writer.fail())
        {
            cout << "\nERROR: Cannot open file:\n" << str_file << endl;
            return false;
        }

        writer << "ITEM: TIMESTEP\n";
        writer << i << "\n";

        writer << "ITEM: NUMBER OF ATOMS\n";
        //writer << int(3 * Nmol)  << "\n";
        writer << Nmon << "\n";

        writer << "ITEM: BOX BOUNDS pp pp pp\n";

        strcpy(str_tmp, "%.4f %.4f\n");
        sprintf(str_end, str_tmp, -b_size, b_size);

        writer << str_end;
        writer << str_end;
        writer << str_end;
        
        writer << "ITEM: ATOMS id mol type x y z vx vy vz fx fy fz tqx tqy tqz omegax omegay omegaz mux muy muz radius\n"; //add velocity later

        for (int j = 0; j < Nmon; j++)
        {
            double x = 1.0e9 * pos[i * Nmon + j].x;
            double y = 1.0e9 * pos[i * Nmon + j].y;
            double z = 1.0e9 * pos[i * Nmon + j].z;

            int cl_id=0;
            double vx = 0, vy = 0, vz = 0;
            double fx = 0, fy = 0, fz = 0;
            double tx = 0, ty = 0, tz = 0;

            double ox = 0, oy = 0, oz = 0;
            double mx = 0, my = 0, mz = 0;

            if(cluserIDs!=0)
                cl_id = cluserIDs[i * Nmon + j];

            if (vel != 0)
            {
                double len_vel = vec3D_length(vel[i * Nmon + j]);
                if (len_vel > 0)
                {
                    vx = vel[i * Nmon + j].x / len_vel;
                    vy = vel[i * Nmon + j].y / len_vel;
                    vz = vel[i * Nmon + j].z / len_vel;
                }
            }

            if (force != 0)
            {
                double len_force = vec3D_length(force[i * Nmon + j]);
                if (len_force > 0)
                {
                    fx = force[i * Nmon + j].x / len_force * pow(len_force,1.0/8.0);
                    fy = force[i * Nmon + j].y / len_force * pow(len_force, 1.0 / 8.0);
                    fz = force[i * Nmon + j].z / len_force * pow(len_force, 1.0 / 8.0);
                }
            }

            if (torque != 0)
            {
                double len_torque = vec3D_length(torque[i * Nmon + j]);
                if (len_torque > 0)
                {
                    tx = torque[i * Nmon + j].x / len_torque;
                    ty = torque[i * Nmon + j].y / len_torque;
                    tz = torque[i * Nmon + j].z / len_torque;
                }
            }

            if (omega != 0)
            {
                double len_omega = vec3D_length(omega[i * Nmon + j]);
                if (len_omega > 0)
                {
                    ox = omega[i * Nmon + j].x / len_omega;
                    oy = omega[i * Nmon + j].y / len_omega;
                    oz = omega[i * Nmon + j].z / len_omega;
                }
            }

            if (mag != 0)
            {
                double len_mag = vec3D_length(mag[i * Nmon + j]);
                if (len_mag != 0)
                {
                    //mx = mag[i * Nmon + j].x / max_mag;
                    //my = mag[i * Nmon + j].y / max_mag;
                    //mz = mag[i * Nmon + j].z / max_mag;

                    mx = mag[i * Nmon + j].x / len_mag * pow(len_mag, 1.0 / 8.0);
                    my = mag[i * Nmon + j].y / len_mag * pow(len_mag, 1.0 / 8.0);
                    mz = mag[i * Nmon + j].z / len_mag * pow(len_mag, 1.0 / 8.0);
                }
            }

            int mat_id = matID[j] + 1;

            double r = 1.0e9 * amon[j];

            strcpy(str_tmp, "%d %d %d %.5f %.5f %.5f %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5f\n");
            sprintf(str_end, str_tmp, j, cl_id, mat_id, x, y, z, vx, vy, vz, fx, fy, fz, tx, ty, tz, ox, oy, oz, mx, my, mz, r);
            writer << str_end;
        }

        writer.close();
    }

    return true;
}

string CPipeline::seperateString(string& str)
{
    string::size_type pos1 = 0, pos2 = 0;
    string ret = "";
    int len = -1;

    if (str.find_first_of("\"") != string::npos)
    {
        pos1 = str.find("\"");
        pos2 = str.find("\"", pos1 + 1);

        len = int(pos2 - pos1 + 1);

        if (len < 0)
            return ret;

        ret = str.substr(pos1, len);
        str.erase(pos1, len);
    }

    while (ret.find("\"") != string::npos)
    {
        pos1 = ret.find("\"");
        ret.erase(pos1, 1);
    }

    return ret;
}

dlist CPipeline::parseValues(string& str)
{
    int pos;
    dlist values;
    string v;

    formatLine(str);

    if (str.size() == 0)
        return values;

    while (str.find(" ") != string::npos)
    {
        pos = int(str.find(" "));
        v = str.substr(0, pos);
        str.erase(0, pos + 1);
        values.push_back(atof(v.c_str()));
    }

    values.push_back(atof(str.c_str()));

    return values;
}

bool CPipeline::isMatID(int id)
{
    int len = int(lst_matID.size());

    for (int i = 0; i < len; i++)
    {
        if (lst_matID[i] == i)
            return true;
    }

    return false;
}

bool CPipeline::prepareData(vec3D*& pos, vec3D*& vel, vec3D*& omega_tot, vec3D*& mag, double*& amon, double*& mass, double*& moment, int*& matIDs, int& Nmon)
{
    ifstream reader_A, reader_B;

    vlist lst_pos_A, lst_pos_B;
    ilist lst_matID_A, lst_matID_B;
    dlist lst_amon_A, lst_amon_B;

    string line_A, line_B;
    int line_counter_A = 0, line_counter_B = 0;

    agg_filename_A = path_A;
    agg_filename_B = path_B;

    reader_A.open(agg_filename_A.c_str());

    if (reader_A.fail())
    {
        cout << "\nERROR: Cannot read aggregate A from :\n\t" << agg_filename_A << "\n" << flush;
        return false;
    }

    reader_B.open(agg_filename_B.c_str());

    if (reader_B.fail())
    {
        cout << "\nERROR: Cannot read aggregate B from:\n\t" << agg_filename_B << "\n" << flush;
        return false;
    }

    a_mon_min = 1e200;
    a_mon_max = 0;

    while (getline(reader_A, line_A))
    {
        line_counter_A++;

        if (line_counter_A == 1)
        {
            dlist values = parseValues(line_A);

            Nmon_A = int(values[0]);
            a_eff_A = 1e-9 * values[2];
        }

        if (line_counter_A > 5)
        {
            dlist values = parseValues(line_A);

            if (values.size() == 0)
                continue;

            if (line_counter_A % 50 == 0)
                cout << "Reading aggregate A: " << 100.0 * float(line_counter_A) / float(Nmon_A) << "                 \r" << flush;

            vec3D tmp_pos;

            tmp_pos.x = 1e-9 * values[0];
            tmp_pos.y = 1e-9 * values[1];
            tmp_pos.z = 1e-9 * values[2];
            lst_pos_A.push_back(tmp_pos);

            double a_mon = 1e-9 * values[4];
            lst_amon_A.push_back(a_mon);

            if (a_mon_min > a_mon)
                a_mon_min = a_mon;

            if (a_mon_max < a_mon)
                a_mon_max = a_mon;

            double distance = sqrt(tmp_pos.x * tmp_pos.x + tmp_pos.y * tmp_pos.y + tmp_pos.z * tmp_pos.z) + a_mon;

            if (a_out_A < distance)
                a_out_A = distance;

            int mat_id = int(values[6]-1);

            if (!isMatID(mat_id))
            {
                cout << "ERROR: Material id " << mat_id << " in aggregate A line " << line_counter_A << " does not match the IDs in the command file!   \n" << flush;
                return false;
            }

            lst_matID_A.push_back(mat_id);
        }
    }

    while (getline(reader_B, line_B))
    {
        line_counter_B++;

        if (line_counter_B == 1)
        {
            dlist values = parseValues(line_B);

            Nmon_B = int(values[0]);
            a_eff_B = 1e-9 * values[2];
        }

        if (line_counter_B > 5)
        {
            dlist values = parseValues(line_B);

            if (values.size() == 0)
                continue;

            if (line_counter_B % 50 == 0)
                cout << "Reading aggregate B: " << 100.0 * float(line_counter_B) / float(Nmon_B) << "                 \r" << flush;

            vec3D tmp_pos;

            tmp_pos.x = 1e-9 * values[0];
            tmp_pos.y = 1e-9 * values[1];
            tmp_pos.z = 1e-9 * values[2];
            lst_pos_B.push_back(tmp_pos);

            double a_mon = 1e-9 * values[4];
            lst_amon_B.push_back(a_mon);

            if (a_mon_min > a_mon)
                a_mon_min = a_mon;

            if (a_mon_max < a_mon)
                a_mon_max = a_mon;

            double distance = sqrt(tmp_pos.x * tmp_pos.x + tmp_pos.y * tmp_pos.y + tmp_pos.z * tmp_pos.z) + a_mon;

            if (a_out_B < distance)
                a_out_B = distance;

            int mat_id = int(values[6]-1);

            if (!isMatID(mat_id))
            {
                cout << "ERROR: Material id " << mat_id << " in aggregate B line " << line_counter_B << " does not match the IDs in the command file!   \n" << flush;
                return false;
            }

            lst_matID_B.push_back(mat_id);
        }
    }

    cout << CLR_LINE;

    Nmon = Nmon_A + Nmon_B;

    pos = new vec3D[Nmon];
    vel = new vec3D[Nmon];
    mag = new vec3D[Nmon];
    omega_tot = new vec3D[Nmon];
    matIDs = new int[Nmon];
    amon = new double[Nmon];

    moment = new double[Nmon];
    mass = new double[Nmon];
    
    for (int i = 0; i < Nmon_A; i++)
    {
        pos[i].x = lst_pos_A[i].x + pos_A.x;
        pos[i].y = lst_pos_A[i].y + pos_A.y;
        pos[i].z = lst_pos_A[i].z + pos_A.z;

        amon[i] = lst_amon_A[i];

        vec3D r;

        r.x = 0;
        r.y = lst_pos_A[i].y;
        r.z = lst_pos_A[i].z; //todo: r needs to be perp. to ang_A

        vec3D vel_tan = vec3D_cross(ang_A, r);
         
        vel[i].x = vel_A.x + vel_tan.x;
        vel[i].y = vel_A.y + vel_tan.y;
        vel[i].z = vel_A.z + vel_tan.z;

        omega_tot[i] = ang_A;

        int mat_id = lst_matID_A[i];
        matIDs[i] = mat_id;

        double rho = lst_rho[mat_id];
        
        mass[i] = 4. / 3. * PI * rho * amon[i] * amon[i] * amon[i];
        moment[i] = 2. / 5. * mass[i] * amon[i] * amon[i];

        double chi = lst_chi[mat_id];
        double len_Bext = vec3D_length(Bext);

        if (abs(chi) != 0)
        {
            if (abs(chi) > LIMIT_FER) //ferromagnetic
            {
                if (len_Bext > 0.) //set initial direction to Bext
                {
                    mag[i].x = lst_Msat[mat_id] * Bext.x / len_Bext;
                    mag[i].y = lst_Msat[mat_id] * Bext.y / len_Bext;
                    mag[i].z = lst_Msat[mat_id] * Bext.z / len_Bext;

                    //mag[i].x = 0.984807753 * lst_Msat[mat_id];
                    //mag[i].y = 0;
                    //mag[i].z = 0.173648178 * lst_Msat[mat_id];
                }
                else
                {
                    double len_omega = vec3D_length(omega_tot[i]);
                    if (len_omega > 0.) //set initial direction to omega
                    {
                        mag[i].x = lst_Msat[mat_id] * omega_tot[i].x / len_omega;
                        mag[i].y = lst_Msat[mat_id] * omega_tot[i].y / len_omega;
                        mag[i].z = lst_Msat[mat_id] * omega_tot[i].z / len_omega;
                    }
                    else //set initial direction to x-direction
                    {
                        mag[i].x = lst_Msat[mat_id];
                        mag[i].y = 0.;
                        mag[i].z = 0.;
                    }
                }
            }
            else
            {
                vec3D M_ind, M_Bar;
                double chi_fac = chi / (chi + 1.);

                M_ind.x = chi_fac * Bext.x / mu0;
                M_ind.y = chi_fac * Bext.y / mu0;
                M_ind.z = chi_fac * Bext.z / mu0;

                M_Bar.x = chi * omega_tot[i].x / PROD_BARR;
                M_Bar.y = chi * omega_tot[i].y / PROD_BARR;
                M_Bar.z = chi * omega_tot[i].z / PROD_BARR;

                mag[i].x = M_ind.x + M_Bar.x;
                mag[i].y = M_ind.y + M_Bar.y;
                mag[i].z = M_ind.z + M_Bar.z;

                double len_mag = vec3D_length(mag[i]);

                if (len_mag> lst_Msat[mat_id])
                {
                    mag[i].x = lst_Msat[mat_id] * mag[i].x / len_mag;
                    mag[i].y = lst_Msat[mat_id] * mag[i].y / len_mag;
                    mag[i].z = lst_Msat[mat_id] * mag[i].z / len_mag;
                }

                //mag[i].x = 0.984807753 * lst_Msat[mat_id];
                //mag[i].y = 0;
                //mag[i].z = 0.173648178 * lst_Msat[mat_id];
            }
        }
        else
        {
            mag[i].x = 0;
            mag[i].y = 0;
            mag[i].z = 0;
        }
    }

    for (int i = Nmon_A; i < Nmon; i++)
    {
        int index = i - Nmon_A;
        pos[i].x = lst_pos_B[index].x + pos_B.x;
        pos[i].y = lst_pos_B[index].y + pos_B.y;
        pos[i].z = lst_pos_B[index].z + pos_B.z;

        amon[i] = lst_amon_B[index];

        vec3D r;

        r.x = 0;
        r.y = lst_pos_B[index].y;
        r.z = lst_pos_B[index].z;

        vec3D vel_tan = vec3D_cross(ang_B, r);

        vel[i].x = vel_B.x + vel_tan.x;
        vel[i].y = vel_B.y + vel_tan.y;
        vel[i].z = vel_B.z + vel_tan.z;

        omega_tot[i] = ang_B;

        int mat_id = lst_matID_B[index];
        matIDs[i] = mat_id;

        double rho = lst_rho[mat_id];

        mass[i] = 4. / 3. * PI * rho * amon[i] * amon[i] * amon[i];
        moment[i] = 2. / 5. * mass[i] * amon[i] * amon[i];

        double chi = lst_chi[mat_id];
        double len_Bext = vec3D_length(Bext);

        if (abs(chi) != 0)
        {
            if (abs(chi) > LIMIT_FER) //ferromagnetic
            {
                if (len_Bext > 0.) //set initial direction to Bext
                {
                    mag[i].x = lst_Msat[mat_id] * Bext.x / len_Bext;
                    mag[i].y = lst_Msat[mat_id] * Bext.y / len_Bext;
                    mag[i].z = lst_Msat[mat_id] * Bext.z / len_Bext;

                    //mag[i].x = 0.984807753 * lst_Msat[mat_id];
                    //mag[i].y = 0;
                    //mag[i].z = 0.173648178 * lst_Msat[mat_id];
                }
                else
                {
                    double len_omega = vec3D_length(omega_tot[i]);
                    if (len_omega > 0.) //set initial direction to omega
                    {
                        mag[i].x = lst_Msat[mat_id] * omega_tot[i].x / len_omega;
                        mag[i].y = lst_Msat[mat_id] * omega_tot[i].y / len_omega;
                        mag[i].z = lst_Msat[mat_id] * omega_tot[i].z / len_omega;
                    }
                    else //set initial direction to x-direction
                    {
                        mag[i].x = lst_Msat[mat_id];
                        mag[i].y = 0.;
                        mag[i].z = 0.;
                    }
                }
            }
            else
            {
                vec3D M_ind, M_Bar;
                double chi_fac = chi / (chi + 1.);

                M_ind.x = chi_fac * Bext.x / mu0;
                M_ind.y = chi_fac * Bext.y / mu0;
                M_ind.z = chi_fac * Bext.z / mu0;

                M_Bar.x = chi * omega_tot[i].x / PROD_BARR;
                M_Bar.y = chi * omega_tot[i].y / PROD_BARR;
                M_Bar.z = chi * omega_tot[i].z / PROD_BARR;

                mag[i].x = M_ind.x + M_Bar.x;
                mag[i].y = M_ind.y + M_Bar.y;
                mag[i].z = M_ind.z + M_Bar.z;

                double len_mag = vec3D_length(mag[i]);

                if (len_mag > lst_Msat[mat_id])
                {
                    mag[i].x = lst_Msat[mat_id] * mag[i].x / len_mag;
                    mag[i].y = lst_Msat[mat_id] * mag[i].y / len_mag;
                    mag[i].z = lst_Msat[mat_id] * mag[i].z / len_mag;
                }

                //mag[i].x = 0.984807753 * lst_Msat[mat_id];
                //mag[i].y = 0;
                //mag[i].z = 0.173648178 * lst_Msat[mat_id];
            }
        }
        else
        {
            mag[i].x = 0;
            mag[i].y = 0;
            mag[i].z = 0;
        }
    }

    //find smallest time step
    if (time_step == 0)
    {
        time_step = 1e200;

        for (int i = 0; i < Nmon; i++)
        {
            for (int j = 0; j < Nmon; j++)
            {
                if (i == j)
                    continue;

                int mat_id_A = matIDs[i];
                int mat_id_B = matIDs[j];

                double mass_A = mass[i];
                double mass_B = mass[j];
                double mass = (mass_A * mass_B) / (mass_A + mass_B);

                double a_mon_A = amon[i];
                double a_mon_B = amon[j];
                double R = (a_mon_A * a_mon_B) / (a_mon_A + a_mon_B);

                double nu_A = lst_nu[mat_id_A];
                double nu_B = lst_nu[mat_id_B];

                double E_A = lst_E[mat_id_A];
                double E_B = lst_E[mat_id_B];

                double Es = (1 - nu_A * nu_A) / E_A + (1 - nu_B * nu_B) / E_B;
                Es = 1. / Es;

                double gamma_A = lst_gamma[mat_id_A];
                double gamma_B = lst_gamma[mat_id_B];
                double gamma = gamma_A + gamma_B - 2. / (1. / gamma_A + 1. / gamma_B);

                double r0 = pow(9 * PI * gamma * R * R / Es, 1. / 3.);
                double delta_c = 0.5 * r0 * r0 / (R * pow(6., 1. / 3.));
                double F_c = 3.0 * PI * gamma * R;

                double tc = sqrt(mass * delta_c / F_c);

                //double tc = 0.95 * (pow(R, 7. / 6.) * sqrt(lst_rho[mat_id_A])) / (pow(gamma, 1. / 6.) * pow(Es, 1. / 3.));

                if (time_step > tc)
                    time_step = tc;
            }
        }

        int Nmat = int(lst_matName.size());

        for (int i = 0; i < Nmat; i++)
        {
            double tss = lst_tss[i];
            double tsl = lst_tsl[i];

            if (tss > 0)
            {
                if (time_step > tss)
                    time_step = tss;
            }

            if (tsl > 0)
            {
                if (time_step > tsl)
                    time_step = tsl;
            }
        }/**/

        time_step = 0.005 * time_step;
    }

    return true;
}

void CPipeline::printParameters()
{
    cout << SEP_LINE;
    cout << "\nSimulation parameters:\n";
    cout << "  - Nr. of iterations: " << N_iter << "\n";
    cout << "  - command file :\n\t" << cmd_filename << "\n";
    cout << "  - start time   : " << time_start << " sec\n";
    cout << "  - stop time    : " << time_stop << " sec\n";
    cout << "  - time step    : " << time_step << " sec\n\n";

    cout << "  - position A: (" << pos_A.x << ", " << pos_A.y << ", " << pos_A.z << ") m\n";
    cout << "  - position B: (" << pos_B.x << ", " << pos_B.y << ", " << pos_B.z << ") m\n\n";
    
    cout << "  - velocity A: (" << vel_A.x << ", " << vel_A.y << " " << vel_A.z << ") m/s\n";
    cout << "  - velocity B: (" << vel_B.x << ", " << vel_B.y << " " << vel_B.z << ") m/s\n\n";

    cout << "  - ang. velocity A: (" << ang_A.x << ", " << ang_A.y << ", " << ang_A.z << ") m\n";
    cout << "  - ang. velocity B: (" << ang_B.x << ", " << ang_B.y << ", " << ang_B.z << ") m\n\n";

    if (Bext.x + Bext.y + Bext.z == 0)
        cout << "  - ext. mag. field: none\n";
    else
        cout << "  - ext. mag. field: (" << Bext.x << ", " << Bext.y << ", " << Bext.z << ") T\n";

    cout << "  - dust temp.     : " << T_dust << " K \n\n";

    cout << SEP_LINE;
    cout << "\nMetarial parameters:\n";

    cout << "  - surface energy : " << min_gamma << " - " << max_gamma << " J m^-1     \n";
    cout << "  - Young's modulus: " << min_E << " - " << max_E << " Pa    \n";
    cout << "  - Poisson's ratio: " << min_nu << " - " << max_nu << "     \n";
    cout << "  - density        : " << min_rho << " - " << max_rho << " kg    \n";
    cout << "  - crit. roll.    : " << min_xi << " - " << max_xi << " m    \n";
    cout << "  - visc. damp time: " << min_Tvis << " - " << max_Tvis << " sec    \n";

    cout << "\n";
    if (mat_type == MAT_TYPE_MAG)
    {
        if (min_tss != MIN_DEF && max_tss != MAX_DEF)
            cout << "  - spin-spin rel. time   : " << min_tss << " - " << max_tss << " sec    \n";

        if (min_tsl != MIN_DEF && max_tsl != MAX_DEF)
            cout << "  - spin-lattice rel. time: " << min_tsl << " - " << max_tsl << " sec    \n";

        if (min_Msat != MIN_DEF && max_Msat != MAX_DEF)
            cout << "  - sat. magnetization    : " << min_Msat << " - " << max_Msat << " A m^-1    \n";

        if (min_chi != MIN_DEF && max_chi != MAX_DEF)
            cout << "  - mag. suceptibility    : " << min_chi << " - " << max_chi << "     \n";

        if (min_Tc != MIN_DEF && max_Tc != MAX_DEF)
            cout << "  - Curie temperature     : " << min_Tc << " - " << max_Tc << " K    \n";

        cout << "\n";
    }

    if (save_pos || save_vel || save_force || save_torque || save_cluster)
    {
        cout << SEP_LINE;
        cout << "\nOutput files every " << N_save << "-th time step: \n";

        if (save_pos)
            cout << "  - positions :\n\t" << path_binary << "pos.dat\n";

        if (save_vel)
            cout << "  - velocities:\n\t" << path_binary << "vel.dat\n";

        if (save_force)
            cout << "  - forces    :\n\t" << path_binary << "force.dat\n";

        if (save_torque)
            cout << "  - torques   :\n\t" << path_binary << "torque.dat\n";

        if (save_cluster)
            cout << "  - clusters  :\n\t" << path_binary << "cluster.dat\n";
    }

    cout << SEP_LINE;
    cout << "\nAggregate parameters:   \n";
    cout << "   - Nr. of monomers: " << Nmon_A+ Nmon_B << "\n";
    cout << "   - monomer radius : " << a_mon_min << " - " << a_mon_max << " m\n\n";

    cout << "   - file A:\n\t" << agg_filename_A << "\n";
    cout << "   - Nr. of monomers  A: " << Nmon_A << "\n";
    cout << "   - effective radius A: " << a_eff_A << " m\n";
    cout << "   - outer radius     A: " << a_out_A << " m\n\n";
        
    cout << "   - file B:\n\t" << agg_filename_B << "\n";
    cout << "   - Nr. of monomers  B: " << Nmon_B << "\n";
    cout << "   - effective radius B: " << a_eff_B << " m\n";
    cout << "   - outer radius     B: " << a_out_B << " m\n";
    cout << SEP_LINE << flush;
}

double CPipeline::getMin(dlist& lst)
{
    double min = *std::max_element(lst.begin(), lst.end());
    return min;
}

double CPipeline::getMax(dlist & lst)
{
    double max = *std::max_element(lst.begin(), lst.end());
    return max;
}

bool CPipeline::init(int argc, const char** argv)
{
    cout << SEP_LINE;
    cout << PROG_ID;
    cout << SEP_LINE << flush;

#ifdef DEBUG
    cmd_filename = "/home/ilion/0/pzuern/development/TestFiles/squeeze/cmd_file";
#elif RELEASE
    if (argc != 2)
    {
        cout << "\nERROR: Wrong number of arguments!                     \n";
        cout << "       DUST COLLIDER requires only the path of a command file!            \n";
        cout << SEP_LINE;
        return false;
    }
    cmd_filename = argv[1];
#endif

    return true;
}


bool CPipeline::createPath(string path)
{
    /*cout << path << endl;

    path = "F:\\test\\";
    cout << path << endl;*/

#ifdef _WIN32
    // Windows-specific 
    if (_mkdir(path.c_str()) == 0) 
    {
        cout << "Directory created successfully on Windows:\n\t" << path << std::endl;
        return true;
    }
    else 
    {
        if (errno == EEXIST)
        {
            cout << "Directory exists:\n\t" << path << std::endl;
            return true;
        }
        else
        {
            if (errno == ENOENT)
            {
                cout << "ERROR: Path does not exist:\n\t" << path << std::endl;
                return false;

                /*stringstream str_strem;
                string segment;
                vector<std::string> seglist;

                str_strem << path;
                
                while (getline(str_strem, segment, SEP))
                {
                    seglist.push_back(segment);
                }

                //for(int i=0;i< seglist.size(), i++)*/


            }
            else
            {
                cout << "Failed to create directory on Windows! " << errno << std::endl;
                cout << "Directory:\n\t" << path << std::endl;

                return false;
            }
        }
    }
#elif __linux__
    // Linux-specific 
    if (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) == 0) 
    {
        cout << "Directory created successfully on Linux:\n\t" << path << std::endl;
        return true;
    }
    else 
    {
        if (errno == EEXIST)
        {
            cout << "Directory exists:\n\t" << path << std::endl;
            return true;
        }
        else
            cout << "Failed to create directory on Linux: " << strerror(errno) << std::endl;
    }
#else
    cout << "Unsupported OS" << std::endl;
#endif
    return false;
}

bool CPipeline::parse()
{
    ifstream reader(cmd_filename.c_str());

    string line, cmd;
    string::size_type pos = 0;

    int line_counter = 0;

    if (reader.fail())
    {
        cout << "\nERROR: Cannot open command file:\n\t" << cmd_filename << endl;
        return false;
    }

    while (getline(reader, line))
    {
        line_counter++;

        formatLine(line);

        if (line.compare("") == 0)
            continue;

        if (line.c_str()[0] == '#')
            continue;


        pos = line.find(">");

        if (pos != string::npos)
        {
            cmd = line.substr(0, pos + 1);
            line.erase(0, pos + 2);
        }
        else
        {
            cout << "ERROR: Cannot identify tag \"" << line << "\" in line " << line_counter << "!    \n" << flush;
            return false;
        }

        if (cmd.compare("") == 0)
            continue;

        if (cmd.compare(" ") == 0)
            continue;

        if (!parseLine(cmd, line))
        {
            cout << "ERROR: Cannot parse command \"" << cmd << "\" in line " << line_counter << endl << flush;
            return false;
        }
    }
        
    return true;
}

bool CPipeline::writeHeader()
{
    string str_file = path_binary;
    str_file += "header.txt";

    int Nmat = int(lst_gamma.size());

    ofstream writer(str_file.c_str());

    if (writer.fail())
    {
        cout << "\nERROR: Cannot write header file:\n\t" << str_file << endl;
        return false;
    }

    writer << "#N_iter\tNmon_A\tNmon_B\tN_save\tN_mat\ttime_step [s]\n";
    writer << N_iter << "\t" << Nmon_A << "\t" << Nmon_B << "\t" << N_save << "\t" << Nmat << "\t" << time_step << "\n";

    //further header information can be written here...

    writer.close();
    return true;
}

bool CPipeline::writeBinaryVec(string name_file, const vec3D* data, ullong N)
{
    string path_tmp = path_binary + name_file;
    ullong len_array = N * sizeof(vec3D);

    ofstream bin_writer(path_tmp.c_str(), ios::binary);

    if (bin_writer.fail())
    {
        cout << "\nERROR: Cannot write binary file:\n\t";
        cout << path_tmp << "\n";
        return false;
    }

    cout << "   - Writing binary file:\n\t";
    cout << path_tmp << "\n";

    bin_writer.write((const char*)data, len_array);

    return true;
}

bool CPipeline::writeBinaryDouble(string name_file, const double* data, ullong N)
{
    string path_tmp = path_binary + name_file;
    ullong len_array = N * sizeof(double);

    ofstream bin_writer(path_tmp.c_str(), ios::out | ios::binary);

    if (bin_writer.fail())
    {
        cout << "\nERROR: Cannot write binary file:\n\t";
        cout << path_tmp << "\n";
        return false;
    }

    cout << "   - Writing binary file:\n\t";
    cout << path_tmp << "\n";

    bin_writer.write((const char*)data, len_array);

    return true;
}


bool CPipeline::writeBinaryInt(string name_file, const int* data, ullong N)
{
    string path_tmp = path_binary + name_file;
    ullong len_array = N * sizeof(int);

    ofstream bin_writer(path_tmp.c_str(), ios::out | ios::binary);

    if (bin_writer.fail())
    {
        cout << "\nERROR: Cannot write binary file:\n\t";
        cout << path_tmp << "\n";
        return false;
    }

    cout << "   - Writing binary file:\n\t";
    cout << path_tmp << "\n";

    bin_writer.write((const char*)data, len_array);

    return true;
}

bool CPipeline::parseLine(string cmd, string data)
{
    if (cmd.compare("<path_results>") == 0)
    {
        path_results = seperateString(data);

        path_binary = path_results + "binary" + SEP;
        path_ovito = path_results + "ovito" + SEP;
        path_plots = path_results + "plots" + SEP;

        return true;
    }

    if (cmd.compare("<path_A>") == 0)
    {
        path_A = seperateString(data);
        return true;
    }

    if (cmd.compare("<path_B>") == 0)
    {
        path_B = seperateString(data);
        return true;
    }

    if (cmd.compare("<pos_A>") == 0)
    {
        dlist values = parseValues(data);

        if (values.size() != 3)
        {
            cout << "ERROR: Wrong amount of coordinates defined in command " << cmd << "!      \n";
            return false;
        }

        pos_A.x = values[0];
        pos_A.y = values[1];
        pos_A.z = values[2];

        return true;
    }

    if (cmd.compare("<pos_B>") == 0)
    {
        dlist values = parseValues(data);

        if (values.size() != 3)
        {
            cout << "ERROR: Wrong amount of coordinates defined in command " << cmd << "!      \n";
            return false;
        }

        pos_B.x = values[0];
        pos_B.y = values[1];
        pos_B.z = values[2];

        return true;
    }

    if (cmd.compare("<B_ext>") == 0)
    {
        dlist values = parseValues(data);

        if (values.size() != 3)
        {
            cout << "ERROR: Wrong amount of coordinates defined in command " << cmd << "!      \n";
            return false;
        }

        Bext.x = values[0];
        Bext.y = values[1];
        Bext.z = values[2];

        return true;
    }


    if (cmd.compare("<vel_A>") == 0)
    {
        dlist values = parseValues(data);

        if (values.size() != 3)
        {
            cout << "ERROR: Wrong amount of coordinates defined in command " << cmd << "!      \n";
            return false;
        }

        vel_A.x = values[0];
        vel_A.y = values[1];
        vel_A.z = values[2];

        return true;
    }

    if (cmd.compare("<vel_B>") == 0)
    {
        dlist values = parseValues(data);

        if (values.size() != 3)
        {
            cout << "ERROR: Wrong amount of coordinates defined in command " << cmd << "!      \n";
            return false;
        }

        vel_B.x = values[0];
        vel_B.y = values[1];
        vel_B.z = values[2];

        return true;
    }

    if (cmd.compare("<ang_A>") == 0)
    {
        dlist values = parseValues(data);

        if (values.size() != 3)
        {
            cout << "ERROR: Wrong amount of coordinates defined in command " << cmd << "!      \n";
            return false;
        }

        ang_A.x = values[0];
        ang_A.y = values[1];
        ang_A.z = values[2];

        return true;
    }

    if (cmd.compare("<ang_B>") == 0)
    {
        dlist values = parseValues(data);

        if (values.size() != 3)
        {
            cout << "ERROR: Wrong amount of coordinates defined in command " << cmd << "!      \n";
            return false;
        }

        ang_B.x = values[0];
        ang_B.y = values[1];
        ang_B.z = values[2];

        return true;
    }

    if (cmd.compare("<T_dust>") == 0)
    {
        double t = double(atof(data.c_str()));
        T_dust = t;

        return true;
    }

    if (cmd.compare("<time_start>") == 0)
    {
        double t = double(atof(data.c_str()));
        time_start = t;

        return true;
    }

    if (cmd.compare("<time_stop>") == 0)
    {
        double t = double(atof(data.c_str()));
        time_stop = t;

        return true;
    }

    if (cmd.compare("<time_step>") == 0)
    {
        double t = double(atof(data.c_str()));
        time_stop = t;

        return true;
    }

    if (cmd.compare("<N_iter>") == 0)
    {
        ullong i = ullong(atof(data.c_str()));
        N_iter = i;

        return true;
    }



    if (cmd.compare("<N_save>") == 0)
    {
        ullong i = ullong(atof(data.c_str()));
        N_save = i;

        return true;
    }

    if (cmd.compare("<save_pos>") == 0)
    {
        int s = int(atof(data.c_str()));
        
        if (s == 1)
            save_pos = true;

        return true;
    }

    if (cmd.compare("<save_vel>") == 0)
    {
        int s = int(atof(data.c_str()));

        if (s == 1)
            save_vel = true;

        return true;
    }

    if (cmd.compare("<save_force>") == 0)
    {
        int s = int(atof(data.c_str()));

        if (s == 1)
            save_force = true;

        return true;
    }

    if (cmd.compare("<save_torque>") == 0)
    {
        int s = int(atof(data.c_str()));

        if (s == 1)
            save_torque = true;

        return true;
    }

    if (cmd.compare("<save_cluster>") == 0)
    {
        int s = int(atof(data.c_str()));

        if (s == 1)
            save_cluster = true;

        return true;
    }

    if (cmd.compare("<save_omega>") == 0)
    {
        int s = int(atof(data.c_str()));

        if (s == 1)
            save_omega = true;

        return true;
    }

    if (cmd.compare("<save_mag>") == 0)
    {
        int s = int(atof(data.c_str()));

        if (s == 1)
            save_mag = true;

        return true;
    }

    if (cmd.compare("<save_ovito>") == 0)
    {
        int s = int(atof(data.c_str()));

        if (s == 1)
            save_ovito = true;

        return true;
    }

    if (cmd.compare("<material id = >") == 0)
    {
        string str_m = seperateString(data);
        string str_id = seperateString(data);
        int ID = int(atof(str_id.c_str()));
        //formatLine()

        dlist values = parseValues(data);

        int length = int(values.size());

        if (values.size() == 6)
        {
            lst_matName.push_back(str_m);
            lst_matID.push_back(ID - 1);

            lst_gamma.push_back(values[0]);
            lst_E.push_back(values[1]);
            lst_nu.push_back(values[2]);
            lst_rho.push_back(values[3]);
            lst_xi.push_back(values[4]);
            lst_Tvis.push_back(values[5]);

            lst_tss.push_back(0);
            lst_tsl.push_back(0);
            lst_Msat.push_back(0);
            lst_chi.push_back(0);
            lst_Tc.push_back(0);

            return true;
        }

        if (values.size() == 11)
        {
            mat_type = MAT_TYPE_MAG;

            lst_matName.push_back(str_m);
            lst_matID.push_back(ID - 1);

            lst_gamma.push_back(values[0]);
            lst_E.push_back(values[1]);
            lst_nu.push_back(values[2]);
            lst_rho.push_back(values[3]);
            lst_xi.push_back(values[4]);
            lst_Tvis.push_back(values[5]);

            lst_tss.push_back(values[6]);
            lst_tsl.push_back(values[7]);
            lst_Msat.push_back(values[8]);
            lst_chi.push_back(values[9]);
            lst_Tc.push_back(values[10]);

            return true;
        }        

        cout << "ERROR: Wrong amount of material constants in command file!      \n";
        return false;
    }

    return false;
}

void CPipeline::prepareMaterial(material*& mat, int& Nmat)
{
    Nmat = int(lst_matName.size());

    mat = new material[Nmat];
    for (int i = 0; i < Nmat; i++)
    {
        mat[i].gamma = lst_gamma[i];
        mat[i].E = lst_E[i];
        mat[i].nu = lst_nu[i];
        mat[i].rho = lst_rho[i];
        mat[i].xi = lst_xi[i];
        mat[i].tvis = lst_Tvis[i];

        mat[i].tss = lst_tss[i];
        mat[i].tsl = lst_tsl[i];
        mat[i].Msat = lst_Msat[i];
        mat[i].chi = lst_chi[i];
        mat[i].Tc = lst_Tc[i];
        
        //T_dust == -1: Impact of dust temp. is not considdered
        if (T_dust != -1)
        {
            double corr;

            // Curie's law
            corr = CHI_20 / T_dust;
            //mat[i].chi *= corr;

            if (abs(mat[i].Msat) > 0 && mat[i].Tc > 0)
            {
                //Bloch T^3/2 law
                corr = 1.0 - pow(T_dust / mat[i].Tc, 1.5);

                if (corr < 0.0)
                    corr = 0.0;

                //mat[i].Msat *= corr;
            }

            // reduced surface energy Bogdan+ 2020
            corr = SLOPE * T_dust + INTERCEPT;
            //mat[i].gamma *= corr;
        }

        if (min_gamma > mat[i].gamma)
            min_gamma = mat[i].gamma;

        if (max_gamma < mat[i].gamma)
            max_gamma = mat[i].gamma;


        if (min_E > mat[i].E)
            min_E = mat[i].E;

        if (max_E < mat[i].E)
            max_E = mat[i].E;


        if (min_nu > mat[i].nu)
            min_nu = mat[i].nu;

        if (max_nu < mat[i].nu)
            max_nu = mat[i].nu;


        if (min_rho > mat[i].rho)
            min_rho = mat[i].rho;

        if (max_rho < mat[i].rho)
            max_rho = mat[i].rho;


        if (min_xi > mat[i].xi)
            min_xi = mat[i].xi;

        if (max_xi < mat[i].xi)
            max_xi = mat[i].xi;


        if (min_Tvis > mat[i].tvis)
            min_Tvis = mat[i].tvis;

        if (max_Tvis < mat[i].tvis)
            max_Tvis = mat[i].tvis;

        if (mat[i].tss != 0)
        {
            if (min_tss > mat[i].tss)
                min_tss = mat[i].tss;

            if (max_tss < mat[i].tss)
                max_tss = mat[i].tss;
        }

        if (mat[i].tsl != 0)
        {
            if (min_tsl > mat[i].tsl)
                min_tsl = mat[i].tsl;

            if (max_tsl < mat[i].tsl)
                max_tsl = mat[i].tsl;
        }

        if (mat[i].Msat != 0)
        {
            if (min_Msat > mat[i].Msat)
                min_Msat = mat[i].Msat;

            if (max_Msat < mat[i].Msat)
                max_Msat = mat[i].Msat;
        }

        if (mat[i].chi != 0)
        {
            if (min_chi > mat[i].chi)
                min_chi = mat[i].chi;

            if (max_chi < mat[i].chi)
                max_chi = mat[i].chi;
        }

        if (mat[i].Tc != 0)
        {
            if (min_Tc > mat[i].Tc)
                min_Tc = mat[i].Tc;

            if (max_Tc < mat[i].Tc)
                max_Tc = mat[i].Tc;
        }
    }
}

bool CPipeline::checkParameters()
{
    int len = int(lst_matID.size());

    if (len == 0)
    {
        cout << "ERROR: No materials defined!  \n";
        return false;
    }

    if (N_save <= 0)
    {
        cout << "ERROR: Command \"<N_save>\" needs to be larger than zero!  \n";
        return false;
    }

    if (lst_matID[0] != 0)
    {
        cout << "ERROR: Lowest material ID has to be 1!  \n";
        return false;
    }

    for (int i = 0; i < len; i++)
    {
        if (lst_matID[i]!=i)
        {
            cout << "ERROR: Material IDs not incremented by 1!  \n";
            return false;
        }
    }

    for (int i = 0; i < len; i++)
    {
        if ( (abs(lst_Msat[i]) == 0 && lst_Tc[i] > 0) || (abs(lst_Msat[i]) > 0 && lst_Tc[i] == 0))
        {
            cout << "ERROR: Spont. magnetization and Curie temperature connot both be zero!  \n";
            return false;
        }
    }

    if (mat_type == MAT_TYPE_NONE)
    {
        if (vec3D_length(Bext) > 0)
        {
            cout << "WARNING: None of the material can be magnetized but an external B-field is defined! \n";
        }
    }

    if (mat_type == MAT_TYPE_MAG)
    {
        if (vec3D_length(Bext) == 0)
        {
            cout << "WARNING: Some material can be magnetized but no external B-field is defined! \n";
        }
    }

    if (path_results.length() > 5)
    {
        if (!createPath(path_results))
            return false;

        if (!createPath(path_binary))
            return false;

        if (!createPath(path_ovito))
            return false;

        if (!createPath(path_plots))
            return false;
    }
    else
    {
        cout << "ERROR: Invalid path for simulation results!  \n";
        return false;
    }

    if (path_A.length() <= 5)
    {
        cout << "ERROR: Invalid path for aggregate A!  \n";
        return false;
    }

    if (path_B.length() <= 5)
    {
        cout << "ERROR: Invalid path for aggregate B!  \n";
        return false;
    }
     
    //TODO: check overlap between entire aggregates

    double len_pos_A = sqrt(pos_A.x * pos_A.x + pos_A.y * pos_A.y + pos_A.z * pos_A.z);
    double len_pos_B = sqrt(pos_B.x * pos_B.x + pos_B.y * pos_B.y + pos_B.z * pos_B.z);

    if(len_pos_A + len_pos_B == 0)
    {
        cout << "ERROR: Both aggregate postions are at the center! \n";
        return false;
    }

    if (T_dust <= 0 && T_dust != -1)
    {
        cout << "ERROR: Dust temperature needs to be larger than zero! \n";
        return false;
    }

    if (T_dust == -1)
    {
        cout << "WARNING: Impact of dust temperature is ignored! \n";
    }

    double len_vel_A = sqrt(vel_A.x * vel_A.x + vel_A.y * vel_A.y + vel_A.z * vel_A.z);
    double len_vel_B = sqrt(vel_B.x * vel_B.x + vel_B.y * vel_B.y + vel_B.z * vel_B.z);

    if (len_vel_A + len_vel_B == 0)
    {
        cout << "WARNING: Both aggregate velocities are zero! \n";
        //return false;
    }

    double len_ang_A = sqrt(ang_A.x * ang_A.x + ang_A.y * ang_A.y + ang_A.z * ang_A.z);
    double len_ang_B = sqrt(ang_B.x * ang_B.x + ang_B.y * ang_B.y + ang_B.z * ang_B.z);

    if (len_ang_A + len_ang_B == 0)
    {
        cout << "WARNING: Non rotating aggregates! \n";
    }

    if (N_iter <= 0)
    {
        cout << "ERROR: Invalid number of iterations! \n";
        return false;
    }

    if (time_start < 0)
    {
        cout << "WARNING: Invalid starting time! \n";
        //return false;
    }

    if (time_stop <= 0)
    {
        cout << "WARNING: Invalid stopping time! \n";
        //return false;
    }

    //TODO: check time step?
    return true;
}