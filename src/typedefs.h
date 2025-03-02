#include <vector>
using namespace std;

#ifndef TYPEDEFS
#define TYPEDEFS

#define BLOCK_SIZE 256
#define SOFTENING  1.0e-9f

#define MIN_DEF  1e100
#define MAX_DEF -1e100

#define mu0       1.25663706127e-6  // vacuum permeability N A^-2
#define GAMMA_E   1.76085963052e11  // gyromagnetic ratio rad s^-1 T^-1
#define CHI_20    293.15            // rescaling of susceptibility in K
#define PROD_BARR (mu0 * GAMMA_E)   // product to calculate Barnett effect

//reduced surface energy Bogdan+ 2020
#define INTERCEPT 1.0
#define SLOPE	 -7.018E-04

#define LIMIT_FER 0.1  //ferromagnetic limit

//material type
#define MAT_TYPE_NONE 0
#define MAT_TYPE_MAG  1

typedef struct
{
	double x, y, z;
} vec3D;

typedef struct
{
	double e0, e1, e2, e3;
} quat;

typedef struct
{
    double gamma; // surface energy
    double E;     // Young's modulus
    double nu;    // Poisson ratio
    double rho;   // density
    double xi;    // critical rolling displacement
    double tvis;  // viscous damping time scale

    double tss;   // spin-spin relaxation time
    double tsl;   // spin-lattice relaxation time
    double Msat;  // saturation magnetization
    double chi;   // magnetic susceptibility
    double Tc;    // Curie temperature
} material;

typedef unsigned long long ullong;

typedef vector<double> dlist; 
typedef vector<int> ilist;
typedef vector<vec3D> vlist;
typedef vector<string> strlist;

#define PI 3.1415926535897932
#define PIx4 (4*3.1415926535897932)

#define PROG_ID  "                              DUST COLLIDER    V0.00.00                              \n"
#define SEP_LINE "*************************************************************************************\n"
#define CLR_LINE "                                                                                     \r"



#ifdef _WIN32
#define SEP '\\'
#elif __linux__
#define SEP '/'
#endif



#endif
