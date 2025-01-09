#include <vector>
using namespace std;

#ifndef TYPEDEFS
#define TYPEDEFS

/**
 * @brief A makro to determine if the code compiles the GPU versions of most function.
 */
#define RUN_ON_GPU

/**
 * @brief UNKNOWN PURPOSE
 */
#define BLOCK_SIZE 256

// TODO: Move this to where it is actually used...
/**
 * @brief A starting value for a maximum number algorithm.
 */
#define MIN_DEF  1e100

// TODO: Move this to where it is actually used...
/**
 * @brief A starting value for a minimum number algorithm.
 */
#define MAX_DEF -1e100

/**
 * A makro for pi.
 */
#define PI 3.1415926535897932

/**
 * A makro for four times pi.
 */
#define PIx4 (4*3.1415926535897932)

/**
 * @brief A makro for the vacuum permeability [N * A**2]
 */
#define mu0       1.25663706127e-6

/**
 * @brief A makro for the gyromagnetic ratio [rad s**-1 T**-1]
 */
#define GAMMA_E   1.76085963052e11

/** 
 * @brief A makro for the ice point temperature [K] 
 */
#define CHI_20    293.15

/**
 * @brief A makro for the product of vacuum permeability and gyromagnetic ratio used in Barnett effect calculation.
 */
#define PROD_BARR (mu0 * GAMMA_E)

//reduced surface energy Bogdan+ 2020
/**
 * @brief A makro for the intercept value for the reduced surface energy (Bogdan et al 2020)
 */
#define INTERCEPT 1.0

/**
 * @brief A makro for the slope value for the reduced surface energy (Bogdan et al 2020)
 */
#define SLOPE	 -7.018E-04

/**
 * @brief A parameter that determines the softening strenght as per (citation here).
 */
#define SOFTENING  1.0e-9f

/**
 * @brief A makro that determins which materials will be classified as ferromagnetic.
 */
#define LIMIT_FER 0.1

/**
 * @brief A makro that contains the ID of nonmagnetic materials (???).
 */
#define MAT_TYPE_NONE 0

/**
 * @brief A makro that contains the ID of magnetic materials (???).
 */
#define MAT_TYPE_MAG  1

/**
 * @brief A struct representing a 3D vector.
 * 
 * A struct with fields `x`, `y` and `z` representing the components of a 3D vector.
 */
typedef struct
{
    double x, y, z;
} vec3D;

/**
 * @brief A struct representing a quaternion.
 * 
 * A struct with the fields `e0`, `e1`, `e2` and `e3` representing the components of a quaternion.
 */
typedef struct
{
    double e0, e1, e2, e3;
} quat;

typedef unsigned long long ullong;

/**
 * A standard container for double values.
 */
typedef vector<double> dlist;

/**
 * A standard container for integer values.
 */
typedef vector<int> ilist;

/**
 * A standard container for 3D vectors.
 */
typedef vector<vec3D> vlist;

/**
 * A standard container for strings.
 */
typedef vector<string> strlist;

/**
 * @brief A struct containing the all material parameters for a single material.
 * 
 * @param gamma: surface energy.
 * @param E: Young's modulus
 * @param nu: Poissons ratio
 * @param rho: density
 * @param xi: critical rolling displacement
 * @param tvis: viscous damping time scale
 * 
 * @param tss: spin-spin relaxation time
 * @param tsl: spin-lattice relaxation time
 * @param Msat: saturation magnetization
 * @param chi: magnetic susceptibility
 * @param Tc: Curie temperature
 */
typedef struct
{
    double gamma;
    double E;
    double nu;
    double rho;
    double xi;
    double tvis;

    double tss;
    double tsl;
    double Msat;
    double chi;
    double Tc;
} material;

#ifdef VERSION
    #define PROG_ID "                         DUST COLLIDER    " VERSION "                              \n"
#else
    #define PROG_ID "                         DUST COLLIDER    undefined version                        \n"
#endif

/**
 * A makro for a seperator line for console printing.
 */
#define SEP_LINE "*************************************************************************************\n"

/**
 * A makro for a clear line for console printing
 */
#define CLR_LINE "                                                                                     \r"

#ifdef _WIN32
    /**
     * A makro for the path seperator \\.
     */
    #define SEP '\\'
#elif __linux__
    /** 
     * A makro for the path seperator /.
     */
    #define SEP '/'
#endif

#endif
