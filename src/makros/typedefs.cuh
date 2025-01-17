#include <vector>
using namespace std;

#pragma once

// FIXME: Replace with double3
/**
 * @brief A struct representing a 3D vector.
 * 
 * A struct with fields `x`, `y` and `z` representing the components of a 3D vector.
 */
typedef struct
{
    double x, y, z;
} vec3D;

// FIXME: Replace with double4
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