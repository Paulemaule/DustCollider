#pragma once

///////////////////////// PHYISCAL CONSTANTS /////////////////////////

// A makro for pi.
#define PI 3.1415926535897932

// A makro for 4 * pi
#define PIx4 (4*3.1415926535897932)

// A makro for the vacuum permeability [N * A**2]
#define mu0       1.25663706127e-6

// A makro for the gyromagnetic ratio [rad s**-1 T**-1]
#define GAMMA_E   1.76085963052e11

// A makro for the ice point temperature [K] 
#define CHI_20    293.15

// A makro for the product of vacuum permeability and gyromagnetic ratio used in Barnett effect calculation.
#define PROD_BARR (mu0 * GAMMA_E)

///////////////////////// MODEL PARAMETERS /////////////////////////

// A makro for the intercept value for the reduced surface energy (Bogdan et al 2020)
#define INTERCEPT 1.0

// A makro for the slope value for the reduced surface energy (Bogdan et al 2020)
#define SLOPE	 -7.018E-04

// A parameter that determines the softening strenght (citation here).
#define SOFTENING  1.0e-9f

///////////////////////// MAGNETIC MODELLING /////////////////////////

// TODO: Rework this.
// A makro that determins which materials will be classified as ferromagnetic.
#define LIMIT_FER 0.1

// TODO: Rework this.
// A makro used to track if there are magnetic materials and if magnetic modelling should be active.
#define MAT_TYPE_NONE 0

// TODO: Rework this.
// A makro used to track if there are magnetic materials and if magnetic modelling should be active.
#define MAT_TYPE_MAG  1