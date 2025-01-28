#pragma once

/**
 * @brief A struct of pointers to host memory of arrays containing static monomer properties.
 */
typedef struct {
    double* mass;
    double* magnetic_susceptibiliy;
} hostMaterialProperties;

/**
 * @brief A struct of pointers to device memory of arrays containing static monomer properties.
 */
typedef struct {
    double* mass;
    double* magnetic_susceptibiliy;
} deviceMaterialProperties;