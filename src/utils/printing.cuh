#pragma once

#include <ctime>
#include <iostream>
#include <sstream>
#include <string>

#include "config.cuh"   // VERBOSITY, referenced by the PRINT_LOG macro below

///////////////////////// PRINT FORMAT /////////////////////////

#define SEP_LINE    "*************************************************************************************"
#define CLR_LINE    "                                                                                     "

///////////////////////// PRINTING /////////////////////////

/**
 * @brief Converts a duration in nanoseconds into a string of the form HHHH:MM:SS.mm.
 * 
 * @param duration: The duration that is to be converted into a string.
 * @param buffer: A buffer the resulting string is to be written into, needs to have length .
 */
void ns_to_time_string(const long duration, char* buffer, size_t buffer_size) {
    if (buffer_size < 14) {
        Logger::error("The supplied buffer was too small. Buffer needs to be at least 14 chars long.");
    }

    long remaining = duration;

    long hours = remaining / 3'600'000'000'000;
    remaining %= 3'600'000'000'000;

    long minutes = remaining / 60'000'000'000;
    remaining %= 60'000'000'000;

    double seconds = static_cast<double>(remaining) / 1'000'000'000.0;
    
    std::snprintf(buffer, buffer_size, "%04ld:%02ld:%05.02f", hours, minutes, seconds);
}

/**
 * @brief Returns the current wall-clock time as a "YYYY-MM-DD HH:MM:SS" string.
 */
inline std::string current_time_string() {
    std::time_t now = std::time(nullptr);
    char buffer[20];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return buffer;
}

///////////////////////// DEBUGGING /////////////////////////
void print_double (double* array, int start, int stop) {
    for (int i = start; i < stop; i++) {
        printf("%4d: %12.3e | ", i, array[i]);
    }
    printf("\n");
}

void print_double3 (double3* array, int start, int stop) {
    for (int i = start; i < stop; i++) {
        printf("%4d: %12.3e | ", i, array[i].x);
    }
    printf("\n");
    for (int i = start; i < stop; i++) {
        printf("%4d: %12.3e | ", i, array[i].y);
    }
    printf("\n");
    for (int i = start; i < stop; i++) {
        printf("%4d: %12.3e | ", i, array[i].z);
    }
    printf("\n");
}

void print_double4 (double4* array, int start, int stop) {
    for (int i = start; i < stop; i++) {
        printf("%4d: %12.3e | ", i, array[i].w);
    }
    printf("\n");
    for (int i = start; i < stop; i++) {
        printf("%4d: %12.3e | ", i, array[i].x);
    }
    printf("\n");
    for (int i = start; i < stop; i++) {
        printf("%4d: %12.3e | ", i, array[i].y);
    }
    printf("\n");
    for (int i = start; i < stop; i++) {
        printf("%4d: %12.3e | ", i, array[i].z);
    }
    printf("\n");
}