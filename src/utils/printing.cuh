#pragma once

#include <sstream>

///////////////////////// PRINT FORMAT /////////////////////////

#define SEP_LINE    "*************************************************************************************"
#define CLR_LINE    "                                                                                     "

///////////////////////// PRINTING /////////////////////////

#define PRINT_SEP_LINE() {std::cout << SEP_LINE << std::endl << std::flush;}
#define PRINT_CLR_LINE() {std::cout << CLR_LINE << std::endl << std::flush;}

#ifndef VERSION
    #define VERSION "undefined vesion"
#endif

// A makro that will print the programms headline including name and version ID.
#define PRINT_HEADLINE() {                                                                                      \
    std::string formatted_title = std::string("****** ") + "DUST COLLIDER - " + VERSION + " " + SEP_LINE;       \
    formatted_title = formatted_title.substr(0, 85);                                                            \
    PRINT_SEP_LINE();                                                                                           \
    std::cout << formatted_title << std::endl << std::flush;                                                    \
    PRINT_SEP_LINE();                                                                                           \
}

// A makro that will print a title line with the specified title.
#define PRINT_TITLE(title) {                                                            \
    std::string formatted_title = std::string("****** ") + title + " " + SEP_LINE;      \
    formatted_title = formatted_title.substr(0, 85);                                    \
    std::cout << formatted_title << std::endl << std::flush;                            \
}

// A makro that will print a log message to console, with verbosity control.
#define PRINT_LOG(message, level) {                                             \
    if (VERBOSITY >= level) {                                                   \
        std::string formatted_message = std::string("> ") + message;            \
        std::cout << formatted_message << std::endl << std::flush;              \
    }                                                                           \
}

// A makro that will print an error message to console, without throwing an exception.
#define PRINT_ERROR(message) {                                                  \
    std::ostringstream oss; \
    oss << "ERROR: " << message << "\n      at " << __FILE__ << ":" << __LINE__;           \
    std::cout << oss.str() << std::endl << std::flush;                  \
}

/**
 * @brief Converts a duration in nanoseconds into a string of the form HHHH:MM:SS.mm.
 * 
 * @param duration: The duration that is to be converted into a string.
 * @param buffer: A buffer the resulting string is to be written into, needs to have length .
 */
void ns_to_time_string(const long duration, char* buffer, size_t buffer_size) {
    if (buffer_size < 14) {
        PRINT_ERROR("The supplied buffer was too small. Buffer needs to be at least 14 chars long.");
    }

    long remaining = duration;

    long hours = remaining / 3'600'000'000'000;
    remaining %= 3'600'000'000'000;

    long minutes = remaining / 60'000'000'000;
    remaining %= 60'000'000'000;

    double seconds = static_cast<double>(remaining) / 1'000'000'000.0;
    
    std::snprintf(buffer, buffer_size, "%04ld:%02ld:%05.02f\n", hours, minutes, seconds);
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