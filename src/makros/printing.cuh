#pragma once

///////////////////////// FORMAT /////////////////////////

#define SEP_LINE    "*************************************************************************************"
#define CLR_LINE    "                                                                                     "

///////////////////////// PRINTING /////////////////////////

#define PRINT_SEP_LINE() {std::cout << SEP_LINE << std::endl << std::flush;}
#define PRINT_CLR_LINE() {std::cout << CLR_LINE << std::endl << std::flush;}

#ifndef VERSION
    #define VERSION "undefined vesion"
#endif

// A makro that will print a log message to console, with verbosity control.
#define PRINT_LOG(message, level) {                                             \
    if (VERBOSITY >= level) {                                                   \
        std::string formatted_message = std::string("> ") + message;            \
        std::cout << formatted_message << std::endl << std::flush;              \
    }                                                                           \
}

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
