#pragma once

#include <iostream>
#include <string>
#include <format>
#include <utility>

/**
 * @brief Contains logging functionality.
 * 
 * All output is written to std::cout.
 * When the code was compiled in the test build (-DTEST) all regular output is disabled.
 */
namespace Logger {
    // Width of a separator/title rule and the prefix used for header lines.
    inline const std::string SEP_LINE       = "************************************************************";
    inline const std::string HEADER_INDENT  = "******";

    /**
     * @brief Writes a single line to the output stream.
     *
     * @param content The string that is to be printed.
     */
    inline void writeLn(const std::string& content) {
        #ifndef TEST
            std::cout << content << std::endl << std::flush;
        #endif
    }

    /**
     * @brief Prints a boxed title (rule / title line / rule).
     * 
     * If the title is too large for the box it will still be printed in full,
     * the box formatting will just not look as good tho.
     */
    inline void title(const std::string& title) {
        std::string title_line = (HEADER_INDENT + " " + title);

        // Ensure that the title is not cut as it may contain critical information
        if (title_line.size() <= SEP_LINE.size()) {
            title_line = (title_line + " " + SEP_LINE).substr(0, SEP_LINE.size());
        }

        writeLn(SEP_LINE);
        writeLn(title_line);
        writeLn(SEP_LINE);
    }

    /**
     * @brief Prints a single header line.
     */
    inline void header(const std::string& header) {
        std::string header_line = (HEADER_INDENT + " " + header + " " + SEP_LINE).substr(0, SEP_LINE.size());

        writeLn(header_line);
    }

    /**
     * @brief Prints an empty line.
     */
    inline void lineBreak() {
        writeLn("");
    }

    /**
     * @brief Prints a seperator
     */
    inline void seperator() {
        writeLn(SEP_LINE);
    }

    /**
     * @brief Logs an info line, describing the current action.
     * 
     * The function uses the std::format syntax.
     *
     * @param  fmt  Compile-time std::format format string.
     * @param  args Values substituted into the `{}` placeholders.
     */
    template <class... Args>
    inline void log(std::format_string<Args...> fmt, Args&&... args) {
        writeLn("[>] " + std::format(fmt, std::forward<Args>(args)...));
    }

    /**
     * @brief Logs a warning line.
     *
     * The function uses the std::format syntax.
     *
     * @param  fmt  Compile-time std::format format string.
     * @param  args Values substituted into the `{}` placeholders.
     */
    template <class... Args>
    inline void warn(std::format_string<Args...> fmt, Args&&... args) {
        writeLn("[W] " + std::format(fmt, std::forward<Args>(args)...));
    }

    /**
     * @brief Logs an error line, prefixed with "[E]".
     *
     * The function uses the std::format syntax.
     * 
     * @param  fmt  Compile-time std::format format string.
     * @param  args Values substituted into the `{}` placeholders.
     */
    template <class... Args>
    inline void error(std::format_string<Args...> fmt, Args&&... args) {
        writeLn("[E] " + std::format(fmt, std::forward<Args>(args)...));
    }

    /**
     * @brief Prints formatted content.
     *
     * The function uses the std::format syntax.
     * 
     * @param  fmt  Compile-time std::format format string.
     * @param  args Values substituted into the `{}` placeholders.
     */
    template <class... Args>
    inline void print(std::format_string<Args...> fmt, Args&&... args) {
        writeLn(std::format(fmt, std::forward<Args>(args)...));
    }
}
