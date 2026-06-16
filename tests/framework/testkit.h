#pragma once

#include <cstdio>
#include <cmath>

/**
 * Minimal testing framework for the dust collider test suite.
 *
 * Usage
 * -----
 * Each test is a void function. At the end of main() call testkit::summary()
 * and return its value as the process exit code.
 *
 *   void test_foo() {
 *       CHECK(1 + 1 == 2);
 *       CHECK_APPROX(computed, expected, 1e-9);
 *   }
 *
 *   int main() {
 *       RUN_SUITE("foo", test_foo);
 *       return testkit::summary();
 *   }
 *
 * Define TESTKIT_NO_COLOR before including to suppress ANSI color codes.
 */

// ANSI color strings — empty when TESTKIT_NO_COLOR is defined
#ifndef TESTKIT_NO_COLOR
#  define _TC_RST  "\033[0m"
#  define _TC_SUIT "\033[1;36m"   // bold cyan  — suite headers
#  define _TC_PASS "\033[32m"     // green      — pass lines
#  define _TC_FAIL "\033[31m"     // red        — fail lines
#  define _TC_DIM  "\033[2m"      // dim        — file/line info
#else
#  define _TC_RST  ""
#  define _TC_SUIT ""
#  define _TC_PASS ""
#  define _TC_FAIL ""
#  define _TC_DIM  ""
#endif

namespace testkit {

inline int& _pass() { static int n = 0; return n; }
inline int& _fail() { static int n = 0; return n; }

inline void _check(bool ok, const char* expr, const char* file, int line) {
    if (ok) {
        _pass()++;
    } else {
        _fail()++;
        printf("    " _TC_FAIL "✗" _TC_RST "  %s\n"
               "    " _TC_DIM  "   at %s:%d" _TC_RST "\n",
               expr, file, line);
    }
}

inline bool _approx(double a, double b, double rel_tol) {
    double scale = 1.0 + std::fabs(a) + std::fabs(b);
    return std::fabs(a - b) <= rel_tol * scale;
}

inline int summary() {
    int total = _pass() + _fail();
    if (_fail() == 0)
        printf("\n  " _TC_PASS "✓  All %d checks passed." _TC_RST "\n\n", total);
    else
        printf("\n  " _TC_FAIL "✗  %d / %d checks FAILED." _TC_RST "\n\n", _fail(), total);
    return _fail() > 0 ? 1 : 0;
}

} // namespace testkit

// Assert that expr is true.
#define CHECK(expr) \
    testkit::_check(!!(expr), #expr, __FILE__, __LINE__)

// Assert that |a - b| <= tol * (1 + |a| + |b|).
#define CHECK_APPROX(a, b, tol) \
    testkit::_check( \
        testkit::_approx((double)(a), (double)(b), (tol)), \
        #a " ~ " #b, __FILE__, __LINE__)

// Run a named suite, print a header, and report per-suite pass/fail count.
#define RUN_SUITE(name, fn) \
    do { \
        printf("\n  " _TC_SUIT "▸ %s" _TC_RST "\n", (name)); \
        int _p0 = testkit::_pass(), _f0 = testkit::_fail(); \
        (fn)(); \
        int _n  = (testkit::_pass() + testkit::_fail()) - (_p0 + _f0); \
        int _nf = testkit::_fail() - _f0; \
        if (_nf == 0) \
            printf("  " _TC_PASS "✓  %d checks" _TC_RST "\n", _n); \
        else \
            printf("  " _TC_FAIL "✗  %d / %d checks failed" _TC_RST "\n", _nf, _n); \
    } while (0)
