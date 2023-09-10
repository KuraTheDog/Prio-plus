#ifndef PRIOBASEUTILS_H
#define PRIOBASEUTILS_H

#include <emp-tool/emp-tool.h>  // for timing

// shorthand for timing functions
#define clock_start emp::clock_start
inline float sec_from(time_point<high_resolution_clock> start) {
  return (((float)emp::time_from(start)) / CLOCKS_PER_SEC);
}

[[noreturn]] inline void error_exit(const char* const msg) {
  perror(msg);
  exit(EXIT_FAILURE);
}

// Floor log, on nonnegative integers only. Ceiling
// Can also do math.h or cmath
// Note Log2(0) = 1
#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll(((X) | 1))))

inline bool isPowerOfTwo(const size_t n) {
  return n && !(n & (n - 1));
}

// todo: possibly connection code, from utils test connect and others?

#endif
