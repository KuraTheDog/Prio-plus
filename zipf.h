#ifndef ZIPF_H
#define ZIPF_H

#include <algorithm>
#include <iostream>
#include <random>
#include <set>
// #include <vector>

/* Currently returns values from 1 to support
Weight of value n is 1/n^exponent
*/
class ZipF {
  const uint64_t support;   // How many possible values
  const double exponent;        // Weight. 1 = base. Larger = more skewed

  const double normalizer;      // 1 / total weight

  std::random_device rd;
  std::mt19937 gen;

  // TODO: Future work. |values| = support random values up to bits large.
  // Requires sampling S distinct random items out of 2^bits
  // Implemnting is overkill for now.
  // std::vector<unsigned int> const values;

  // Sum 1/n^exponent from 1 to support
  double compute_normalizer(const uint64_t support, const double exponent) {
    double ans = 0.0;
    for (uint64_t i = 1; i <= support; i++) {
      ans += 1.0 / pow((double) i, exponent);
    }
    return 1.0 / ans;
  }

public:

  ZipF(const uint64_t support, const double exponent)
  : support(support)
  , exponent(exponent)
  , normalizer(compute_normalizer(support, exponent))
  , gen(rd())
  {
  }

  ~ZipF() {}

  uint64_t sample() {
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double z = dis(gen);
    // double z = rand() / RAND_MAX;  // Low set of possible values, not optimal
    // std::cout << "Sample z = " << z << "\n";

    double sum_prob = 0;
    for (uint64_t i = 1; i <= support; i++) {
      sum_prob += normalizer / pow((double) i, exponent);
      // std::cout << "  sum prob(1 to " << i << ") = " << sum_prob << "\n";
      if (sum_prob >= z) {
        // std::cout << "returning " << i << std::endl;
        return i;
      }
    }

    // std::cout << "returning 0" << std::endl;
    return 0;  // Should't really happen
  }
};

#endif
