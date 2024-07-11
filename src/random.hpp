/**
 * \file random.hpp
 *   \author Marcus Ritt <mrpritt@inf.ufrgs.br>
 *
 * Random number generator.
 */
#pragma once

#include <random>
#include <set>
#include <vector>

extern std::mt19937 rng;

unsigned setupRandom(unsigned seed = 0);

inline double getRandom() {
  std::uniform_real_distribution<> U;
  return U(rng);
}

inline std::vector<unsigned> sample_floyd(unsigned k, unsigned N) {
  std::uniform_int_distribution<> dis;
  std::set<unsigned> S;
  for (unsigned u = N - k; u != N; u++) {
    auto sample = dis(rng) % u;
    S.insert(S.find(sample) == S.end() ? sample : u);
  }
  return std::vector<unsigned>(S.begin(), S.end());
}
