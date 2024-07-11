/**
 * \file random.cpp
 *   \author Marcus Ritt <mrpritt@inf.ufrgs.br>
 */
#include "random.hpp"

#include <time.h>

#include <fstream>
using namespace std;

mt19937 rng;

unsigned setupRandom(unsigned seed) {
  if (seed == 0) {
    seed = time(0);
    ifstream f("/dev/urandom");
    if (f.good()) {
      f.read((char *)(&seed), sizeof(unsigned int));
    }
  }
  rng.seed(seed);
  srand48(seed);
  return seed;
}
