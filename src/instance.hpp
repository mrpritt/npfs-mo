/**
 * \file instance.hpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#pragma once

#include <iostream>
#include <vector>

#include "boost/multi_array.hpp"

using Job = uint16_t;
using Time = uint32_t;
const Time infinite_time = std::numeric_limits<int32_t>::max();
const Time uinf = std::numeric_limits<unsigned>::max();

struct Instance {
  unsigned n;                    // number of jobs
  unsigned m;                    // number of machines
  double r;                      // missing operations rate
  unsigned oeff;                 // number of effective operations
  boost::multi_array<Time, 2> p; // processing times, job j=1:n, machine i=1:m+1

  Instance(unsigned n = 0, unsigned m = 0) : n(n), m(m), p(boost::extents[n + 1][m + 2]) {}

  // create from input stream
  Instance(std::istream &in);

  // read from stream (Henneberg & Neufeld's format)
  void read_hn(std::istream &in);

  // reverse job order
  void reverse();

  void compute_auxiliary_data();

  std::vector<Time> totalTimes() const;
  Time totalTime() const;
  double avgTime() const { return totalTime() / static_cast<double>(numEffectiveOperations()); }
  unsigned numOperations() const { return n * m; }
  unsigned numEffectiveOperations() const { return oeff; }
  unsigned numPseudojobs() const;

  unsigned firstOperation(Job j) const { unsigned o = 1; while (o <= m && p[j][o] == 0) ++o; return o; }
  unsigned nextOperation(unsigned i, Job j) const { unsigned o = i + 1; while (o <= m && p[j][o] == 0) ++o; return o; }
  unsigned lastOperation(Job j) const { unsigned o = m; while (p[j][o] == 0 && o>0) --o; return o; }
};

unsigned kendall_tau(const std::vector<Job>&, const std::vector<Job>&);
