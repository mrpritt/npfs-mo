/**
 * \file solution.hpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 *
 * Representation of solutions.
 */
#pragma once

#include <cstdint>
#include <limits>

#define FMT_HEADER_ONLY
#include "fmt/color.h"
#include "fmt/format.h"

#include "instance.hpp"

// result for reporting
struct Result {
  Time ms, ft;
  double tfound;
  std::string to_string() { return fmt::format("{} {} {}", ms, ft, tfound); }
};

struct Solution {
  Solution() : n(0), m(0), of(infinite_time) {}
  Solution(const Instance &I) : n(I.n), m(I.m), of(infinite_time) {}
  unsigned n, m;
  Time of; // primary objective (makespan or flowtime)

  Solution(Solution &&other) { this->swap(other); }

  Solution(const Solution &other) {
    n = other.n;
    m = other.m;
    of = other.of;
  }

  Solution &operator=(Solution other) {
    this->swap(other);
    return *this;
  }

  void swap(Solution &other) {
    using std::swap;
    swap(n, other.n);
    swap(m, other.m);
    swap(of, other.of);
  }
};

// solution for a permutation flow shop
struct PSolution : public Solution {
  typedef Solution Base;

  std::vector<Job> π; // permutation of jobs, 1-based
  unsigned fbegin;

  PSolution(const Instance &I) : Base(I), π(I.n + 1, 0), fbegin(I.n + 1) { std::iota(π.begin(), π.end(), 0); }

  PSolution(const Instance &I, const std::vector<Job> &π) : Base(I), π(π), fbegin(I.n + 1) {
    assert(π.size() == I.n + 1);
    assert(valid_permutation());
  }

  PSolution(PSolution &&other) { this->swap(other); }

  PSolution(const PSolution &other) : Base(other) {
    π = other.π;
    fbegin = other.fbegin;
  }

  PSolution &operator=(PSolution other) {
    this->swap(other);
    return *this;
  }

  void swap(PSolution &other) {
    Base::swap(other);
    using std::swap;
    π.swap(other.π);
    std::swap(fbegin, other.fbegin);
  }

  void read(std::istream &);
  void write(std::ostream &);

  // compute completion times `C`, update makespan, return flowtime
  std::pair<Time, Time> compute_ms_ft_mo(const Instance &I, std::vector<Time> &C);
  std::pair<Time, Time> compute_ms_ft_mo(const Instance &I) {
    std::vector<Time> C(I.m + 1, 0);
    return compute_ms_ft_mo(I, C);
  }
  std::pair<Time, Time> evaluateNPSset(const Instance &I, bool = false) const;

  Result getResult(const Instance &I);
  Time getMakespan(const Instance &I) { return compute_ms_ft_mo(I).first; }
  Time getFlowtime(const Instance &I) { return compute_ms_ft_mo(I).second; }

  std::string to_string() const;

  unsigned computeBufferspace(const Instance &I) const;

  void compute_completion_times(const Instance &, boost::multi_array<Time, 2> &, boost::multi_array<Time, 2> &) const;
  void compute_completion_times(const Instance &, boost::multi_array<Time, 2> &) const;

  bool isValid() const { return of < infinite_time; }

  bool valid_permutation() const;
};

// solution for a non-permutation flow shop
struct NPSolution : public Solution {
  typedef Solution Base;

  boost::multi_array<Job, 2> π; // permutations of jobs, i=1:m, j=1:n
  unsigned fbegin;              // by convention: fbegin,... all stored straight.

  NPSolution(const Instance &I) : Base(I), π(boost::extents[I.m + 1][I.n + 1]), fbegin(I.n + 1) {
    for (auto i = 1u; i <= I.m; ++i)
      for (auto j = 1u; j <= I.n; ++j)
        π[i][j] = j;
  }

  NPSolution(const Instance &I, const PSolution &P) : Base(Solution(P)), π(boost::extents[I.m + 1][I.n + 1]), fbegin(I.n + 1) {
    for (auto i = 1u; i <= I.m; ++i)
      for (auto j = 1u; j <= I.n; ++j)
        π[i][j] = P.π[j];
  }

  NPSolution(NPSolution &&other) { this->swap(other); }

  NPSolution(const NPSolution &other) : Base(other) {
    π.resize(boost::extents[other.π.shape()[0]][other.π.shape()[1]]);
    π = other.π;
    fbegin = other.fbegin;
  }

  NPSolution &operator=(NPSolution other) {
    this->swap(other);
    return *this;
  }

  void swap(NPSolution &other) {
    Base::swap(other);
    using std::swap;
    π.resize(boost::extents[other.π.shape()[0]][other.π.shape()[1]]);
    swap(π, other.π);
    swap(fbegin, other.fbegin);
  }

  void read(std::istream &);
  void write(std::ostream &);

  // compute completion times `C`, return makespan and flowtime
  std::pair<Time, Time> compute_ms_ft_mo(const Instance &I, std::vector<Time> &C);
  std::pair<Time, Time> compute_ms_ft_mo(const Instance &I) {
    std::vector<Time> C(I.m + 1, 0);
    return compute_ms_ft_mo(I, C);
  }

  Time getMakespan(const Instance &I) { return compute_ms_ft_mo(I).first; }
  Time getFlowtime(const Instance &I) { return compute_ms_ft_mo(I).second; }

  unsigned computeBufferspace(const Instance &I) const;
  double computeJRI(const Instance &I) const;
  unsigned jobIndex(unsigned i, Job j) const {
    unsigned k = 1;
    while (k <= n && π[i][k] != j)
      ++k;
    assert(k <= n);
    return k;
  }

  void compute_completion_times(const Instance &, boost::multi_array<Time, 2> &, boost::multi_array<Time, 2> &, std::vector<unsigned> &, std::vector<Time> &) const;
  void compute_completion_times(const Instance &, boost::multi_array<Time, 2> &) const;

  bool isValid() const { return of < infinite_time; }
};
