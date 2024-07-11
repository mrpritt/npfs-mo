/**
 * \file heuristics.hpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#pragma once

#include "logging.hpp"
#include "solution.hpp"

struct IGAOptions {
  unsigned dc;
  double timelimit;
  int iterlimit;
  double T;
  double alpha;

  IGAOptions() : dc(8), timelimit(30), iterlimit(5000), T(0), alpha(0.2353) {}

  bool stop(unsigned steps) const;
};

struct SSolution {
  std::vector<Job> π;
  Time of;
  double tfound;
  SSolution() : of(infinite_time) {}
  SSolution(std::vector<Job> π, Time of, double tfound) : π(π), of(of), tfound(tfound) {}
};

struct EPSolution : public PSolution {
  using Base = PSolution;
  const Instance &I;
  boost::multi_array<Time, 2> h, t;
  bool of_makespan;
  double tfound;
  SSolution so;

  EPSolution(const Instance &I) : Base(I), I(I), h(boost::extents[m + 1][n + 1]), t(boost::extents[m + 2][n + 1]), of_makespan(true) {}
  EPSolution(const Instance &I, const Base &S) : Base(S), I(I), of_makespan(true) {}

  void update_heads(unsigned, unsigned);
  void update_heads_flowtimes(unsigned, unsigned, std::vector<Time>&);
  void update_tails(unsigned, unsigned);
  bool makespan_valid(Time);
  bool flowtime_valid(Time);
  Result getResultPO();
  Result getResultSO();

  void clear();
  void store_so() {
    auto [ms, ft] = compute_ms_ft_mo(I);
    auto sof = of_makespan ? ft : ms;
    if (sof < so.of)
      so = SSolution{π, sof, run::elapsed()};
  }
  void totalTimeOrder();
  void insert_all() { if (of_makespan) insert_all_ms(); else insert_all_ft(); }
  void insert_all_ms();
  void insert_all_ft();
  void shuffle_free();
  void remove(unsigned);
  bool shift_step();
  unsigned shift_ls();
  void iga_perturb(unsigned);
  unsigned iga(const IGAOptions &);
};

struct ENPSolution : public NPSolution {
  using Base = NPSolution;
  const Instance &I;
  double tfound;

  boost::multi_array<unsigned, 2> ρ;
  ENPSolution *S₀;

  ENPSolution(const Instance &I) : Base(I), I(I), tfound(0.0), ρ(boost::extents[m + 1][n + 1]), S₀(nullptr) { compute_ρ(); }
  ENPSolution(const Instance &I, const Base &S) : Base(S), I(I), tfound(0.0), ρ(boost::extents[m + 1][n + 1]), S₀(nullptr) { compute_ρ(); }
  ENPSolution(const Instance &I, const EPSolution &S) : Base(I,PSolution(S)), I(I), tfound(0.0), ρ(boost::extents[m + 1][n + 1]), S₀(nullptr) { compute_ρ(); }

  ENPSolution(ENPSolution &&other) : Base(other), I(other.I) { this->swap(other); }

  ENPSolution(const ENPSolution &other) : Base(other), I(other.I), tfound(other.tfound), S₀(nullptr) {
    ρ.resize(boost::extents[other.ρ.shape()[0]][other.ρ.shape()[1]]);
    ρ = other.ρ;
  }

  ENPSolution &operator=(ENPSolution other) {
    this->swap(other);
    return *this;
  }

  void swap(ENPSolution &other) {
    Base::swap(other);
    using std::swap;
    swap(tfound, other.tfound);
    ρ.resize(boost::extents[other.ρ.shape()[0]][other.ρ.shape()[1]]);
    swap(ρ, other.ρ);
  }

  void compute_ρ();

  Result getResultPO();
  void store_so() {}
  void push_po() { S₀ = new ENPSolution(*this); }
  void pop_po() { *this = *S₀; delete S₀; S₀ = nullptr; }

  void clear();
  void insert_all() { insert_all_ft(); }
  void insert_all_ft();
  void remove(Job);
  bool shift_step();
  unsigned shift_ls();
};
