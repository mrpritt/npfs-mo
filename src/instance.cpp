/**
 * \file instance.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#include "instance.hpp"

#include <set>
#include <string>
using namespace std;

#define FMT_HEADER_ONLY
#include "fmt/color.h"
#include "fmt/format.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h"

Instance::Instance(istream &in) { read_hn(in); }

void Instance::read_hn(istream &in) {
  string t1, t2, t3, t4, tj;
  in >> t1 >> m >> t2 >> n >> t3 >> t4 >> r;
  p.resize(boost::extents[n + 1][m + 1]);
  assert(t1 == "numberMachines" && t2 == "numberJobs" && t3 == "missing" && t4 == "Operations");
  for (auto i = 1u; i <= m; ++i)
    for (auto j = 1u; j <= n; ++j) {
      in >> tj >> p[j][i];
      assert(tj == fmt::format("t_{}_{}", i - 1, j - 1));
    }
  compute_auxiliary_data();
}

void Instance::reverse() {
  for (auto i = 1u; i != m / 2; ++i)
    for (auto j = 1u; j != n; ++j)
      std::swap(p[j][i], p[j][m - i + 1]);
  compute_auxiliary_data();
}

void Instance::compute_auxiliary_data() {
  oeff = 0;
  for (unsigned j = 1; j <= n; ++j)
    for (unsigned i = 1; i <= m; ++i)
      if (p[j][i] != 0)
        ++oeff;
}

vector<Time> Instance::totalTimes() const {
  vector<Time> T(n + 1, 0);
  for (unsigned j = 1; j <= n; ++j)
    for (unsigned i = 1; i <= m; ++i)
      T[j] += p[j][i];
  return T;
}

Time Instance::totalTime() const {
  Time T = 0;
  for (unsigned j = 1; j <= n; ++j)
    for (unsigned i = 1; i <= m; ++i)
      T += p[j][i];
  return T;
}

unsigned Instance::numPseudojobs() const {
  unsigned pj = 0;
  for (unsigned j = 1; j <= n; ++j) {
    bool missing = true;
    for (unsigned i = 1; i <= m; ++i) {
      if (missing && p[j][i] != 0)
        ++pj;
      missing = p[j][i] == 0;
    }
  }
  return pj;
}

void compute_inverse(const vector<Job> &π, vector<Job> &π⁻) {
  assert(π.size() > 0);
  π⁻.resize(*max_element(π.begin(), π.end()) + 1);
  for (unsigned k = 1, ke = π.size(); k != ke; ++k)
    π⁻[π[k]] = k;
}

unsigned kendall_tau(const vector<Job> &π, const vector<Job> &σ) {
  assert(set(π.begin(), π.end()) == set(σ.begin(), σ.end()));
  if (π.size() == 0)
    return 0;

  vector<Job> π⁻, σ⁻;
  compute_inverse(π, π⁻);
  compute_inverse(σ, σ⁻);

  const unsigned ke = π.size();
  unsigned kt = 0;
  for (unsigned k1 = 1; k1 != ke; ++k1) {
    const Job j1 = π[k1];
    for (unsigned k2 = k1 + 1; k2 != ke; ++k2) {
      const Job j2 = π[k2];
      if ((π⁻[j1] < π⁻[j2] && σ⁻[j1] > σ⁻[j2]) || (π⁻[j1] > π⁻[j2] && σ⁻[j1] < σ⁻[j2]))
        kt++;
    }
  }
  return kt;
}
