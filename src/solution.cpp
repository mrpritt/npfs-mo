/**
 * \file solution.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#include "solution.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "boost/multi_array.hpp"
using namespace boost;

#include "holes.hpp"

void PSolution::read(istream &in) {
  string line;
  auto j = 1u;
  while (getline(in,line)) {
    if (line[0] == '#') continue;
    istringstream iss(line);
    Time t;
    while (iss >> t)
      π[j++] = t + 1;
  }
}

void PSolution::write(ostream &out) {
  for (auto j = 1u; j <= n; ++j)
    fmt::print(out, "{} ", π[j] - 1);
  fmt::print(out, "\n");
}

Result PSolution::getResult(const Instance& I) {
  auto [ms, ft] = compute_ms_ft_mo(I);
  return Result{ms, ft, 0.0};
}

pair<Time, Time> PSolution::compute_ms_ft_mo(const Instance &I, vector<Time> &C) {
  assert(I.n == n && I.m == m && C.size() == m + 1);
  Time ft = 0, ms = 0;
  fill(C.begin(), C.end(), 0);
  for (unsigned j = 1; j != fbegin; ++j) {
    Time Cj = 0;
    const auto jb = π[j];
    for (unsigned i = 1; i <= m; ++i)
      if (I.p[jb][i] > 0)
        C[i] = Cj = max(C[i], Cj) + I.p[jb][i];
    ft += Cj;
    ms = max(ms, Cj);
  }
  return { ms, ft };
}

void PSolution::compute_completion_times(const Instance &I, multi_array<Time,2> &Ci, multi_array<Time,2> &Cj) const {
  for (unsigned j = 0; j < fbegin; ++j)
    Ci[0][j] = Cj[0][j] = 0;
  for (unsigned i = 0; i <= m; ++i)
    Ci[i][0] = Cj[i][0] = 0;

  for (unsigned k = 1; k < fbegin; ++k) {
    const auto jb = π[k];
    for (unsigned i = 1; i <= m; ++i) {
      Cj[i][k] = Cj[i - 1][k];
      Ci[i][k] = Ci[i][k - 1];
      if (I.p[jb][i] > 0)
	Ci[i][k] = Cj[i][k] = max(Ci[i][k], Cj[i][k]) + I.p[jb][i];
    }
  }
}

void PSolution::compute_completion_times(const Instance &I, multi_array<Time,2> &C) const {
  multi_array<Time,2> Ci(extents[m + 1][n + 1]);
  multi_array<Time,2> Cj(extents[m + 1][n + 1]);
  compute_completion_times(I, Ci, Cj);

  for (unsigned i = 1; i <= m; ++i)
    for (unsigned k = 1; k != fbegin; ++k) {
      const auto j = π[k];
      C[i][j] = Cj[i][k];
    }
}

string PSolution::to_string() const { return fmt::format("<{};of={}>", fmt::join(π.begin() + 1, π.begin() + fbegin, ","), of); }

bool PSolution::valid_permutation() const {
#if !defined(NDEBUG)
  vector<Job> σ = π;
  sort(σ.begin(), σ.end());
  for (auto j = 1u; j <= n; ++j)
    if (σ[j] != j)
      return false;
#endif
  return true;
}

std::pair<Time, Time> PSolution::evaluateNPSset(const Instance &I, bool smallest) const {
  vector<hlist> h(m);
  Time ms = 0, ft = 0;

  for (auto j = 1u; j <= n; ++j) {
    Time Ct = 0;
    for (auto i = 1u; i <= m; ++i) {
      const auto p = I.p[π[j]][i];
      if (p == 0)
        continue;
      auto e = smallest ? h[i - 1].smallest(p, Ct) : h[i - 1].earliest(p, Ct);
      const auto s = e->start();
      if (Ct <= s) {
        h[i - 1].reduce(e, p);
        Ct = s;
      } else {
        h[i - 1].cut(e, Ct, p);
      }
      Ct += p;
    }
    ms = max(ms, Ct);
    ft += Ct;
  }
  return {ms, ft};
}

void NPSolution::read(istream &in) {
  string line;
  auto i = 1u, j = 1u;
  while (getline(in,line)) {
    if (line[0] == '#') continue;
    istringstream iss(line);
    Time t;

    while (iss >> t) {
      π[i][j++] = t + 1;
      if (j > n) {
	j = 1;
	i++;
      }
    }
  }
}

void NPSolution::write(ostream &out) {
  for (auto i = 1u; i <= m; ++i) {
    for (auto j = 1u; j <= n; ++j)
      fmt::print(out, "{} ", π[i][j] - 1);
    fmt::print(out, "\n");
  }
}

pair<Time, Time> NPSolution::compute_ms_ft_mo(const Instance &I, vector<Time> &C) {
  assert(I.n == n && I.m == m && C.size() == m + 1);
  Time ms = 0;
  vector<Time> Cj(n + 1, 0);
  for (unsigned i = 1; i <= m; ++i) {
    Time Ci = 0;
    for (unsigned k = 1; k != fbegin; ++k) {
      const auto j = π[i][k];
      if (I.p[j][i] > 0)
        Ci = Cj[j] = max(Ci, Cj[j]) + I.p[j][i];
    }
    ms = max(ms, Ci);
  }
  return {ms, accumulate(Cj.begin() + 1, Cj.end(), 0)};
}

double NPSolution::computeJRI(const Instance& I) const {
  multi_array<Time,2> Cj(extents[m + 1][n + 1]);
  compute_completion_times(I, Cj);

  struct Event { Time t; Job j; };
  vector<vector<Event>> buffer(m + 1, vector<Event>());

  for(unsigned i = 1; i < m; ++i)
    for(unsigned k = 1; k <= n; ++k) {
      Job j = π[i][k];
      if (I.p[j][i] == 0) {
	if (i == 1) {
	  unsigned fm = I.firstOperation(j);
	  if (fm <= m) {
	    buffer[fm].push_back({0, j});
	  }
	}
	continue;
      }

      unsigned nm = I.nextOperation(i, j);
      if (nm == m + 1) continue;

      unsigned nji = k + 1;
      while (nji <= n && I.p[π[i][nji]][i] == 0) nji++;
      if (nji > n) {
	buffer[nm].push_back({Cj[I.lastOperation(j)][j], j});
	continue;
      }
      unsigned nj = π[i][nji];

      Time Sⱼ = Cj[nm][j] - I.p[j][nm];
      Time Cⱼ = std::min(Cj[i][nj] - I.p[nj][i], Sⱼ);

      if (Cⱼ < Sⱼ) {
	buffer[nm].push_back({Cⱼ, j});
      } else {
	buffer[nm].push_back({Sⱼ, j});
      }
    }

  double jri = 0.0;
  for(unsigned i = 2; i <= m; ++i) {
    sort(buffer[i].begin(), buffer[i].end(), [&](const Event& e, const Event & f) { return e.t < f.t || (e.t == f.t && jobIndex(i, e.j) < jobIndex(i, f.j)); });
    std::vector<Job> bv, πv;
    for (auto be : buffer[i])
      bv.push_back(be.j);
    for(unsigned k = 1; k <= n; ++k)
      if (I.p[π[i][k]][i] != 0) πv.push_back(π[i][k]);
    jri += kendall_tau(bv, πv);
  }
  return jri / ( n * (m - 1) );
}

unsigned NPSolution::computeBufferspace(const Instance& I) const {
  multi_array<Time,2> Cj(extents[m + 1][n + 1]);
  compute_completion_times(I, Cj);

  struct Event { Time t; unsigned i; int d;
    bool operator<(const Event& e) const { return t < e.t; }
  };
  std::vector<Event> e;

  for(unsigned i = 1; i < m; ++i)
    for(unsigned k = 1; k < n; ++k) {
      unsigned j = π[i][k];
      if (I.p[j][i] == 0) continue;

      unsigned nm = I.nextOperation(i, j);
      if (nm == m + 1) continue;

      unsigned nji = k + 1;
      while (nji <= n && I.p[π[i][nji]][i] == 0) nji++;
      if (nji > n) continue;
      unsigned nj = π[i][nji];

      Time Sⱼ = Cj[nm][j] - I.p[j][nm];
      Time Cⱼ = std::min(Cj[i][nj] - I.p[nj][i], Sⱼ);
      if (Cⱼ < Sⱼ) {
	e.push_back({Cⱼ, nm,  1});
	e.push_back({Sⱼ, nm, -1});
      }
    }

  unsigned mbs = 0;
  std::sort(e.begin(),e.end());
  unsigned bs = 0;
  Time lt = 0;
  for(auto ev : e) {
    assert(int(bs) + ev.d >= 0);
    bs += ev.d;
    if (ev.t != lt) {
      mbs = std::max(mbs,bs);
      lt = ev.t;
    }
  }
  mbs = std::max(mbs,bs);
  return mbs;
}

void NPSolution::compute_completion_times(const Instance &I, multi_array<Time,2> &Ci, multi_array<Time,2> &Cj, vector<unsigned>& last, vector<Time>& Cjb) const {
  for (unsigned j = 0; j != fbegin; ++j)
    Ci[0][j] = Cj[0][j] = 0;
  for (unsigned i = 0; i <= m; ++i)
    Ci[i][0] = Cj[i][0] = 0;

  for (unsigned i = 1; i <= m; ++i) {
    for (unsigned k = 1; k != fbegin; ++k) {
      const auto j = π[i][k];
      Cj[i][k] = Cjb[j];
      Ci[i][k] = Ci[i][k - 1];
      if (I.p[j][i] > 0) {
	Ci[i][k] = Cj[i][k] = Cjb[j] = max(Ci[i][k], Cj[i][k]) + I.p[j][i];
	last[j] = i;
      }
    }
  }
}

void NPSolution::compute_completion_times(const Instance &I, multi_array<Time,2> &C) const {
  multi_array<Time,2> Ci(extents[m + 1][n + 1]);
  multi_array<Time,2> Cj(extents[m + 1][n + 1]);
  vector<unsigned> last(n + 1, 0);
  vector<Time> Cjb(n + 1, 0);

  compute_completion_times(I, Ci, Cj, last, Cjb);

  for (unsigned i = 1; i <= m; ++i)
    for (unsigned k = 1; k != fbegin; ++k) {
      const auto j = π[i][k];
      C[i][j] = Cj[i][k];
    }
}
