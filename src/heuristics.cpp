/**
 * \file heuristics.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#include "heuristics.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "helpers.hpp"
#include "random.hpp"

void EPSolution::clear() { fbegin = 1; }

void EPSolution::totalTimeOrder() {
  auto T = I.totalTimes();
  sort(π.begin() + 1, π.begin() + fbegin, [&T](Job i, Job j) { return T[i] > T[j] || (T[i] == T[j] && i < j); });
}

void EPSolution::update_heads(unsigned kb, unsigned ke) {
  for (unsigned k = kb; k != ke; ++k) {
    Time Ck = 0;
    for (unsigned i = 1; i <= m; ++i)
      if (I.p[π[k]][i] == 0)
        h[i][k] = h[i][k - 1];
      else
        h[i][k] = Ck = std::max(Ck, h[i][k - 1]) + I.p[π[k]][i];
  }
}

void EPSolution::update_heads_flowtimes(unsigned kb, unsigned ke, vector<Time> &ft) {
  for (unsigned k = kb; k != ke; ++k) {
    Time Ck = 0;
    for (unsigned i = 1; i <= m; ++i) {
      if (I.p[π[k]][i] == 0)
        h[i][k] = h[i][k - 1];
      else
        h[i][k] = Ck = std::max(Ck, h[i][k - 1]) + I.p[π[k]][i];
      ft[k] = ft[k - 1] + Ck;
    }
  }
}

void EPSolution::update_tails(unsigned kb, unsigned ke) {
  for (unsigned k = kb; k != ke; ++k) {
    Time Ck = 0;
    for (unsigned i = m; i >= 1; --i)
      if (I.p[π[fbegin - k]][i] == 0)
        t[i][k] = t[i][k - 1];
      else
        t[i][k] = Ck = std::max(Ck, t[i][k - 1]) + I.p[π[fbegin - k]][i];
  }
}

bool EPSolution::makespan_valid(Time Cm) { return Cm == compute_ms_ft_mo(I).first; }
bool EPSolution::flowtime_valid(Time Cf) { return Cf == compute_ms_ft_mo(I).second; }

void EPSolution::insert_all_ms() {
  update_heads(1, fbegin);
  update_tails(1, fbegin);

  Time Cm = 0;
  for (auto πend = π.size(); fbegin != πend;) {
    const unsigned jb = π[fbegin];

    unsigned bp = 0, Ip = uinf;
    Cm = infinite_time;
    for (unsigned k = fbegin; k >= 1; --k) {
      Time Cmaxk = 0, Cj = 0, Cjk = 0, Ik = 0, Ci = 0;
      for (unsigned i = 1; i <= m; ++i) {
        if (I.p[jb][i] > 0)
          Ci = Cj = std::max(Cj, h[i][k - 1]) + I.p[jb][i];
        else
          Ci = h[i][k - 1];
        Cmaxk = std::max(Cmaxk, Ci + t[i][fbegin - k]);
        if (Cmaxk > Cm)
          break;
        if (k < fbegin) {
          if (I.p[π[k]][i] > 0) {
            Cjk = std::max(Ci, Cjk) + I.p[π[k]][i];
            assert(Cjk >= h[i][k]);
            Ik += Cjk - h[i][k];
          } else {
            assert(Ci >= h[i][k]);
            Ik += Ci - h[i][k];
          }
        } else if (I.p[jb][i] > 0) {
          assert(Cj >= h[i][k - 1]);
          Ik += Cj - h[i][k - 1];
        }
      }
      if (Cmaxk < Cm || (Cmaxk == Cm && Ik < Ip)) {

        Cm = Cmaxk;
        bp = k;
        Ip = Ik;
      }
    }

    assert(bp <= fbegin);
    if (bp < fbegin)
      rotate(π.begin() + bp, π.begin() + fbegin, π.begin() + fbegin + 1);
    fbegin++;

    if (fbegin != πend) {
      update_heads(bp, fbegin);
      assert(bp + 1 <= fbegin);
      update_tails(fbegin - bp, fbegin);
    }
    assert(makespan_valid(Cm));
  }
  assert(makespan_valid(Cm));
  of = Cm;
}

void EPSolution::insert_all_ft() {
  vector<Time> ftk(I.n + 1, 0);
  update_heads_flowtimes(1, fbegin, ftk);
  update_tails(1, fbegin);

  Time Cm = 0, Cf = 0;
  for (auto πend = π.size(); fbegin != πend;) {
    const unsigned jb = π[fbegin];

    unsigned bp = 0;
    Cm = Cf = infinite_time;

    for (unsigned k = fbegin; k >= 1; --k) {
      Time fk = ftk[k - 1];
      vector<Time> C(m + 1, 0);

      Time Cj = 0;
      for (unsigned i = 1; i <= m; ++i) {
        if (I.p[jb][i] > 0)
          C[i] = Cj = std::max(Cj, h[i][k - 1]) + I.p[jb][i];
        else
          C[i] = h[i][k - 1];
      }
      fk += Cj;
      if (fk > Cf)
        break;

      for (unsigned l = k; l != fbegin; ++l) {
        Cj = 0;
        for (unsigned i = 1; i <= m; ++i) {
          if (I.p[π[l]][i] > 0)
            C[i] = Cj = std::max(Cj, C[i]) + I.p[π[l]][i];
        }
        fk += Cj;
        if (fk > Cf)
          break;
      }

      Time Cmaxk = *max_element(C.begin(), C.end());
      if (fk < Cf || (fk == Cf && Cmaxk < Cm)) {
        Cf = fk;
        Cm = Cmaxk;
        bp = k;
      }
    }

    assert(bp <= fbegin);
    if (bp < fbegin)
      rotate(π.begin() + bp, π.begin() + fbegin, π.begin() + fbegin + 1);
    fbegin++;

    if (fbegin != πend) {
      update_heads_flowtimes(bp, fbegin, ftk);
      assert(bp + 1 <= fbegin);
      update_tails(fbegin - bp, fbegin);
    }

    assert(flowtime_valid(Cf));
  }
  assert(flowtime_valid(Cf));
  of = Cf;
}

void EPSolution::shuffle_free() { shuffle(π.begin() + fbegin, π.end(), rng); }

void EPSolution::remove(unsigned k) {
  rotate(π.begin() + k, π.begin() + k + 1, π.begin() + fbegin);
  fbegin--;
}

bool EPSolution::shift_step() {
  Time of_ = of;
  for (auto j = 1u; j <= n; ++j) {
    remove(j);
    insert_all();
    store_so();
  }
  assert(of <= of_);
  return of < of_;
}

unsigned EPSolution::shift_ls() {
  store_so();
  unsigned steps = 0;
  while (shift_step()) {
    tfound = run::elapsed();
    steps++;
    vprint(3, "{}\n", of);
  }
  return steps;
}

void EPSolution::iga_perturb(unsigned dc) {
  assert(fbegin - 1 > dc);
  auto is = sample_floyd(dc, fbegin - 1);
  transform(is.begin(), is.end(), is.begin(), [](auto e) { return e + 1; });
  is.push_back(fbegin);

  for (auto i = 0u; i != dc; ++i)
    rotate(π.begin() + is[i] - i, π.begin() + is[i] + 1, π.begin() + is[i + 1]);
  fbegin -= dc;
  shuffle_free();
  insert_all();
}

unsigned EPSolution::iga(const IGAOptions &opt) {

  unsigned steps = 0;
  SSolution bs{π, of, run::elapsed()};

  vprint(2, "IGA starts {} {}\n", of, so.of);
  double last_report = run::elapsed();
  while (!opt.stop(steps)) {
    vprint(3, "IGA has {} {}\n", of, so.of);
    SSolution ps{π, of, 0};
    iga_perturb(opt.dc);
    shift_ls();
    if (of < bs.of) {
      bs = SSolution{π, of, run::elapsed()};
      Time ft = getFlowtime(I);
      vprint(2, "* {:4.1f} {} {} {}\n", run::elapsed(), of, ft, steps);
      store_so();
      last_report = run::elapsed();
    } else if (!(of < ps.of || getRandom() < exp(-double(of - ps.of) / opt.T))) {
      π = ps.π;
      of = ps.of;
    }
    steps++;
    if (verbose(2) && run::elapsed() > last_report + 1) {
      last_report = run::elapsed();
      Time ft = getFlowtime(I);
      fmt::print(". {:4.1f} {} {} {}\n", run::elapsed(), of, ft, steps);
      store_so();
    }
  }
  π = bs.π;
  of = bs.of;
  tfound = bs.tfound;
  return steps;
}

Result EPSolution::getResultPO() {
  auto [ms, ft] = compute_ms_ft_mo(I);
  return Result{ms, ft, tfound};
}

Result EPSolution::getResultSO() {
  PSolution Sf{I, so.π};
  auto [ms, ft] = Sf.compute_ms_ft_mo(I);
  return Result{ms, ft, so.tfound};
}

bool IGAOptions::stop(unsigned steps) const {
  if (iterlimit > 0 && int(steps) > iterlimit)
    return true;
  if (timelimit > 0 && run::elapsed() > timelimit)
    return true;
  return false;
}

void ENPSolution::compute_ρ() {
  for (auto i = 1u; i <= I.m; ++i)
    for (auto k = 1u; k <= I.n; ++k)
      ρ[i][π[i][k]] = k;
}

void ENPSolution::clear() { fbegin = 1; }

struct NPMove {
  unsigned k;
  int i;
  Time Csum, Cmax;

  void update(const NPMove &cm) {
    if (cm.Csum < Csum || (cm.Csum == Csum && cm.Cmax < Cmax)) {
      *this = cm;
    }
  }
};

void ENPSolution::insert_all_ft() {
  NPMove bm{0, 0, infinite_time, infinite_time};

  for (auto πend = π.shape()[1]; fbegin != πend;) {
    for (unsigned k = fbegin++; k > 0; --k) {
      auto [Cmax, Csum] = compute_ms_ft_mo(I);
      bm.update(NPMove{k, 0, Csum, Cmax});

      if (k == 1)
        break;

      for (unsigned l = m; l != 1; --l) {
        ::swap(π[l][k - 1], π[l][k]);
        auto [Cmax, Csum] = compute_ms_ft_mo(I);
        bm.update(NPMove{k, -int(l), Csum, Cmax});
      }
      ::swap(π[1][k - 1], π[1][k]);
    }
    for (unsigned i = 1; i <= m; ++i)
      rotate(π[i].begin() + 1, π[i].begin() + 2, π[i].begin() + fbegin);
    fbegin--;

    for (unsigned k = fbegin++; k > 1; --k) {
      for (unsigned l = 1; l < m; ++l) {
        ::swap(π[l][k - 1], π[l][k]);
        auto [Cmax, Csum] = compute_ms_ft_mo(I);
        bm.update(NPMove{k, int(l + 1), Csum, Cmax});
      }
      ::swap(π[m][k - 1], π[m][k]);
    }
    for (unsigned i = 1; i <= m; ++i)
      rotate(π[i].begin() + 1, π[i].begin() + 2, π[i].begin() + fbegin);
    fbegin--;

    if (bm.i < 0) {
      for (auto i = 1u; i < unsigned(-bm.i); ++i) {
        rotate(π[i].begin() + bm.k, π[i].begin() + fbegin, π[i].begin() + fbegin + 1);
        for (auto k = bm.k; k != fbegin + 1; ++k)
          ρ[i][π[i][k]] = k;
      }
      for (auto i = -unsigned(bm.i); i <= m; ++i) {
        rotate(π[i].begin() + bm.k - 1, π[i].begin() + fbegin, π[i].begin() + fbegin + 1);
        for (auto k = bm.k - 1; k != fbegin + 1; ++k)
          ρ[i][π[i][k]] = k;
      }
    } else if (bm.i == 0) {
      for (auto i = 1u; i <= I.m; ++i) {
        rotate(π[i].begin() + bm.k, π[i].begin() + fbegin, π[i].begin() + fbegin + 1);
        for (auto k = bm.k; k != fbegin + 1; ++k)
          ρ[i][π[i][k]] = k;
      }
    } else {
      for (auto i = 1u; i < unsigned(bm.i); ++i) {
        rotate(π[i].begin() + bm.k - 1, π[i].begin() + fbegin, π[i].begin() + fbegin + 1);
        for (auto k = bm.k - 1; k != fbegin + 1; ++k)
          ρ[i][π[i][k]] = k;
      }
      for (auto i = unsigned(bm.i); i <= m; ++i) {
        rotate(π[i].begin() + bm.k, π[i].begin() + fbegin, π[i].begin() + fbegin + 1);
        for (auto k = bm.k; k != fbegin + 1; ++k)
          ρ[i][π[i][k]] = k;
      }
    }
    of = bm.Csum;
    fbegin++;
  }
}

void ENPSolution::remove(Job j) {
  fbegin--;
  for (auto i = 1u; i <= I.m; ++i) {
    for (auto ji = ρ[i][j]; ji != fbegin; ++ji) {
      π[i][ji] = π[i][ji + 1];
      ρ[i][π[i][ji]] = ji;
    }
    π[i][fbegin] = j;
    ρ[i][j] = fbegin;
  }
}

bool ENPSolution::shift_step() {
  Time of_ = of;
  push_po();
  for (auto j = 1u; j <= n; ++j) {
    remove(j);
    insert_all();
    store_so();
  }
  if (of >= of_)
    pop_po();
  assert(of <= of_);
  return of < of_;
}

unsigned ENPSolution::shift_ls() {
  store_so();
  unsigned steps = 0;
  while (shift_step()) {
    tfound = run::elapsed();
    steps++;
    vprint(3, "{}\n", of);
  }
  return steps;
}

Result ENPSolution::getResultPO() {
  auto [ms, ft] = compute_ms_ft_mo(I);
  return Result{ms, ft, tfound};
}
