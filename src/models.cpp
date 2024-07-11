/**
 * \file models.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#include "models.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "logging.hpp"

map<IloAlgorithm::Status, string> status_string = {{IloAlgorithm::Unknown, "Unknown"}, {IloAlgorithm::Feasible, "Feasible"}, {IloAlgorithm::Optimal, "Optimal"}, {IloAlgorithm::Infeasible, "Infeasible"}, {IloAlgorithm::Unbounded, "Unbounded"}, {IloAlgorithm::InfeasibleOrUnbounded, "InfeasibleOrUnbounded"}, {IloAlgorithm::Error, "Error"}};

string to_string(IloAlgorithm::Status status) { return status_string[status]; }

void setupSolver(IloCplex solver, const ModelOptions &m) {
  IloEnv env = solver.getEnv();
  if (!verbose(2)) {
    solver.setOut(env.getNullStream());
    solver.setWarning(env.getNullStream());
  }
  solver.setParam(IloCplex::Param::Threads, 1);
  solver.setParam(IloCplex::Param::TimeLimit, m.timelimit);
}

void MPFSMO::addVars() {
  const auto Csum = I.totalTime();
  for (unsigned i = 0u; i != I.m; ++i)
    for (unsigned j = 0u; j != I.n; ++j) {
      C.add(IloNumVar(env, I.p[j + 1][i + 1], IloInfinity, ILOFLOAT));
      C[C.getSize() - 1].setName(fmt::format("C[{},{}]", i + 1, j + 1).c_str());
      obj.setLinearCoef(C[C.getSize() - 1], 1.0 / (Csum * Csum));
    }
  idx_c = C.getSize();
  C.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
  C[idx_c].setName("Csum");
  obj.setLinearCoef(C[idx_c], 1.0);

  for (unsigned i = 0u; i != I.m; ++i)
    for (unsigned j1 = 0; j1 != I.n; ++j1)
      for (unsigned j2 = 0; j2 != I.n; ++j2)
        if (j1 != j2) {
          x.add(IloNumVar(env, 0.0, 1.0, ILOINT));
          x[x.getSize() - 1].setName(fmt::format("x[{},{},{}]", i + 1, j1 + 1, j2 + 1).c_str());
          idx[{i, j1, j2}] = x.getSize() - 1;
        }
}

void MPFSMO::addCompletion(Time ub_cmax) {
  IloRangeArray completion(env);
  int back = -1;
  for (unsigned i = 0; i != I.m; ++i) {
    unsigned mach = i;
    for (unsigned j = 0; j != I.n; ++j) {
      Time pij = I.p[j + 1][i + 1];
      if (pij == 0)
        continue;
      const double Mij = ub_cmax;

      for (unsigned j2 = 0; j2 != I.n; j2++) {
        if (j2 == j)
          continue;
        completion.add(IloRange(env, -Mij + pij, IloInfinity, fmt::format("C{}{}a", i + 1, j + 1).c_str()));
        back++;
        completion[back].setLinearCoef(C[i * I.n + j], 1.0);
        completion[back].setLinearCoef(C[i * I.n + j2], -1.0);
        completion[back].setLinearCoef(x[idx[{mach, j2, j}]], -Mij);
      }

      unsigned i2 = i;
      while (i2 > 0) {
        i2--;
        if (I.p[j + 1][i2 + 1] > 0) {
          completion.add(IloRange(env, pij, IloInfinity, fmt::format("C{}{}b", i + 1, j + 1).c_str()));
          back++;
          completion[back].setLinearCoef(C[i * I.n + j], 1.0);
          completion[back].setLinearCoef(C[i2 * I.n + j], -1.0);
          break;
        }
      }
    }
  }

  completion.add(IloRange(env, 0.0, 0.0, "Csum"));
  back++;
  completion[back].setLinearCoef(C[idx_c], 1.0);
  for (unsigned j = 0; j != I.n; ++j) {
    int i = I.lastOperation(j + 1);
    if (i-- == 0)
      continue;
    completion[back].setLinearCoef(C[i * I.n + j], -1.0);
  }

  model.add(completion);
}

void MPFSMO::addPrecedences() {
  IloRangeArray precedences(env);
  int back = -1;

  for (unsigned i = 0u; i != I.m; ++i)
    for (unsigned j1 = 0; j1 != I.n; ++j1)
      for (unsigned j2 = j1 + 1; j2 != I.n; ++j2) {
        precedences.add(IloRange(env, 1.0, 1.0, fmt::format("P{}{}", j1 + 1, j2 + 1).c_str()));
        back++;
        precedences[back].setLinearCoef(x[idx[{i, j1, j2}]], 1.0);
        precedences[back].setLinearCoef(x[idx[{i, j2, j1}]], 1.0);
      }

  model.add(precedences);
}

void MPFSMO::build(Time ub_cmax) {
  if (ub_cmax == infinite_time)
    ub_cmax = I.totalTime();
  addVars();
  addCompletion(ub_cmax);
  addPrecedences();
}

void MPFSMO::setSolution(IloCplex solver, const PSolution &S) {
  IloNumArray v(env, x.getSize());
  for (unsigned i = 0, ie = x.getSize(); i != ie; ++i)
    v[i] = 0.0;
  const auto ie = S.π.size();

  for (auto i1 = 1u; i1 != ie; ++i1)
    for (auto i2 = i1 + 1; i2 != ie; ++i2) {
      auto j1 = S.π[i1] - 1;
      auto j2 = S.π[i2] - 1;
      for (unsigned i = 0u; i != I.m; ++i)
        v[idx[{i, j1, j2}]] = 1.0;
    }
  solver.addMIPStart(x, v);
  v.end();
}

void MPFSMO::setSolution(IloCplex solver, const NPSolution &S) {
  IloNumArray v(env, x.getSize());
  for (unsigned i = 0, ie = x.getSize(); i != ie; ++i)
    v[i] = 0.0;

  for (unsigned i = 0u; i != I.m; ++i)
    for (auto i1 = 1u; i1 <= I.n; ++i1)
      for (auto i2 = i1 + 1; i2 <= I.n; ++i2) {
        auto j1 = S.π[i + 1][i1] - 1;
        auto j2 = S.π[i + 1][i2] - 1;
        v[idx[{i, j1, j2}]] = 1.0;
      }
  solver.addMIPStart(x, v);
  v.end();
}

NPSolution MPFSMO::getSolution(IloCplex solver) {
  IloNumArray value(env);
  solver.getValues(x, value);
  NPSolution S(I);
  for (unsigned i = 0u; i != I.m; ++i) {
    unsigned mach = i;
    vector<Job> π(I.n + 1);
    iota(π.begin(), π.end(), 0);
    sort(π.begin() + 1, π.end(), [&](Job j1, Job j2) { return value[idx[{mach, j1 - 1, j2 - 1}]] > 0.5; });
    for (unsigned j = 0u; j != I.n + 1; ++j)
      S.π[i + 1][j] = π[j];
  }
  S.of = solver.getObjValue();

  return S;
}

ModelStat MPFSMO::getStatistics() const { return mstat; }

pair<NPSolution, IloAlgorithm::Status> MPFSMO::solve(const PSolution &S) {
  IloCplex solver(model);
  if (m.exportModel)
    solver.exportModel("model.lp");
  mstat = ModelStat{unsigned(solver.getNrows()), unsigned(solver.getNcols()), unsigned(solver.getNNZs())};
  setupSolver(solver, m);
  if (S.isValid())
    setSolution(solver, S);
  solver.solve();
  IloAlgorithm::Status status = solver.getStatus();
  if (status == IloAlgorithm::Feasible || status == IloAlgorithm::Optimal)
    return {getSolution(solver), status};
  else
    return {NPSolution(I, S), status};
}
