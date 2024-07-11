/**
 * \file models.hpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 *
 * Implementation of mathematical models.
 *
 */
#pragma once

#include <iostream>
#include <map>
#include <string>
#include <tuple>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "instance.hpp"
#include "solution.hpp"

void setupSolver(IloCplex solver);
string to_string(IloAlgorithm::Status status);

struct ModelOptions {
  bool exportModel;
  unsigned timelimit;

  ModelOptions() : timelimit(30) {}
};

struct ModelStat {
  unsigned rows, cols, nnz;
  std::string to_string() { return fmt::format("{} {} {}", rows, cols, nnz); }
};

struct MPFSMO {
  const Instance &I;
  IloEnv env;
  IloModel model;
  IloObjective obj;

  ModelOptions m;
  ModelStat mstat;

  map<tuple<unsigned, unsigned, unsigned>, unsigned> idx;
  unsigned idx_c;

  IloNumVarArray C, x;

  MPFSMO(const Instance &I, const ModelOptions &m) : I(I), model(env), obj(IloAdd(model, IloMinimize(env))), m(m), C(env), x(env) {}
  void addVars();
  void addCompletion(Time = infinite_time);
  void addPrecedences();

  void build(Time);

  void setSolution(IloCplex, const PSolution &);
  void setSolution(IloCplex, const NPSolution &);
  NPSolution getSolution(IloCplex);
  ModelStat getStatistics() const;

  std::pair<NPSolution, IloAlgorithm::Status> solve(const PSolution &);
};
