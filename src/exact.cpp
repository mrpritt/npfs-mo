/**
 * \file exact.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <filesystem>
namespace fs = std::filesystem;

using namespace std;

#include "helpers.hpp"
#include "heuristics.hpp"
#include "logging.hpp"
#include "models.hpp"
#include "options.hpp"
#include "random.hpp"

struct SolverOptions : public standardOptions, public ModelOptions {
  string solution;
};

int main(int argc, char *argv[]) {
  SolverOptions opt;
  std_description desc("Options", opt);
  // clang-format off
  desc.add_options()
    ("timelimit",    po::value<unsigned>(&opt.timelimit)->default_value(0),   "Time limit (seconds).")
    ;

  po::options_description out("Output options", get_terminal_width());
  out.add_options()
    ("exportModel",  po::bool_switch(&opt.exportModel)->default_value(false),      "Export model.")
    ("solution",     po::value<string>(&opt.solution)->default_value("/dev/null"), "File to write solution to.")
    ;
  // clang-format on

  desc.add(out);

  po::positional_options_description pod;
  pod.add("instance", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("instance")) {
    if (!vm.count("instance"))
      fmt::print("No instance given.\n\n");
    cout << desc << endl;
    return 0;
  }
  opt.seed = setupRandom(opt.seed);

  string iname = canonical_name(fs::path(opt.instance).filename());
  ifstream ins(opt.instance);
  if (ins.fail()) {
    fmt::print(cerr, "Failed to open {}\n", opt.instance);
    return 1;
  }
  Instance I(ins);
  ins.close();

  vprint(1, "Instance with {} jobs and {} machines, missing operations rate {}.\n", I.n, I.m, I.r);
  IGAOptions iopt;
  if (opt.timelimit == 0.0)
    opt.timelimit = min(max(double(I.n * I.m), 30.0), 600.0);
  const double iterfactor = 0.1;
  iopt.iterlimit = max(1.0, iterfactor * double(150000) / I.n);
  iopt.dc = min(iopt.dc, (8 * I.n + 9) / 10);
  vprint(1, "Timelimit {}, iteration limit for IGA {}, dc {}.\n", opt.timelimit, iopt.iterlimit, iopt.dc);

  EPSolution S(I);
  S.of_makespan = false;

  vector<Result> results;

  S.totalTimeOrder();
  S.clear();
  S.insert_all();
  S.store_so();
  S.tfound = run::elapsed();

  const double pavg = double(I.totalTime()) / (I.n * I.m);
  iopt.T = iopt.alpha * pavg / 10;
  iopt.timelimit = infinite_time;

  S.iga(iopt);
  results.push_back(S.getResultPO());
  vprint(1, "IGA {}\n", results.back().to_string());

  ModelStat mstat{0, 0, 0};
  MPFSMO m(I, opt);
  m.build(S.getMakespan(I));

  IloAlgorithm::Status status = IloAlgorithm::Unknown;
  NPSolution Sm(I);
  tie(Sm,status) = m.solve(S);
  results.push_back(Result{Sm.getMakespan(I), Sm.of, run::elapsed()});
  vprint(1, "Model results {}\n", results.back().to_string());
  mstat = m.getStatistics();

  fmt::print("INFO {} ", iname);
  for (auto res : results)
    fmt::print("{} ", res.to_string());
  fmt::print("{} ", to_string(status));
  fmt::print("\n");
  fmt::print("STAT {} {}\n", iname, mstat.to_string());

  if (opt.solution != "/dev/null") {
    ofstream sol(opt.solution);
    if (sol.fail()) {
      fmt::print(cerr, "Failed to open {}\n", opt.solution);
      return 1;
    }
    Sm.write(sol);
    sol.close();
  }
}
