/**
 * \file pfsmo.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 *
 * Permutation and non-permutation flow shop scheduling with missing operations.
 */
#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
namespace fs = std::filesystem;
using namespace std;

#include "heuristics.hpp"
#include "helpers.hpp"
#include "instance.hpp"
#include "logging.hpp"
#include "options.hpp"
#include "random.hpp"
#include "solution.hpp"

struct PFSOptions : public standardOptions {
  double timelimit;
  int iterlimit;
  double iterfactor;
  bool model;
  string psolution, npsolution;
  bool flowtime;
  bool npfs;
  string solution, wpsolution;

  PFSOptions() : timelimit(60) {}
};

PFSOptions opt;

void add_results(vector<Time> &v, const pair<Time, Time> &tt) {
  v.push_back(tt.first);
  v.push_back(tt.second);
}

int main(int argc, char *argv[]) {
  run::start.reset();

  
  IGAOptions iopt;
  std_description desc("Options", opt);
  
  desc.add_options()
    ("timelimit",    po::value<double>(&opt.timelimit)->default_value(0.0),    "Time limit for heuristics (seconds; default 5ms/op, negative for none).")
    ("iterlimit",    po::value<int>(&opt.iterlimit)->default_value(0.0),       "Iteration limit for heuristics (default 1.5×10⁵/n, negative for none).")
    ("iterfactor",   po::value<double>(&opt.iterfactor)->default_value(1.0),   "Multiplier for default iteration limit (which has been calibrated for about 5ms/op)")
    ("flowtime",     po::bool_switch(&opt.flowtime)->default_value(false),     "Make flowtime the primary objective.")
    ("npfs",         po::bool_switch(&opt.npfs)->default_value(false),         "Apply NPFS optimizations.")
    ;

  po::options_description iga("IGA options", get_terminal_width());
  iga.add_options()
    ("alpha",        po::value<double>(&iopt.alpha)->default_value(0.234375), "Alpha.")
    ("dc",           po::value<unsigned>(&iopt.dc)->default_value(8),         "D&C jobs.")
    ;

  po::options_description out("Output options", get_terminal_width());
  out.add_options()
    ("psolution",     po::value<string>(&opt.wpsolution)->default_value("/dev/null"), "File to write permutation solution to.")
    ("solution",      po::value<string>(&opt.solution)->default_value("/dev/null"), "File to write last solution to.")
    ;
  
  desc.add(iga).add(out);

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
  if (opt.timelimit == 0.0)
    opt.timelimit = double(I.n * I.m * 5) / 1000;
  if (opt.iterlimit == 0)
    opt.iterlimit = max(1.0, opt.iterfactor * double(150000) / I.n);
  iopt.dc = min(iopt.dc, (8 * I.n + 9) / 10);
  vprint(1, "Timelimit {:.12f}, iteration limit {}, dc {}.\n", opt.timelimit, opt.iterlimit, iopt.dc);

  EPSolution S(I);
  if (opt.flowtime)
    S.of_makespan = !opt.flowtime;

  vprint(1, "Optimizing for {}.\n", opt.flowtime ? "flowtime" : "makespan");
  
  vector<Result> results;
  vector<Time> npsset;

  S.totalTimeOrder();
  S.clear();
  S.insert_all();
  S.store_so();
  S.tfound = run::elapsed();
  results.push_back(S.getResultPO());
  vprint(1, "Construction {} ", results.back().to_string());
  results.push_back(S.getResultSO());
  vprint(1, "{}\n", results.back().to_string());

  add_results(npsset, S.evaluateNPSset(I));

  
  unsigned steps_shift = S.shift_ls();
  results.push_back(S.getResultPO());
  vprint(1, "Local search {} ", results.back().to_string());
  results.push_back(S.getResultSO());
  vprint(1, "{}\n", results.back().to_string());

  
  const double pavg = double(I.totalTime()) / (I.n * I.m);
  iopt.T = iopt.alpha * pavg / 10;
  iopt.timelimit = opt.timelimit - run::elapsed();
  iopt.iterlimit = opt.iterlimit;
  unsigned steps_iga = S.iga(iopt);
  results.push_back(S.getResultPO());
  vprint(1, "IGA {} ", results.back().to_string());
  results.push_back(S.getResultSO());
  vprint(1, "{} ", results.back().to_string());
  vprint(1, "\n");

  if (opt.wpsolution != "/dev/null") {
    ofstream sol(opt.wpsolution);
    if (sol.fail()) {
      fmt::print(cerr, "Failed to open {}\n", opt.wpsolution);
      return 1;
    }
    fmt::print(sol, "# {}\n", S.to_string());
    S.write(sol);
    sol.close();
  }

  add_results(npsset, S.evaluateNPSset(I));
  PSolution Sf{I, S.so.π};
  add_results(npsset, Sf.evaluateNPSset(I));
  
  unsigned steps_shift_np = 0;
  ENPSolution N(I, S);
  if (opt.npfs) {
    steps_shift_np = N.shift_ls();
    results.push_back(N.getResultPO());
    vprint(1, "Local search NP {}\n", results.back().to_string());
  } else
    results.push_back({0, 0, 0});

  fmt::print("PARAM {} {} {}\n", iopt.dc, opt.timelimit, opt.iterlimit);
  fmt::print("INFO {} ", iname);
  for (auto res : results)
    fmt::print("{} ", res.to_string());
  fmt::print("\n");
  fmt::print("STAT {} {} {}\n", steps_shift, steps_iga, run::elapsed());
  fmt::print("NSTAT {}\n", steps_shift_np);
  fmt::print("NPSSET {}\n", fmt::join(npsset.begin(), npsset.end(), " "));

  if (opt.solution != "/dev/null") {
    ofstream sol(opt.solution);
    if (sol.fail()) {
      fmt::print(cerr, "Failed to open {}\n", opt.solution);
      return 1;
    }
    N.write(sol);
    sol.close();
  }
}
