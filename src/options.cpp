/**
 * \file options.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#include "options.hpp"
#include <sys/ioctl.h>

using namespace std;

#include "logging.hpp"

namespace po = boost::program_options;

void add_std_options(po::options_description &desc, standardOptions &opt) {
  desc.add_options()
      // clang-format off
    ("help",							      "Show help.")
    ("version",                                                       "Show version.")
    ("verbose,v",   po::value(&run::verbosec)->zero_tokens(),         "Verbosity. If present, extra output is shown. If -v is repeated, more output is given.")
    ("instance",    po::value<string>(&opt.instance),	              "Instance.")
    ("seed",        po::value<unsigned>(&opt.seed)->default_value(1), "Random seed (0 for random value).")
    ;
  // clang-format on
}

unsigned get_terminal_width() {
  struct winsize wsize;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &wsize);
  return wsize.ws_col;
}

std_description::std_description(const string &caption, standardOptions &opt) : Base(caption, get_terminal_width()) { add_std_options(*this, opt); }
