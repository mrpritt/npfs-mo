/**
 * \file helpers.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#include "helpers.hpp"

#include <regex>
using namespace std;

#define FMT_HEADER_ONLY
#include "fmt/format.h"

string canonical_name(string fname) {
  const regex hn_fname("([\\d\\.]+)_(\\d+)_(\\d+)_(\\d+)\\.txt");
  std::smatch match;
  if (regex_match(fname, match, hn_fname))
    return fmt::format("{} {} {} {}", match[1].str(), match[2].str(), match[3].str(), match[4].str());
  else
    return fname;
}
