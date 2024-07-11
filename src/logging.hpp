/**
 * \file logging.hpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#pragma once

#include <iostream>
#include <vector>

#include <boost/any.hpp>

#include "timer.hpp"

#define FMT_HEADER_ONLY
#include "fmt/format.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h"

struct options_counter {
  int count = 0;
};
inline void validate(boost::any &v, std::vector<std::string> const &xs, options_counter *, long) {
  if (v.empty())
    v = options_counter{1};
  else
    ++boost::any_cast<options_counter &>(v).count;
}

namespace run {
extern timer start;
extern options_counter verbosec;
inline double elapsed() { return start.elapsed(); }
} // namespace run

inline bool verbose(int level = 1) { return run::verbosec.count >= level; }

template <typename S, typename... Args> void vprint(int level, const S &format_str, Args &&...args) {
  if (run::verbosec.count >= level)
    fmt::print(format_str, std::forward<Args>(args)...);
}
