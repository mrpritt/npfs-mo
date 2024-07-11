/**
 * \file holes.hpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#pragma once

#include <iostream>
#include <set>
#include <utility>

#include "instance.hpp"
#include "logging.hpp"

struct Interval {
  Time p, f;
  friend bool operator<(const Interval &i, const Interval &j) { return i.p < j.p || (!(i.p > j.p) && i.f < j.f); }
  friend std::ostream &operator<<(std::ostream &o, const Interval &i) { return o << "[" << i.start() << "," << i.f << "]"; }
  Time duration(Time C = 0) const {
    auto s = std::max(start(), C);
    return s < f ? f - s : 0;
  }
  Time start() const { return f - p; }
};

template <> struct fmt::formatter<Interval> : ostream_formatter {};

struct hlist {
  using element = typename std::set<Interval>::iterator;

  std::set<Interval> h;

  hlist() { h.insert({infinite_time, infinite_time}); }

  unsigned size() const { return h.size(); }
  element earliest(Time p, Time C);
  element smallest(Time p, Time C);
  bool reduce(element, Time d);
  void cut(element, Time C, Time p);

  std::string to_string() const;
};
