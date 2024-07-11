/**
 * \file holes.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 */
#include "holes.hpp"

#include "logging.hpp"

using namespace std;

hlist::element hlist::earliest(Time p, Time C) {
  auto e = h.lower_bound({p, C});
  auto b = --h.end();

  while (e != h.end()) {
    if (p <= e->duration(C) && e->start() < b->start())
      b = e;
    e++;
  }
  return b;
}

hlist::element hlist::smallest(Time p, Time C) {
  auto e = h.lower_bound({p, 0});
  while (e->duration(C) < p)
    e++;
  return e;
}

bool hlist::reduce(element e, Time d) {
  if (e->p <= d) {
    h.erase(e);
    return true;
  } else {
    auto v = *e;
    h.erase(e);
    v.p -= d;
    h.insert(v);
    return false;
  }
}

void hlist::cut(element e, Time C, Time p) {
  auto s = e->start();
  assert(s <= C);
  reduce(e, C + p - s);
  if (s < C)
    h.insert({C - s, C});
}

string hlist::to_string() const {
  string s = "";
  for (auto f : h)
    s += fmt::format("{} ", f);
  return s;
}
