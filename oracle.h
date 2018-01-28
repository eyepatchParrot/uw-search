#ifndef ORACLE_H
#define ORACLE_H

#include "util.h"
#include <vector>
#include <cstdlib>

class Oracle {
  using Index = int;
  const std::vector<Key>& a;
  std::vector<Index> i;
  int j;

  public:
  // note that to ensure fast code, offset is always to the left, which means that offset searches wont have the full offset
  Oracle(const std::vector<Key>& _a, const std::vector<Index>& _i, const Index offset = 0, const bool rnd = false) : a(_a), j(0) {
    // easier to make a single long pipeline instead of a branch
    // There's blocking on the output to be produced, instead of being able to go through the steady state throughput
    // You may end up paying the latency cost multiple times instead of just once with left-deep
    // single relation plan because only one table in FROM clause
    for (Index ti : _i) {
      Index tti = ti;
      if (rnd && ((rand() % 2) == 0)) {
        tti += offset;
      } else {
        tti -= offset;
      }
      unsigned r = a.size() - 1;
      tti = tti < 0 ? 0 : tti;
      tti = tti > r ? r : tti;
      i.push_back(tti);
    }
  }

  Key operator()(const Key x) {
    Index k = i[j++];
    assert(a[k] == x);
    j = j >= i.size() ? 0 : j;
    return a[k] == x ? a[k] : -1;
  }
};

#endif
