#ifndef ADAPTIVESEARCH_H
#define ADAPTIVESEARCH_H

#include <cinttypes>
#include <assert.h>

#include "lin.h"
#include "util.h"
#include "padded_vector.h"

#include <limits>

template <int record_bytes = 8, typename Index = unsigned long>
class adaptivesearch {
  using Vector = PaddedVector<record_bytes>;
  using Linear = LinearUnroll<Vector>;

  const Vector& A;

 public:
  adaptivesearch(const Vector &_a) : A(_a) {}

  __attribute__((always_inline)) Key operator()(const Key x) {
    Index bot, top, next, med, position;
    bool finished = false;
    bot = 0;
    top = A.size() - 1;
    position = -1;
    while((position == -1) && (bot < top)) {
      med = (top - bot) / 2;

      next = (int) ((x - A[bot]) / (A[top] - A[bot]) * (top - bot) + bot);

      if ((x < A[next]) && ((next - bot) > med)) {
        top = next - 1;
        next = (bot + top) / 2;
      } else if ((x > A[next]) && ((top - next) > med)) {
        bot = next + 1;
        next = (bot + top) / 2;
      }

      if (x > A[next]) {
        bot = next + 1;
      } else if (x < A[next]) {
        top = next - 1;
      } else {
        position = next;
      }
    }

    if(A[top] == x) {
      position = top;
    } else {
      position = -top - 1;
    }

    return A[position];
  }

};


#endif //ADAPTIVESEARCH_H
