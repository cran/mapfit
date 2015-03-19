/*
   mixed Erlang distribution
*/

#include <cmath>
#include <memory>

#include "gamma.h"
#include "pgamma.h"

#define ERL_MAX 10

namespace mapfit {

  double erlang_pdf(int a, double b, double t) {
    return b * pow(b*t, a-1) * exp(-b * t) / tfact(a-1);
  }

  double erlang_lpdf(int a, double b, double t) {
    return a * log(b) + (a-1) * log(t) - b * t - lfact(a-1);
  }

  double erlang_cdf(int a, double b, double t) {
    if (a <= ERL_MAX) {
      double w = 1.0;
      for (int i=1; i<=a-1; i++) {
        w += pow(b*t, i) / tfact(i);
      }
      return 1.0 - exp(-b*t) * w;
    } else {
      return p_gamma(a+1.0, b*t, mylgamma(a+1.0));
    }
  }

  double erlang_ccdf(int a, double b, double t) {
    if (a <= ERL_MAX) {
      double w = 1.0;
      for (int i=1; i<=a-1; i++) {
        w += pow(b*t, i) / tfact(i);
      }
      return exp(-b*t) * w;
    } else {
      return q_gamma(a+1.0, b*t, mylgamma(a+1.0));
    }
  }

}
