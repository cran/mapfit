/*
Poisson prob
*/

#include "poisson.h"
#include "sci_spblas.hh"

namespace pois {

  int rightbound(double lambda, double eps) {
    return poisson_rightbound(lambda, eps);
  }

  double pmf(double lambda, const sci::range& range, sci::vector<double>& x) {
    double weight;
    poisson_prob(lambda, range.begin, range.end, x.ptr, &weight);
    return weight;
  }

}



