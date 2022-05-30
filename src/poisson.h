/*
  Poisson
*/

#pragma once

#include <cfloat>
#include "sci_spblas.h"

namespace pois {

  int rightbound(double lambda, double eps = DBL_EPSILON);
  double pmf(double lambda, const sci::range& range, sci::vector<double>& x);

}



