#include <math.h>
#include "pgamma.h"

static double PI = 3.14159265358979324; // pi
static double LOG_PI = 1.14472988584940017; // log(pi)

double d_normal(double x) {
  return 1.0/sqrt(2.0*PI)*exp(-x*x/2.0);
}

double p_normal(double x) {
  if (x >= 0.0)
    return 0.5 * (1 + p_gamma(0.5, 0.5 * x * x, LOG_PI / 2));
  else
    return 0.5 * q_gamma(0.5, 0.5 * x * x, LOG_PI / 2);
//    return 0.5 * (1.0 + erf(x/sqrt(2.0)));
}

double q_normal(double x) {
  if (x >= 0.0)
    return 0.5 * q_gamma(0.5, 0.5 * x * x, LOG_PI / 2);
  else
    return 0.5 * (1 + p_gamma(0.5, 0.5 * x * x, LOG_PI / 2));
//    return 0.5 * (1.0 + erf(-x/sqrt(2.0)));
}



