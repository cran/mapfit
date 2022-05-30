/*
  nmath header
*/

#ifndef PGAMMA_H
#define PGAMMA_H

#ifdef __cplusplus
extern "C" {
#endif

  double p_gamma(double a, double x, double loggamma_a);
  double q_gamma(double a, double x, double loggamma_a);

#ifdef __cplusplus
}
#endif

#endif
