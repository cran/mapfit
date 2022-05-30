/*
  nmath header
*/

#ifndef NGAMMA_H
#define NGAMMA_H

#ifdef __cplusplus
extern "C" {
#endif

  double mylgamma(double x);
  double mytgamma(double x);
  double psi(double x);
  double polygamma(int n, double x);
  double tfact(int s);
  double lfact(int s);

#ifdef __cplusplus
}
#endif

#endif
