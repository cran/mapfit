/*
  gauss inte
*/

#ifndef GAUSS_INTE_H
#define GAUSS_INTE_H

#ifdef __cplusplus
extern "C" {
#endif

  void gauss_inte_w(int n, double *x, double *w, double eps);
  double gauss_inte_fx(int n, double *x, double a, double b, double *fx);
  double gauss_inte_fv(int n, double *w, double c, double *fv);

#ifdef __cplusplus
}
#endif

#endif
