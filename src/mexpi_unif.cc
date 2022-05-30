/*
  ! Description: integral operation for matrix exp form;
  !
  !                    |t
  !        cME = cME + | exp(Q*s) ds
  !                    |0
  !
  !        ME = exp(Q*t)
  !
  !        Q is uniformized to P and qv
  !        t is involved in the Poisson probability vector.
  !        return value is ME
 */

#include <cmath>

#include "poisson.h"
#include "mexp.h"

namespace mexp {

  sci::dmatrix<double>& mexpi_unif(const sci::matrix<double>& P, double qv,
    const sci::range& r, const sci::vector<double>& poivec, double weight,
    const sci::dmatrix<double>& x, sci::dmatrix<double>& y, sci::dmatrix<double>& cy) {

    sci::vector<double>::const_range poi = poivec.alias(r);
    int left = r.begin;
    int right = r.end;
    sci::vector<double> cpoivec(right-left+1);
    sci::vector<double>::range cpoi = cpoivec.alias(r);

    cpoi(right) = 0.0;
    for (int k=right-1; k>=left; k--) {
      cpoi(k) = cpoi(k+1) + poi(k+1);
    }

    int n = P.nrow;
    sci::dmatrix<double> xi(x);
    sci::dmatrix<double> tmp(n,n);
    y = 0.0;
    sci::dscal(qv*weight, cy); // original cy is scaled

    sci::daxpy(poi(left), xi, y);
    sci::daxpy(cpoi(left), xi, cy);
    for (int k=left+1; k<=right; k++) {
      tmp = xi;
      sci::dgemm(sci::mat::N, sci::mat::N, 1.0, P, tmp, 0.0, xi);
      sci::daxpy(poi(k), xi, y);
      sci::daxpy(cpoi(k), xi, cy);
    }
    sci::dscal(1.0/weight, y);
    sci::dscal(1.0/qv/weight, cy);
    return y;
  }

  sci::vector<double>& mexpi_unifvec(sci::mat::trans tr,
    const sci::matrix<double>& P, double qv,
    const sci::range& r, const sci::vector<double>& poivec, double weight,
    const sci::vector<double>& x, sci::vector<double>& y, sci::vector<double>& cy,
    double atol) {

    sci::vector<double>::const_range poi = poivec.alias(r);
    int left = r.begin;
    int right = r.end;
    sci::vector<double> cpoivec(right-left+1);
    sci::vector<double>::range cpoi = cpoivec.alias(r);

    cpoi(right) = 0.0;
    for (int k=right-1; k>=left; k--) {
      cpoi(k) = cpoi(k+1) + poi(k+1);
    }

    int n = P.nrow;
    sci::vector<double> xi(x);
    sci::vector<double> tmp(n);
    y = 0.0;
    sci::dscal(qv*weight, cy);

    sci::daxpy(poi(left), xi, y);
    sci::daxpy(cpoi(left), xi, cy);
    for (int k=left+1; k<=right; k++) {
      tmp = xi;
      sci::dgemv(tr, 1.0, P, tmp, 0.0, xi);
      sci::daxpy(poi(k), xi, y);
      sci::daxpy(cpoi(k), xi, cy);
      if (sci::damax(xi) < atol) break;
    }
    sci::dscal(1.0/weight, y);
    sci::dscal(1.0/qv/weight, cy);
    return y;
  }
}

