/*
! Description: matrix exp with uniformization
!
!        y = exp(Q*t)
!
!        Q is uniformized to P and qv
!        t is involved in the Poisson probability vector.
*/

#include <cmath>

#include "poisson.hh"
#include "mexp.hh"

namespace mexp {

  sci::matrix<double>& mexp_unif(const sci::matrix<double>& P, double qv,
    const sci::range& r, const sci::vector<double>& poivec, double weight,
    const sci::dmatrix<double>& x, sci::dmatrix<double>& y,
    double atol) {

    sci::vector<double>::const_range poi = poivec.alias(r);
    int left = r.begin;
    int right = r.end;

    int n = P.nrow;
    sci::dmatrix<double> xi(x);
    sci::dmatrix<double> tmp(n,n);
    y = 0.0;
    sci::daxpy(poi(left), xi, y);
    for (int k=left+1; k<=right; k++) {
      tmp = xi;
      sci::dgemm(sci::mat::N, sci::mat::N, 1.0, P, tmp, 0.0, xi);
      sci::daxpy(poi(k), xi, y);
      if (damax(xi) < atol) break;
    }
    sci::dscal(1.0/weight, y);
    return y;
  }

  sci::vector<double>& mexp_unifvec(sci::mat::trans tr,
    const sci::matrix<double>& P, double qv,
    const sci::range& r, const sci::vector<double>& poivec, double weight,
    const sci::vector<double>& x, sci::vector<double>& y,
    double atol) {

    sci::vector<double>::const_range poi = poivec.alias(r);
    int left = r.begin;
    int right = r.end;

    int n = P.nrow;
    sci::vector<double> xi(x);
    sci::vector<double> tmp(n);
    y = 0.0;
    sci::daxpy(poi(left), xi, y);
    for (int k=left+1; k<=right; k++) {
      tmp = xi;
      sci::dgemv(tr, 1.0, P, tmp, 0.0, xi);
      sci::daxpy(poi(k), xi, y);
      if (sci::damax(xi) < atol) break;
    }
    sci::dscal(1.0/weight, y);
    return y;
  }

}
