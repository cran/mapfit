/*
  ! Description: stationary vector via GTH algorithm
  !        x = x * exp(Q),   x * Q = 0
 */

#include <cmath>
#include "mexp.hh"

namespace mexp {

  sci::vector<double>& ctmc_gth(const sci::dmatrix<double>& Q, sci::vector<double>& x) {

    int n = x.size;

    double tmp;
    sci::dmatrix<double> A(Q);

    for (int l=n; l>=2; l--) {
      sci::range rr(1,l-1);

      tmp = sci::dasum(A(l,rr));
      for (int j=1; j<=l-1; j++) {
        for (int i=1; i<=l-1; i++) {
          if (i != j) {
            A(i,j) += A(l,j) * A(i,l) / tmp;
          }
        }
      }
      sci::vector<double> tmpv = A(rr,l);
      sci::dscal(1.0/tmp, tmpv);
      A(l,rr) = 0.0;
      A(l,l) = -1.0;
    }

    x(1) = 1.0;
    for (int l=2; l<=n; l++) {
      sci::range rr(1,l-1);
      x(l) = sci::ddot(x(rr), A(rr,l));
    }
    sci::dscal(1.0/dsum(x), x);
    return x;
  }

}
