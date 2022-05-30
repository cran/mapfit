/*
! Description: integral operation for matrix exp form with threads;
!
!        y = exp(trans(Q)*t) * x
!
!        Q is uniformized to P and qv
!        t is involved in the Poisson probability vector.
*/

#include <cmath>

#include "poisson.h"
#include "mexp.h"

namespace mexp {

  sci::vector<double>& mexpi_unifvec_thread(sci::mat::trans trans,
    const sci::array< sci::matrix<double>* >& P, double qv,
    const sci::range& r, const sci::vector<double>& poivec, double weight,
    const sci::vector<double>& x, sci::vector<double>& y, sci::vector<double>& cy,
    double atol) {

    if (r.begin != 0) {
      throw;
    }

    int n = (*(P[0])).nrow;
    int m = P.size;
    int th = pow2i(m-1);
    sci::range xrow(1, n);
    sci::dmatrix<double> xi(n, th);
    sci::dmatrix<double> tmp(n, th);

    sci::vector<double>::const_range poi = poivec.alias(r);
    int left = r.begin;
    int right = r.end;

    sci::vector<double> cpoivec(right-left+1);
    sci::vector<double>::range cpoi = cpoivec.alias(sci::range(r));
    cpoi(right) = 0.0;
    for (int k=right-1; k>=left; k--) {
      cpoi(k) = cpoi(k+1) + poi(k+1);
    }

    th = 1;
    y = 0.0;
    sci::dscal(qv*weight, cy);

    xi(xrow,1) = x;
    int k = 0;
    int i = 0;
    for (i=0; i<m-1; i++) {
      if (k + th >= right) {
        sci::dgemv(sci::mat::N, 1.0, xi(xrow,sci::range(1,k+1)), poi(sci::range(0,k)), 1.0, y);
        sci::dgemv(sci::mat::N, 1.0, xi(xrow,sci::range(1,k+1)), cpoi(sci::range(0,k)), 1.0, cy);
        goto FINAL;
      }
      k += th;
      sci::dmatrix<double> c = xi(xrow,sci::range(th+1,k+1));
      sci::dgemm(trans, sci::mat::N, 1.0, *(P[i]), xi(xrow,sci::range(1,th)), 0.0, c);
      th *= 2;
    }
    sci::dgemv(sci::mat::N, 1.0, xi, poi(sci::range(0,k)), 1.0, y);
    sci::dgemv(sci::mat::N, 1.0, xi, cpoi(sci::range(0,k)), 1.0, cy);

    i = m - 1;
    do {
      if (k + th >= right) goto FINAL;
      tmp = xi;
      sci::dgemm(trans, sci::mat::N, 1.0, *(P[i]), tmp, 0.0, xi);
      if (sci::damax(xi) < atol) goto FINAL2;
      sci::dgemv(sci::mat::N, 1.0, xi, poi(sci::range(k+1,k+th)), 1.0, y);
      sci::dgemv(sci::mat::N, 1.0, xi, cpoi(sci::range(k+1,k+th)), 1.0, cy);
      k += th;
    } while (1);

    FINAL:

    // final block
    {
      int rnum = right - k;
      tmp = xi;
      sci::dmatrix<double> c = xi(xrow,sci::range(1,rnum));
      sci::dgemm(trans, sci::mat::N, 1.0, *(P[i]), tmp(xrow,sci::range(1,rnum)), 0.0, c);
      sci::dgemv(sci::mat::N, 1.0, xi(xrow,sci::range(1,rnum)), poi(sci::range(k+1,k+rnum)), 1.0, y);
      sci::dgemv(sci::mat::N, 1.0, xi(xrow,sci::range(1,rnum)), cpoi(sci::range(k+1,k+rnum)), 1.0, y);
    }

    FINAL2:

    sci::dscal(1.0/weight, y);
    sci::dscal(1.0/qv/weight, cy);

    return y;
  }

}

