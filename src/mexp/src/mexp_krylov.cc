/*
  ! Description: matrix exp form with Krylov subspace;
  !
  !        y = exp(trans(Q)*s) * x
  !
  !        return value is y
 */

#include <cmath>

#include "mexp.hh"

namespace mexp {

  sci::vector<double>& mexp_krylov_pade(sci::mat::trans trans,
    const sci::spmatrix<double>& Q, double t,
    const sci::vector<double>& x, sci::vector<double>& y,
    int m, double& err, int arnoldi_ite, double arnoldi_tol,
    double pade_eps) {

    int n = Q.nrow;
    double beta, err1, err2;
    bool alloc = false;

    sci::dmatrix<double> H(m+2,m+2);
    sci::dmatrix<double> V(n,m+1);
    sci::dmatrix<double> tmpM(m+2,m+2);
    sci::vector<double> tmpx(n);

    sci::dmatrix<double> subH = H(sci::range(1,m+1),sci::range(1,m+1));
    darnoldi(trans, m+1, Q, x, subH, V, beta, arnoldi_ite, arnoldi_tol);
    H(sci::range(1,m+1),m+1) = 0.0;
    H(m+2,m+1) = 1.0;

    tmpM = 0.0;
    sci::daxpy(t, H, tmpM);
    mexp::mexp_pade(tmpM, H, pade_eps);
    sci::dgemv(sci::mat::N, beta, V(sci::range(1,n),sci::range(1,m)),
      H(sci::range(1,m),1), 0.0, y);

    err1 = beta * std::abs(H(m+1,1));
    sci::dgemv(sci::mat::N, 1.0, Q, V(sci::range(1,n),m+1), 0.0, tmpx);
    err2 = beta * std::abs(H(m+2,1)) * sci::dnrm2(tmpx);
    if (err1 > 10.0*err2) {
      err = err2;
    } else if (err1 > err2) {
      err = err1*err2 / (err1 - err2);
    } else {
      err = err1;
    }

    return y;
  }

}

