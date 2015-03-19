/*
! Description: matrix power
!
!        ME = MA^m
!
*/


#include <cmath>
#include "sci_spblas.hh"

namespace mexp {

  sci::dmatrix<double>& mexp_pade(const sci::matrix<double>& MA, sci::dmatrix<double>& ME, double eps) {

    double norma, tolerr, c;
    int i, j, k, q, info;

    int n = MA.nrow;
    sci::dmatrix<double> MN(n,n);
    sci::dmatrix<double> MD(n,n);
    sci::dmatrix<double> MX(n,n);
    sci::dmatrix<double> MT(n,n);

//! scaled A
    dcopy(MA, ME);
    norma = damax(ME);
    j = log(norma)/log(2.0);
    j = (0 < 1+j) ? 1+j : 0;
    sci::dscal(1.0/std::ldexp(1.0, j), ME);

//! find q
    q = 1;
    tolerr = 1.0 / 6.0;
    while (tolerr > eps/norma) {
      tolerr /= 16.0 * (3.0 + 4.0 * q * (2.0 + q));
      q++;
    }

//! initialize
    c = 1.0;
    for (int i=1; i<=n; i++) {
      MD(i,i) = 1.0;
      MN(i,i) = 1.0;
      MX(i,i) = 1.0;
    }

//! make
    for (k=1; k<=q; k++) {
      c *= (q - k + 1.0)/ ((2.0*q - k + 1.0) * k);
      sci::dgemm(sci::mat::N, sci::mat::N, 1.0, ME, MX, 0.0, MT);
      MX = MT;
      sci::daxpy(c, MX, MN);
      if (k % 2 == 0) {
        sci::daxpy(c, MX, MD);
      } else {
        sci::daxpy(-c, MX, MD);
      }
    }
    ME = MN;
    info = sci::dgesv(MD, ME, ME);
    for (k=1; k<=j; k++) {
      sci::dgemm(sci::mat::N, sci::mat::N, 1.0, ME, ME, 0.0, MT);
      ME = MT;
    }
    return ME;
  }

}
