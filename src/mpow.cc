/*
! Description: matrix power
!
!        ME = MA^m
!
*/


#include <cmath>
#include "sci_spblas.h"

namespace mexp {

  sci::dmatrix<double>& mpow(const sci::dmatrix<double>& MA, sci::dmatrix<double>& ME, int m) {
    int info;
    int n = MA.nrow;

    sci::dmatrix<double> MX(MA);
    sci::dmatrix<double> MT(n,n);

    if (m < 0) {
      for (int i=1; i<=n; i++) {
        MT(i,i) = 1.0;
      }
      info = sci::dgesv(MX, MT, MT);
      return mpow(MT, ME, -m);
    }

    ME = 0.0;
    for (int i=1; i<=n; i++) {
      ME(i,i) = 1.0;
    }
    while (m != 0) {
      if (m % 2 == 1) {
        MT = ME;
        sci::dgemm(sci::mat::N, sci::mat::N, 1.0, MX, MT, 0.0, ME);
      }
      m /= 2;
      sci::dgemm(sci::mat::N, sci::mat::N, 1.0, MX, MX, 0.0, MT);
      MX = MT;
    }
    return ME;
  }

  sci::matrix<double>& mpow(const sci::matrix<double>& MA, sci::matrix<double>& ME, int m) {
    switch (PAIR(MA.type(),ME.type())) {
    case (PAIR(DENSE,DENSE)):
        return mpow(dynamic_cast<const sci::dmatrix<double>&>(MA),
          dynamic_cast<sci::dmatrix<double>&>(ME), m);
    default:
        throw;
    }
  }
}
