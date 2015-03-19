/*
! Description: make DTMC probability transition matrix from
!              an infinitesimal CTMC generator
!                 P = I + Q/qv
!                 qv = unif_factor * max|Q|
! Parameters
!      n: size of CTMC (DTMC) kernel
!      Q: CTMC kernel (sparse matrix)
!      P: DTMC transition matrix (out)
!      qv: maximum of CTMC diagonals (out)
!      ufact: uniformization factor
*/

#include <cmath>

#include "sci_spblas.hh"

namespace mexp {

// for dense matrix
  double unif(const sci::dmatrix<double>& Q, sci::dmatrix<double>& P, double ufact) {
    int n = Q.nrow;
    P = Q;
    double qv = 0.0;
    for (int i=1; i<=n; i++) {
      if (std::abs(P(i,i)) > qv) {
        qv = std::abs(P(i,i));
      }
    }
    qv *= ufact;

    sci::dscal(1.0/qv, P);
    for (int i=1; i<=n; i++) {
      P(i,i) += 1.0;
    }
    return qv;
  }

// for sparse matrix
  double unif(const sci::csrmatrix<double>& Q, sci::csrmatrix<double>& P, double ufact) {
    int n = Q.nrow;
    P = Q;

    sci::vector<int> idiag(n);
    sci::array<int>::range rowptr = P.get_rowptr();
    sci::array<int>::range colind = P.get_colind();
    sci::array<double>::range value = P.get_value();

    double qv = 0.0;
    for (int i=1; i<=n; i++) {
      for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
        if (colind[z] == i) {
          idiag(i) = z;
          if (std::abs(value[z]) > qv) {
            qv = std::abs(value[z]);
          }
          break;
        }
      }
    }
    qv *= ufact;

    sci::dscal(1.0/qv, P);
    for (int i=1; i<=n; i++) {
      value[idiag(i)] += 1.0;
    }
    return qv;
  }

  double unif(const sci::cscmatrix<double>& Q, sci::cscmatrix<double>& P, double ufact) {
    int n = Q.nrow;
    P = Q;

    sci::vector<int> idiag(n);
    sci::array<int>::range colptr = P.get_colptr();
    sci::array<int>::range rowind = P.get_rowind();
    sci::array<double>::range value = P.get_value();

    double qv = 0.0;
    for (int i=1; i<=n; i++) {
      for (int z=colptr[i]; z<colptr[i+1]; z++) {
        if (rowind[z] == i) {
          idiag(i) = z;
          if (std::abs(value[z]) > qv) {
            qv = std::abs(value[z]);
          }
          break;
        }
      }
    }
    qv *= ufact;

    sci::dscal(1.0/qv, P);
    for (int i=1; i<=n; i++) {
      value[idiag(i)] += 1.0;
    }
    return qv;
  }

  double unif(const sci::coomatrix<double>& Q, sci::coomatrix<double>& P, double ufact) {
    int n = Q.nrow;
    P = Q;

    sci::vector<int> idiag(n);
    sci::array<int>::range rowind = P.get_rowind();
    sci::array<int>::range colind = P.get_colind();
    sci::array<double>::range value = P.get_value();

    double qv = 0.0;
    for (int z=1; z<=Q.nnz; z++) {
      if (rowind[z] == colind[z]) {
        idiag(rowind[z]) = z;
        if (std::abs(value[z]) > qv) {
          qv = std::abs(value[z]);
        }
      }
    }
    qv *= ufact;

    sci::dscal(1.0/qv, P);
    for (int i=1; i<=n; i++) {
      value[idiag(i)] += 1.0;
    }
    return qv;
  }

  double unif(const sci::matrix<double>& Q, sci::matrix<double>& P, double ufact) {
    switch (PAIR(Q.type(),P.type())) {
    case (PAIR(DENSE,DENSE)):
        return unif(dynamic_cast<const sci::dmatrix<double>&>(Q),
          dynamic_cast<sci::dmatrix<double>&>(P), ufact);
    case (PAIR(CSR,CSR)):
        return unif(dynamic_cast<const sci::csrmatrix<double>&>(Q),
          dynamic_cast<sci::csrmatrix<double>&>(P), ufact);
    case (PAIR(CSC,CSC)):
        return unif(dynamic_cast<const sci::cscmatrix<double>&>(Q),
          dynamic_cast<sci::cscmatrix<double>&>(P), ufact);
    case (PAIR(COO,COO)):
        return unif(dynamic_cast<const sci::coomatrix<double>&>(Q),
          dynamic_cast<sci::coomatrix<double>&>(P), ufact);
    default:
        throw;
    }
  }

}


