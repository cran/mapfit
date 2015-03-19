/*
  ! Description: integral operation for matrix exp form with Krylov subspace;
  !
  !                  |t
  !        cy = cy + | exp(trans(Q)*s) ds * x
  !                  |0
  !
  !        y = exp(trans(Q)*t) * x
  !
  !        return value is y
 */

#include <cmath>
#include <memory>

#include "mexp.hh"

namespace mexp {

  // prototype
  sci::spmatrix<double>* addelem(sci::mat::trans trans, const sci::spmatrix<double>& Q);
  sci::spmatrix<double>* addelem(sci::mat::trans trans, const sci::csrmatrix<double>& Q);
  sci::spmatrix<double>* addelem(sci::mat::trans trans, const sci::cscmatrix<double>& Q);
  sci::spmatrix<double>* addelem(sci::mat::trans trans, const sci::coomatrix<double>& Q);

  sci::vector<double>& mexpi_krylov_pade(sci::mat::trans trans,
    const sci::spmatrix<double>& Q, double t,
    const sci::vector<double>& x, sci::vector<double>& y, sci::vector<double>& cy,
    int m, double& err, int arnoldi_ite, double arnoldi_tol,
    double pade_eps) {

    int n = Q.nrow;
    double beta, err1, err2;
    sci::dmatrix<double> H(m+2,m+2);
    sci::dmatrix<double> V(n+n,m+1);
    sci::dmatrix<double> tmpM(m+2,m+2);
    sci::vector<double> tmpx(n+n);

    sci::vector<double> xi(n+n);
    xi(sci::range(1,n)) = x;

    // make Q2
//    std::unique_ptr< sci::spmatrix<double> > Q2(addelem(trans, Q));
    sci::spmatrix<double>* Q2 = addelem(trans, Q);
    sci::dmatrix<double> subH = H(sci::range(1,m+1),sci::range(1,m+1));
    sci::darnoldi(trans, m+1, *Q2, xi, subH, V, beta, arnoldi_ite, arnoldi_tol);
    H(sci::range(1,m+1),m+1) = 0.0;
    H(m+2,m+1) = 1.0;

    tmpM = 0.0;
    sci::daxpy(t, H, tmpM);
    mexp::mexp_pade(tmpM, H, pade_eps);
    sci::dgemv(sci::mat::N, beta, V(sci::range(1,n+n),sci::range(1,m)),
      H(sci::range(1,m),1), 0.0, xi);

    err1 = beta * std::abs(H(m+1,1));
    sci::dgemv(sci::mat::N, 1.0, *Q2, V(sci::range(1,n+n),m+1), 0.0, tmpx);
    err2 = beta * std::abs(H(m+2,1)) * sci::dnrm2(tmpx);
    if (err1 > 10.0*err2) {
      err = err2;
    } else if (err1 > err2) {
      err = err1*err2 / (err1 - err2);
    } else {
      err = err1;
    }

    y = xi(sci::range(1,n));
    cy += xi(sci::range(n+1,n+n));

    delete Q2;
    return y;
  }

  // add elements
  sci::spmatrix<double>* addelem(sci::mat::trans trans, const sci::spmatrix<double>& Q) {
    switch (Q.type()) {
    case (CSR):
        return addelem(trans, dynamic_cast<const sci::csrmatrix<double>&>(Q));
    case (CSC):
        return addelem(trans, dynamic_cast<const sci::cscmatrix<double>&>(Q));
    case (COO):
        return addelem(trans, dynamic_cast<const sci::coomatrix<double>&>(Q));
    default:
        throw;
    }
  }

  sci::spmatrix<double>* addelem(sci::mat::trans trans, const sci::csrmatrix<double>& Q) {
    int n = Q.nrow;
    sci::csrmatrix<double>* Q2 = new sci::csrmatrix<double>(n+n, n+n, Q.nnz+n);
    sci::array<int>::const_range rowptr = Q.get_rowptr();
    sci::array<int>::const_range colind = Q.get_colind();
    sci::array<double>::const_range value = Q.get_value();
    sci::array<int>::range rowptr2 = Q2->get_rowptr();
    sci::array<int>::range colind2 = Q2->get_colind();
    sci::array<double>::range value2 = Q2->get_value();

    if (trans == sci::mat::N) {
      int z2 = 1;
      for (int i=1; i<=n+n; i++) {
        rowptr2[i] = z2;
        if (i <= n) {
          for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
            colind2[z2] = colind[z];
            value2[z2] = value[z];
            z2++;
          }
        } else {
            colind2[z2] = i-n;
            value2[z2] = 1.0;
            z2++;
        }
      }
      rowptr2[n+n+1] = z2;
    } else {
      int z2 = 1;
      for (int i=1; i<=n+n; i++) {
        rowptr2[i] = z2;
        if (i <= n) {
          for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
            colind2[z2] = colind[z];
            value2[z2] = value[z];
            z2++;
          }
          colind2[z2] = i + n;
          value2[z2] = 1.0;
          z2++;
        }
      }
      rowptr2[n+n+1] = z2;
    }

    return Q2;
  }

  sci::spmatrix<double>* addelem(sci::mat::trans trans, const sci::cscmatrix<double>& Q) {
    int n = Q.nrow;
    sci::cscmatrix<double>* Q2 = new sci::cscmatrix<double>(n+n, n+n, Q.nnz+n);
    sci::array<int>::const_range rowptr = Q.get_colptr();
    sci::array<int>::const_range colind = Q.get_rowind();
    sci::array<double>::const_range value = Q.get_value();
    sci::array<int>::range rowptr2 = Q2->get_colptr();
    sci::array<int>::range colind2 = Q2->get_rowind();
    sci::array<double>::range value2 = Q2->get_value();

    if (trans == sci::mat::T) {
      int z2 = 1;
      for (int i=1; i<=n+n; i++) {
        rowptr2[i] = z2;
        if (i <= n) {
          for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
            colind2[z2] = colind[z];
            value2[z2] = value[z];
            z2++;
          }
        } else {
            colind2[z2] = i-n;
            value2[z2] = 1.0;
            z2++;
        }
      }
      rowptr2[n+n+1] = z2;
    } else {
      int z2 = 1;
      for (int i=1; i<=n+n; i++) {
        rowptr2[i] = z2;
        if (i <= n) {
          for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
            colind2[z2] = colind[z];
            value2[z2] = value[z];
            z2++;
          }
          colind2[z2] = i + n;
          value2[z2] = 1.0;
          z2++;
        }
      }
      rowptr2[n+n+1] = z2;
    }

    return Q2;
  }

  sci::spmatrix<double>* addelem(sci::mat::trans trans, const sci::coomatrix<double>& Q) {
    int n = Q.nrow;
    sci::coomatrix<double>* Q2 = new sci::coomatrix<double>(n+n, n+n, Q.nnz+n);
    sci::array<int>::const_range rowind = Q.get_rowind();
    sci::array<int>::const_range colind = Q.get_colind();
    sci::array<double>::const_range value = Q.get_value();
    sci::array<int>::range rowind2 = Q2->get_rowind();
    sci::array<int>::range colind2 = Q2->get_colind();
    sci::array<double>::range value2 = Q2->get_value();

    if (trans == sci::mat::N) {
      for (int z=1; z<=Q.nnz; z++) {
        rowind2[z] = rowind[z];
        colind2[z] = colind[z];
        value2[z] = value[z];
      }
      for (int z=Q.nnz+1, i=1; z<=Q.nnz+n; z++, i++) {
        rowind2[z] = n + i;
        colind2[z] = i;
        value2[z] = 1.0;
      }
    } else {
      for (int z=1; z<=Q.nnz; z++) {
        rowind2[z] = rowind[z];
        colind2[z] = colind[z];
        value2[z] = value[z];
      }
      for (int z=Q.nnz+1, i=1; z<=Q.nnz+n; z++, i++) {
        rowind2[z] = i;
        colind2[z] = n + i;
        value2[z] = 1.0;
      }
    }

    return Q2;
  }
}

