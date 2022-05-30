/*
  subroutine for map
*/

#include <cmath>
#include <memory>

#include "poisson.h"
#include "mexp.h"

namespace mapblas {

    // ! Description: make DTMC probability transition matrix from
    // !              an infinitesimal CTMC generators of MAP
    // !                 P0 = I + D0/qv, P1 = D1/qv
    // !                 qv = ufact * max|D|

  double unif(const sci::matrix<double>& D0, const sci::matrix<double>& D1,
    sci::matrix<double>& P0, sci::matrix<double>& P1, double ufact) {
    double qv = mexp::unif(D0, P0, ufact);
    sci::dfill(0.0, P1);
    sci::daxpy(1.0/qv, D1, P1);
    return qv;
  }

    // ! Description: vector-matrix operation for matrix exp form;
    // !
    // !        y = x * exp(D*t)
    // !
    // !        where D is a MAP kernel, i.e,
    // !
    // !              | Q0 Q1          |
    // !              |    Q0 Q1       |
    // !          D = |       Q0 Q1    |
    // !              |          .. .. |
    // !              |             Q0 |
    // !
    // !        The size of D is given by u-by-u structured matrix,
    // !        and t is involved in the Poisson probability vector.
    // !
    // ! Parameters
    // !      n: size of CTMC (DTMC) kernel
    // !      nnz: the number of non-zeros in CTMC kernel
    // !      Q0: CTMC kernel (sparse matrix, phase matrix of MAP)
    // !      rowptr0: row pointer vector of Q0
    // !      colind0: column index vector of Q0
    // !      Q1: CTMC kernel (sparse matrix, rate matrix of MAP)
    // !      rowptr1: row pointer vector of Q1
    // !      colind1: column index vector of Q1
    // !      poi: Poisson probability
    // !      right: right bound of Poisson range
    // !      weight: normalizing constnt of a vector poi
    // !      u: the number of arrivals
    // !      x: input vector
    // !      y: output vector

  sci::vector<double>& mexp_unifvec_nforward(
    const sci::matrix<double>& P0, const sci::matrix<double>& P1, double qv,
    int right, const sci::vector<double>& poivec, double weight,
    int u, const sci::vector<double>& x, sci::vector<double>& y, double *work) {

    int n = x.size;
    sci::vector<double>::const_range poi = poivec.alias(sci::range(0,right));
    sci::vector<double> tmp(n);
    // sci::array< sci::vector<double> > xi(u+1, sci::vector<double>(n));
    sci::array< sci::vector<double> > xi(u+1, sci::vector<double>(n, work));
    double *p = work;
    for (int j=0; j<=u; j++, p+=n) {
      xi[j].ptr = p;
    }

    xi[0] = x;
    for (int j=1; j<=u; j++) {
      xi[j] = 0.0;
    }

    y = 0.0;
    sci::daxpy(poi(0), xi[u], y);
    for (int k=1; k<=right; k++) {
      for (int j=u; j>=0; j--) {
        sci::dgemv(sci::mat::T, 1.0, P0, xi[j], 0.0, tmp);
        if (j != 0) {
          sci::dgemv(sci::mat::T, 1.0, P1, xi[j-1], 1.0, tmp);
        }
        xi[j] = tmp;
      }
      sci::daxpy(poi(k), xi[u], y);
    }
    sci::dscal(1.0/weight, y);
    return y;
  }

  sci::vector<double>& mexp_unifvec_NAforward(
    const sci::matrix<double>& P0, const sci::matrix<double>& P1, double qv,
    int right, const sci::vector<double>& poivec, double weight,
    const sci::vector<double>& x, sci::vector<double>& y) {

    int n = x.size;
    sci::vector<double>::const_range poi = poivec.alias(sci::range(0,right));
    sci::vector<double> tmp(n);
    sci::vector<double> xi(x);
    y = 0.0;
    sci::daxpy(poi(0), xi, y);
    for (int k=1; k<=right; k++) {
        sci::dgemv(sci::mat::T, 1.0, P0, xi, 0.0, tmp);
        sci::dgemv(sci::mat::T, 1.0, P1, xi, 1.0, tmp);
        xi = tmp;
        sci::daxpy(poi(k), xi, y);
    }
    sci::dscal(1.0/weight, y);
    return y;
  }

    // ! Description: vector-matrix operation for matrix exp form;
    // !
    // !        y = exp(D*t) * x
    // !
    // !        where D is a MAP kernel, i.e,
    // !
    // !              | Q0 Q1          |
    // !              |    Q0 Q1       |
    // !          D = |       Q0 Q1    |
    // !              |          .. .. |
    // !              |             Q0 |
    // !
    // !        The size of D is given by u-by-u structured matrix,
    // !        and t is involved in the Poisson probability vector.
    // !
    // ! Parameters
    // !      n: size of CTMC (DTMC) kernel
    // !      nnz: the number of non-zeros in CTMC kernel
    // !      Q0: CTMC kernel (sparse matrix, phase matrix of MAP)
    // !      rowptr0: row pointer vector of Q0
    // !      colind0: column index vector of Q0
    // !      Q1: CTMC kernel (sparse matrix, rate matrix of MAP)
    // !      rowptr1: row pointer vector of Q1
    // !      colind1: column index vector of Q1
    // !      poi: Poisson probability
    // !      right: right bound of Poisson range
    // !      weight: normalizing constnt of a vector poi
    // !      u: the number of arrivals
    // !      x: input vector
    // !      y: output vector

 sci::vector<double>& mexp_unifvec_nbackward(
    const sci::matrix<double>& P0, const sci::matrix<double>& P1, double qv,
    int right, const sci::vector<double>& poivec, double weight,
    int u, const sci::vector<double>& x, sci::vector<double>& y, double *work) {

    int n = x.size;
    sci::vector<double>::const_range poi = poivec.alias(sci::range(0,right));
    sci::vector<double> tmp(n);
    sci::array< sci::vector<double> > xi(u+1, sci::vector<double>(n, work));
    double *p = work;
    for (int j=0; j<=u; j++, p+=n) {
      xi[j].ptr = p;
    }

    for (int j=0; j<=u-1; j++) {
      xi[j] = 0.0;
    }
    xi[u] = x;

    y = 0.0;
    sci::daxpy(poi(0), xi[0], y);
    for (int k=1; k<=right; k++) {
      for (int j=0; j<=u; j++) {
        sci::dgemv(sci::mat::N, 1.0, P0, xi[j], 0.0, tmp);
        if (j != u) {
          sci::dgemv(sci::mat::N, 1.0, P1, xi[j+1], 1.0, tmp);
        }
        xi[j] = tmp;
      }
      sci::daxpy(poi(k), xi[0], y);
    }
    sci::dscal(1.0/weight, y);
    return y;
  }

  sci::vector<double>& mexp_unifvec_NAbackward(
    const sci::matrix<double>& P0, const sci::matrix<double>& P1, double qv,
    int right, const sci::vector<double>& poivec, double weight,
    const sci::vector<double>& x, sci::vector<double>& y) {

    int n = x.size;
    sci::vector<double>::const_range poi = poivec.alias(sci::range(0,right));
    sci::vector<double> tmp(n);
    sci::vector<double> xi(x);
    y = 0.0;
    sci::daxpy(poi(0), xi, y);
    for (int k=1; k<=right; k++) {
        sci::dgemv(sci::mat::N, 1.0, P0, xi, 0.0, tmp);
        sci::dgemv(sci::mat::N, 1.0, P1, xi, 1.0, tmp);
        xi = tmp;
        sci::daxpy(poi(k), xi, y);
    }
    sci::dscal(1.0/weight, y);
    return y;
  }

    // ! Description: convolution integral operation for matrix exp form;
    // !
    // !               u   |t
    // !        MH0 = Sum  | exp(D(i)*s) * vf * vb * exp(D(u-i)*(t-s)) ds
    // !              i=0  |0
    // !
    // !              u-1  |t
    // !        MH1 = Sum  | exp(D(i)*s) * vf * vb * exp(D(u-i-1)*(t-s)) ds
    // !              i=0  |0
    // !
    // !        and
    // !
    // !        vf = vf * exp(D(u)*t)
    // !
    // !        where D is a MAP kernel, i.e,
    // !
    // !                 | Q0 Q1          |
    // !                 |    Q0 Q1       |
    // !          D(i) = |       Q0 Q1    |
    // !                 |          .. .. |
    // !                 |             Q0 |
    // !
    // !        The size of D(i) is given by i-by-i structured matrix,
    // !        and t is involved in the Poisson probability vector.
    // !
    // ! Parameters
    // !      n: size of CTMC (DTMC) kernel
    // !      nnz: the number of non-zeros in CTMC kernel
    // !      Q0: CTMC kernel (sparse matrix, phase matrix of MAP)
    // !      rowptr0: row pointer vector of Q0
    // !      colind0: column index vector of Q0
    // !      Q1: CTMC kernel (sparse matrix, rate matrix of MAP)
    // !      rowptr1: row pointer vector of Q1
    // !      colind1: column index vector of Q1
    // !      poi: Poisson probability
    // !      right: right bound of Poisson range
    // !      weight: normalizing constnt of a vector poi
    // !      u: the number of arrivals
    // !      vf: forward vector (inout)
    // !      vb: backward vector (in)
    // !      MH0: convint result (out), (sojourn time and phase transition)
    // !      MH1: convint result (out), (arrial transition)

  sci::vector<double>& mexpc_unifvec_nforward(
    const sci::matrix<double>& P0,
    const sci::matrix<double>& P1,
    double qv,
    int right,
    const sci::vector<double>& poivec,
    double weight,
    int u,
    const sci::vector<double>& x,
    const sci::vector<double>& y,
    sci::vector<double>& z,
    sci::matrix<double>& H0,
    sci::matrix<double>& H1,
    double *work) {

    sci::vector<double>::const_range poi = poivec.alias(sci::range(0,right));
    int n = x.size;
    sci::vector<double> tmp(n);

    sci::array< sci::vector<double> > xi(u+1, sci::vector<double>(n, work));
    sci::array< sci::array< sci::vector<double> > >
      vc(right+1, sci::array< sci::vector<double> >(u+1, sci::vector<double>(n, work)));

    double *p = work;
    for (int j=0; j<=u; j++, p+=n) {
      xi[j].ptr = p;
    }
    for (int l=0; l<=right; l++) {
      for (int j=0; j<=u; j++, p+=n) {
        vc[l][j].ptr = p;
      }
    }

    // clear
    for (int l=0; l<=right; l++) {
      for (int j=0; j<=u; j++) {
        vc[l][j] = 0.0;
      }
    }

    // forward and backward
    sci::daxpy(poi(right), y, vc[right][u]);
    for (int j=0; j<=u-1; j++) {
      vc[right][j] = 0.0;
    }

    for (int l=right-1; l>=1; l--) {
      for (int j=0; j<=u; j++) {
        sci::dgemv(sci::mat::N, 1.0, P0, vc[l+1][j], 0.0, vc[l][j]);
        if (j != u) {
          sci::dgemv(sci::mat::N, 1.0, P1, vc[l+1][j+1], 1.0, vc[l][j]);
        }
      }
      sci::daxpy(poi(l), y, vc[l][u]);
    }

    // compute H
    xi[0] = x;
    for (int j=1; j<=u; j++) {
      xi[j] = 0.0;
    }
    z = 0.0;
    // H = 0.0;
    sci::dscal(qv*weight, H0); // original H is scaled
    sci::dscal(qv*weight, H1); // original H is scaled
    sci::daxpy(poi(0), xi[u], z);

    for (int j=0; j<=u; j++) {
      sci::dger(1.0, xi[j], vc[1][j], H0);
      if (j != u) {
        sci::dger(1.0, xi[j], vc[1][j+1], H1);
      }
    }

    for (int l=1; l<=right-1; l++) {
      for (int j=u; j>=0; j--) {
        tmp = xi[j];
        sci::dgemv(sci::mat::T, 1.0, P0, tmp, 0.0, xi[j]);
        if (j != 0) {
          sci::dgemv(sci::mat::T, 1.0, P1, xi[j-1], 1.0, xi[j]);
        }
      }
      sci::daxpy(poi(l), xi[u], z);
      for (int j=0; j<=u; j++) {
        sci::dger(1.0, xi[j], vc[l+1][j], H0);
        if (j != u) {
          sci::dger(1.0, xi[j], vc[l+1][j+1], H1);
        }
      }
    }
    sci::dscal(1.0/weight, z);
    sci::dscal(1.0/qv/weight, H0);
    sci::dscal(1.0/qv/weight, H1);
    return z;
  }

  sci::vector<double>& mexpc_unifvec_NAforward(
    const sci::matrix<double>& P0,
    const sci::matrix<double>& P1,
    double qv,
    int right,
    const sci::vector<double>& poivec,
    double weight,
    const sci::vector<double>& x,
    const sci::vector<double>& y,
    sci::vector<double>& z,
    sci::matrix<double>& H0,
    sci::matrix<double>& H1,
    double *work) {

    sci::vector<double>::const_range poi = poivec.alias(sci::range(0,right));
    int n = x.size;

    sci::array< sci::vector<double> > vc(right+1, sci::vector<double>(n, work));
    double *p = work;
    for (int l=0; l<=right; l++, p+=n) {
        vc[l].ptr = p;
    }
    // clear
    for (int l=0; l<=right; l++) {
        vc[l] = 0.0;
    }

    // forward and backward
    sci::daxpy(poi(right), y, vc[right]);
    for (int l=right-1; l>=1; l--) {
        sci::dgemv(sci::mat::N, 1.0, P0, vc[l+1], 0.0, vc[l]);
        sci::dgemv(sci::mat::N, 1.0, P1, vc[l+1], 1.0, vc[l]);
        sci::daxpy(poi(l), y, vc[l]);
    }

    // compute H
    sci::vector<double> tmp(n);
    sci::vector<double> xi(x);
    z = 0.0;
    sci::dscal(qv*weight, H0); // original H is scaled
    sci::dscal(qv*weight, H1); // original H is scaled
    sci::daxpy(poi(0), xi, z);
    sci::dger(1.0, xi, vc[1], H0);
    sci::dger(1.0, xi, vc[1], H1);
    for (int l=1; l<=right-1; l++) {
        tmp = xi;
        sci::dgemv(sci::mat::T, 1.0, P0, tmp, 0.0, xi);
        sci::dgemv(sci::mat::T, 1.0, P1, tmp, 1.0, xi);
        sci::daxpy(poi(l), xi, z);
        sci::dger(1.0, xi, vc[l+1], H0);
        sci::dger(1.0, xi, vc[l+1], H1);
    }
    sci::dscal(1.0/weight, z);
    sci::dscal(1.0/qv/weight, H0);
    sci::dscal(1.0/qv/weight, H1);
    return z;
  }

    // ! Description: convolution integral operation for matrix exp form;
    // !
    // !               u   |t
    // !        MH0 = Sum  | exp(D(i)*s) * vf * vb * exp(D(u-i)*(t-s)) ds
    // !              i=0  |0
    // !
    // !              u-1  |t
    // !        MH1 = Sum  | exp(D(i)*s) * vf * vb * exp(D(u-i-1)*(t-s)) ds
    // !              i=0  |0
    // !
    // !        and
    // !
    // !        vb = exp(D(u)*t) * vb
    // !
    // !        where D is a MAP kernel, i.e,
    // !
    // !                 | Q0 Q1          |
    // !                 |    Q0 Q1       |
    // !          D(i) = |       Q0 Q1    |
    // !                 |          .. .. |
    // !                 |             Q0 |
    // !
    // !        The size of D(i) is given by i-by-i structured matrix,
    // !        and t is involved in the Poisson probability vector.
    // !
    // ! Parameters
    // !      n: size of CTMC (DTMC) kernel
    // !      nnz: the number of non-zeros in CTMC kernel
    // !      Q0: CTMC kernel (sparse matrix, phase matrix of MAP)
    // !      rowptr0: row pointer vector of Q0
    // !      colind0: column index vector of Q0
    // !      Q1: CTMC kernel (sparse matrix, rate matrix of MAP)
    // !      rowptr1: row pointer vector of Q1
    // !      colind1: column index vector of Q1
    // !      poi: Poisson probability
    // !      right: right bound of Poisson range
    // !      weight: normalizing constnt of a vector poi
    // !      u: the number of arrivals
    // !      vf: forward vector (in)
    // !      vb: backward vector (inout)
    // !      MH0: convint result (out), (sojourn time and phase transition)
    // !      MH1: convint result (out), (arrial transition)

  sci::vector<double>& mexpc_unifvec_nbackward(
    const sci::matrix<double>& P0,
    const sci::matrix<double>& P1,
    double qv,
    int right,
    const sci::vector<double>& poivec,
    double weight,
    int u,
    const sci::vector<double>& x,
    const sci::vector<double>& y,
    sci::vector<double>& z,
    sci::matrix<double>& H0,
    sci::matrix<double>& H1,
    double *work) {

    sci::vector<double>::const_range poi = poivec.alias(sci::range(0,right));
    int n = x.size;
    sci::vector<double> tmp(n);

    sci::array< sci::vector<double> > xi(u+1, sci::vector<double>(n, work));
    sci::array< sci::array< sci::vector<double> > >
      vc(right+1, sci::array< sci::vector<double> >(u+1, sci::vector<double>(n, work)));

    double *p = work;
    for (int j=0; j<=u; j++, p+=n) {
      xi[j].ptr = p;
    }
    for (int l=0; l<=right; l++) {
      for (int j=0; j<=u; j++, p+=n) {
        vc[l][j].ptr = p;
      }
    }

    // clear
    for (int l=0; l<=right; l++) {
      for (int j=0; j<=u; j++) {
        vc[l][j] = 0.0;
      }
    }

    // forward and backward
    sci::daxpy(poi(right), x, vc[right][0]);
    for (int j=1; j<=u; j++) {
      vc[right][j] = 0.0;
    }

    for (int l=right-1; l>=1; l--) {
      for (int j=u; j>=0; j--) {
        sci::dgemv(sci::mat::T, 1.0, P0, vc[l+1][j], 0.0, vc[l][j]);
        if (j != 0) {
          sci::dgemv(sci::mat::T, 1.0, P1, vc[l+1][j-1], 1.0, vc[l][j]);
        }
      }
      sci::daxpy(poi(l), x, vc[l][0]);
    }

    // compute H
    xi[u] = y;
    for (int j=0; j<=u-1; j++) {
      xi[j] = 0.0;
    }
    z = 0.0;
    // H = 0.0;
    sci::dscal(qv*weight, H0); // original H is scaled
    sci::dscal(qv*weight, H1); // original H is scaled
    sci::daxpy(poi(0), xi[0], z);

    for (int j=0; j<=u; j++) {
      sci::dger(1.0, vc[1][j], xi[j], H0);
      if (j != u) {
        sci::dger(1.0, vc[1][j], xi[j+1], H1);
      }
    }

    for (int l=1; l<=right-1; l++) {
      for (int j=0; j<=u; j++) {
        tmp = xi[j];
        sci::dgemv(sci::mat::N, 1.0, P0, tmp, 0.0, xi[j]);
        if (j != u) {
          sci::dgemv(sci::mat::N, 1.0, P1, xi[j+1], 1.0, xi[j]);
        }
      }
      sci::daxpy(poi(l), xi[0], z);
      for (int j=0; j<=u; j++) {
        sci::dger(1.0, vc[l+1][j], xi[j], H0);
        if (j != u) {
          sci::dger(1.0, vc[l+1][j], xi[j+1], H1);
        }
      }
    }
    sci::dscal(1.0/weight, z);
    sci::dscal(1.0/qv/weight, H0);
    sci::dscal(1.0/qv/weight, H1);
    return z;
  }

  sci::vector<double>& mexpc_unifvec_NAbackward(
    const sci::matrix<double>& P0,
    const sci::matrix<double>& P1,
    double qv,
    int right,
    const sci::vector<double>& poivec,
    double weight,
    const sci::vector<double>& x,
    const sci::vector<double>& y,
    sci::vector<double>& z,
    sci::matrix<double>& H0,
    sci::matrix<double>& H1,
    double *work) {

    sci::vector<double>::const_range poi = poivec.alias(sci::range(0,right));
    int n = x.size;

    sci::array< sci::vector<double> > vc(right+1, sci::vector<double>(n, work));
    double *p = work;
    for (int l=0; l<=right; l++, p+=n) {
        vc[l].ptr = p;
    }
    // clear
    for (int l=0; l<=right; l++) {
        vc[l] = 0.0;
    }

    // forward and backward
    sci::daxpy(poi(right), x, vc[right]);
    for (int l=right-1; l>=1; l--) {
        sci::dgemv(sci::mat::T, 1.0, P0, vc[l+1], 0.0, vc[l]);
        sci::dgemv(sci::mat::T, 1.0, P1, vc[l+1], 1.0, vc[l]);
        sci::daxpy(poi(l), x, vc[l]);
    }

    // compute H
    sci::vector<double> tmp(n);
    sci::vector<double> xi(y);
    z = 0.0;
    sci::dscal(qv*weight, H0); // original H is scaled
    sci::dscal(qv*weight, H1); // original H is scaled
    sci::daxpy(poi(0), xi, z);
    sci::dger(1.0, vc[1], xi, H0);
    sci::dger(1.0, vc[1], xi, H1);
    for (int l=1; l<=right-1; l++) {
        tmp = xi;
        sci::dgemv(sci::mat::N, 1.0, P0, tmp, 0.0, xi);
        sci::dgemv(sci::mat::N, 1.0, P1, tmp, 1.0, xi);
        sci::daxpy(poi(l), xi, z);
        sci::dger(1.0, vc[l+1], xi, H0);
        sci::dger(1.0, vc[l+1], xi, H1);
    }
    sci::dscal(1.0/weight, z);
    sci::dscal(1.0/qv/weight, H0);
    sci::dscal(1.0/qv/weight, H1);
    return z;
  }
}

