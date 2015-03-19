/*
  ! Description: convolution integral operation for matrix exp form;
  !
  !                           |t
  ! transH(MH) = transH(MH) + | exp(transQ(Q)*s) * x * y' * exp(transQ(Q)*(t-s)) ds
  !                           |0
  !
  !        and
  !
  !        z = exp(transQ(Q)*t) * x
  !
  !        t is involved in the Poisson probability vector.
  !        qv is an uniformed parameter
  !        return value is z
 */

#include <cmath>

#include "poisson.hh"
#include "mexp.hh"

namespace mexp {

  sci::vector<double>& mexpc_unif(
    sci::mat::trans transQ,
		sci::mat::trans transH,
		const sci::matrix<double>& P, double qv,
    const sci::range& r, const sci::vector<double>& poivec, double weight,
    const sci::vector<double>& x, const sci::vector<double>& y,
    sci::vector<double>& z, sci::matrix<double>& H) {

    sci::mat::trans itransQ;
    if (transQ == sci::mat::T) {
      itransQ = sci::mat::N;
    } else if (transQ == sci::mat::N) {
      itransQ = sci::mat::T;
    } else {
      throw;
    }

    int n = P.nrow;

    sci::vector<double>::const_range poi = poivec.alias(r);
    int left = r.begin;
    int right = r.end;

    sci::array< sci::vector<double> > vcarray(right-left, sci::vector<double>(n));
    sci::array< sci::vector<double> >::range vc = vcarray.alias(sci::range(left+1,right));
    sci::daxpy(poi(right), y, vc[right]);
    for (int l=right-1; l>=left+1; l--) {
      sci::dgemv(itransQ, 1.0, P, vc[l+1], 0.0, vc[l]);
      sci::daxpy(poi(l), y, vc[l]);
    }

    sci::vector<double> tmp(n);
    sci::vector<double> xi(x);
    z = 0.0;
    // H = 0.0;
    sci::dscal(qv*weight, H); // original H is scaled
    sci::daxpy(poi(left), xi, z);
    sci::dger(transH, 1.0, xi, vc[left+1], H);
    for (int l=left+1; l<=right-1; l++) {
      tmp = xi;
      sci::dgemv(transQ, 1.0, P, tmp, 0.0, xi);
      sci::daxpy(poi(l), xi, z);
      sci::dger(transH, 1.0, xi, vc[l+1], H);
    }
    sci::dscal(1.0/weight, z);
    sci::dscal(1.0/qv/weight, H);
    return z;
  }

}
