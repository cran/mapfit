/*
! Description: stationary vector via GS algorithm from intial vector x
!        Q * y = 0
*/

#include <cmath>
#include "mexp.hh"

namespace mexp {

  int ctmc_gs(const sci::csrmatrix<double>& Q,
    const sci::vector<double>& x, sci::vector<double>& y, int maxiter, double eps) {

    int n = x.size;
    sci::vector<double> prevy(x);
    sci::vector<double> zero(n);
    y = x;

    for (int k=1; k<=maxiter; k++) {
      sci::dgsstep(sci::mat::T, 1.0, Q, zero, y);
      sci::dscal(1.0/sci::dsum(y), y);
      if (sci::dnrm2(y-prevy)/sci::dnrm2(y) < eps) {
        return 0;
      }
      prevy = y;
    }
//    std::cout << "Did not converge Gauss-Seidel algorithm." << std::endl;
    return -1;
  }
  
  int ctmc_gs(const sci::cscmatrix<double>& Q,
    const sci::vector<double>& x, sci::vector<double>& y, int maxiter, double eps) {

    int n = x.size;
    sci::vector<double> prevy(x);
    sci::vector<double> zero(n);
    y = x;

    for (int k=1; k<=maxiter; k++) {
      sci::dgsstep(sci::mat::T, 1.0, Q, zero, y);
      sci::dscal(1.0/sci::dsum(y), y);
      if (sci::dnrm2(y-prevy)/sci::dnrm2(y) < eps) {
        return 0;
      }
      prevy = y;
    }
//    std::cout << "Did not converge Gauss-Seidel algorithm." << std::endl;
    return -1;
  }

  int ctmc_gs(const sci::coomatrix<double>& Q,
    const sci::vector<double>& x, sci::vector<double>& y, int maxiter, double eps) {

    int n = x.size;
    sci::vector<double> prevy(x);
    sci::vector<double> zero(n);
    y = x;

    for (int k=1; k<=maxiter; k++) {
      sci::dgsstep(sci::mat::T, 1.0, Q, zero, y);
      sci::dscal(1.0/sci::dsum(y), y);
      if (sci::dnrm2(y-prevy)/sci::dnrm2(y) < eps) {
        return 0;
      }
      prevy = y;
    }
//    std::cout << "Did not converge Gauss-Seidel algorithm." << std::endl;
    return -1;
  }
}

