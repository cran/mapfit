/*
mexp
*/

#pragma once

#include <cfloat>
#include "sci_spblas.h"

namespace mexp {
// util
  int pow2i(int m);

// unif
	double unif(const sci::matrix<double>& Q, sci::matrix<double>& P, double ufact = 1.01);
	// double unif(const sci::dmatrix<double>& Q, sci::dmatrix<double>& P, double ufact = 1.01);
	// double unif(const sci::csrmatrix<double>& Q, sci::csrmatrix<double>& P, double ufact = 1.01);
	// double unif(const sci::cscmatrix<double>& Q, sci::cscmatrix<double>& P, double ufact = 1.01);
	// double unif(const sci::coomatrix<double>& Q, sci::coomatrix<double>& P, double ufact = 1.01);

	sci::dmatrix<double>& mpow(const sci::dmatrix<double>& MA, sci::dmatrix<double>& ME, int m);
	sci::matrix<double>& mpow(const sci::matrix<double>& MA, sci::matrix<double>& ME, int m);

	sci::dmatrix<double>& mexp_pade(const sci::matrix<double>& MA, sci::dmatrix<double>& ME,
		double eps = DBL_EPSILON);

	sci::matrix<double>& mexp_unif(const sci::matrix<double>& P, double qv,
		const sci::range& r, const sci::vector<double>& poivec, double weight,
		const sci::dmatrix<double>& x, sci::dmatrix<double>& y,
		double atol = DBL_EPSILON);

	sci::vector<double>& mexp_unifvec(sci::mat::trans tr,
		const sci::matrix<double>& P, double qv,
		const sci::range& r, const sci::vector<double>& poivec, double weight,
		const sci::vector<double>& x, sci::vector<double>& y,
		double atol = DBL_EPSILON);

	sci::dmatrix<double>& mexpi_unif(const sci::matrix<double>& P, double qv,
		const sci::range& r, const sci::vector<double>& poivec, double weight,
		const sci::dmatrix<double>& x, sci::dmatrix<double>& y, sci::dmatrix<double>& cy);

	sci::vector<double>& mexpi_unifvec(sci::mat::trans tr,
		const sci::matrix<double>& P, double qv,
		const sci::range& r, const sci::vector<double>& poivec, double weight,
		const sci::vector<double>& x, sci::vector<double>& y, sci::vector<double>& cy,
		double atol = DBL_EPSILON);

	sci::vector<double>& mexpc_unif(
		sci::mat::trans transQ,
		sci::mat::trans transH,
		const sci::matrix<double>& P, double qv,
		const sci::range& r, const sci::vector<double>& poivec, double weight,
		const sci::vector<double>& x, const sci::vector<double>& y,
		sci::vector<double>& z, sci::matrix<double>& H);

	sci::vector<double>& mexp_unifvec_thread(sci::mat::trans trans,
		const sci::array< sci::matrix<double>* >& P, double qv,
		const sci::range& r, const sci::vector<double>& poivec, double weight,
		const sci::vector<double>& x, sci::vector<double>& y,
		double atol = DBL_EPSILON);

	sci::vector<double>& mexpi_unifvec_thread(sci::mat::trans trans,
		const sci::array< sci::matrix<double>* >& P, double qv,
		const sci::range& r, const sci::vector<double>& poivec, double weight,
		const sci::vector<double>& x, sci::vector<double>& y, sci::vector<double>& cy,
		double atol = DBL_EPSILON);

	sci::vector<double>& mexp_krylov_pade(sci::mat::trans trans,
		const sci::spmatrix<double>& Q, double t,
		const sci::vector<double>& x, sci::vector<double>& y,
		int m, double& err, int arnoldi_ite = 1, double arnoldi_tol = DBL_EPSILON,
		double pade_eps = DBL_EPSILON);

	sci::vector<double>& mexpi_krylov_pade(sci::mat::trans trans,
		const sci::spmatrix<double>& Q, double t,
		const sci::vector<double>& x, sci::vector<double>& y, sci::vector<double>& cy,
		int m, double& err, int arnoldi_ite = 1, double arnoldi_tol = DBL_EPSILON,
		double pade_eps = DBL_EPSILON);

  sci::vector<double>& ctmc_gth(const sci::dmatrix<double>& Q, sci::vector<double>& x);

  int ctmc_gs(const sci::csrmatrix<double>& Q,
    const sci::vector<double>& x, sci::vector<double>& y, int maxiter, double eps = DBL_EPSILON);

  int ctmc_gs(const sci::cscmatrix<double>& Q,
    const sci::vector<double>& x, sci::vector<double>& y, int maxiter, double eps = DBL_EPSILON);

  int ctmc_gs(const sci::coomatrix<double>& Q,
    const sci::vector<double>& x, sci::vector<double>& y, int maxiter, double eps = DBL_EPSILON);

}
