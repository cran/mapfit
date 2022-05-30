/*
 phfitlib
*/

#pragma once

#include <cfloat>
#include "sci_spblas.h"

namespace mapfit {

	void phase_annealing(double K, sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::dmatrix<double>& Q, double& theta, double& c);

	void phase_annealing(double K, sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::csrmatrix<double>& Q, double& theta, double& c);

	void phase_annealing(double K, sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::cscmatrix<double>& Q, double& theta, double& c);

	void phase_annealing(double K, sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::coomatrix<double>& Q, double& theta, double& c);

	// estep for genph
	double phase_estep_wtime(const sci::vector<double>& alpha,
		const sci::vector<double>& xi, const sci::matrix<double>& Q,
		const sci::matrix<double>& P, double qv,
		const sci::vector<double>& tdat,
		const sci::vector<double>& wdat,
		double& etotal,
		sci::vector<double>& eb,
		sci::vector<double>& ey,
		sci::vector<double>& ez,
		sci::matrix<double>& en,
		double poi_eps = DBL_EPSILON);

	double phase_estep_group_trunc(
		const sci::vector<double>& alpha,
		const sci::vector<double>& baralpha,
		const sci::vector<double>& xi,
		const sci::vector<double>& one,
		const sci::matrix<double>& Q,
		const sci::matrix<double>& P, double qv,
		const sci::vector<double>& tdat,
		const sci::vector<int>& gdat,
		int gdatlast,
		const sci::vector<int>& idat,
		double& etotal,
		sci::vector<double>& eb,
		sci::vector<double>& ey,
		sci::vector<double>& ez,
		sci::matrix<double>& en,
		double poi_eps = DBL_EPSILON);

	// mstep for genph
	void phase_mstep(const double& etotal, const sci::vector<double>& eb, const sci::vector<double>& ey,
		const sci::vector<double>& ez, const sci::dmatrix<double>& en, sci::vector<double>& alpha,
		sci::vector<double>& xi, sci::dmatrix<double>& Q);

	void phase_mstep(const double& etotal, const sci::vector<double>& eb, const sci::vector<double>& ey,
		const sci::vector<double>& ez, const sci::csrmatrix<double>& en, sci::vector<double>& alpha,
		sci::vector<double>& xi, sci::csrmatrix<double>& Q);

	void phase_mstep(const double& etotal, const sci::vector<double>& eb, const sci::vector<double>& ey,
		const sci::vector<double>& ez, const sci::cscmatrix<double>& en, sci::vector<double>& alpha,
		sci::vector<double>& xi, sci::cscmatrix<double>& Q);

	void phase_mstep(const double& etotal, const sci::vector<double>& eb, const sci::vector<double>& ey,
		const sci::vector<double>& ez, const sci::coomatrix<double>& en, sci::vector<double>& alpha,
		sci::vector<double>& xi, sci::coomatrix<double>& Q);

	void phase_mstep(const double& etotal, const sci::vector<double>& eb, const sci::vector<double>& ey,
		const sci::vector<double>& ez, const sci::matrix<double>& en, sci::vector<double>& alpha,
		sci::vector<double>& xi, sci::matrix<double>& Q);

	void phase_bidiag_to_cf1(sci::vector<double>& alpha, sci::vector<double>& xi, sci::matrix<double>& Q);
	void phase_bidiag_to_cf1(sci::vector<double>& alpha, sci::vector<double>& xi, sci::dmatrix<double>& Q);
	void phase_bidiag_to_cf1(sci::vector<double>& alpha, sci::vector<double>& xi, sci::csrmatrix<double>& Q);
	void phase_bidiag_to_cf1(sci::vector<double>& alpha, sci::vector<double>& xi, sci::cscmatrix<double>& Q);
	void phase_bidiag_to_cf1(sci::vector<double>& alpha, sci::vector<double>& xi, sci::coomatrix<double>& Q);

	// estep for Erlang ph
	double phase_erlang_estep_wtime(
		const sci::vector<double>& alpha,
		const sci::vector<int>& shape,
		const sci::vector<double>& rate,
		const sci::vector<double>& tdat,
		const sci::vector<double>& wdat,
		double& etotal,
		sci::vector<double>& eb,
		sci::vector<double>& ew);

	double phase_erlang_estep_group_trunc(
		const sci::vector<double>& alpha,
		const sci::vector<int>& shape,
		const sci::vector<double>& rate,
		const sci::vector<double>& tdat,
		const sci::vector<int>& gdat,
		int gdatlast,
		const sci::vector<int>& idat,
		double& etotal,
		sci::vector<double>& eb,
		sci::vector<double>& ew);

	void phase_erlang_mstep(const double& etotal, const sci::vector<double>& eb, const sci::vector<double>& ew,
		sci::vector<double>& alpha, sci::vector<int>& shape, sci::vector<double>& rate);
}
