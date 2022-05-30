/*
 phfitlib
*/

#pragma once

#include <cfloat>
#include "sci_spblas.h"

namespace mapblas {

	double unif(const sci::matrix<double>& D0, const sci::matrix<double>& D1,
		sci::matrix<double>& P0, sci::matrix<double>& P1, double ufact = 1.01);	

	sci::vector<double>& mexp_unifvec_nforward(
		const sci::matrix<double>& P0, const sci::matrix<double>& P1, double qv,
		int right, const sci::vector<double>& poivec, double weight,
		int u, const sci::vector<double>& x, sci::vector<double>& y, double *work);

	sci::vector<double>& mexp_unifvec_nbackward(
		const sci::matrix<double>& P0, const sci::matrix<double>& P1, double qv,
		int right, const sci::vector<double>& poivec, double weight,
		int u, const sci::vector<double>& x, sci::vector<double>& y, double *work);

	sci::vector<double>& mexp_unifvec_NAforward(
		const sci::matrix<double>& P0, const sci::matrix<double>& P1, double qv,
		int right, const sci::vector<double>& poivec, double weight,
		const sci::vector<double>& x, sci::vector<double>& y);

	sci::vector<double>& mexp_unifvec_NAbackward(
		const sci::matrix<double>& P0, const sci::matrix<double>& P1, double qv,
		int right, const sci::vector<double>& poivec, double weight,
		const sci::vector<double>& x, sci::vector<double>& y);

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
		double *work);

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
		double *work);

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
		double *work);

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
		double *work);
}

namespace mapfit {

	double map_estep_group(
		const sci::vector<double>& alpha,
		const sci::vector<double>& xi,
		const sci::matrix<double>& D0,
		const sci::matrix<double>& D1,
		const sci::matrix<double>& P0,
		const sci::matrix<double>& P1,
		double qv,
		const sci::vector<double>& tdat,
		const sci::vector<int>& gdat,
		const sci::vector<int>& idat,
		sci::vector<double>& eb,
		sci::vector<double>& ez,
		sci::matrix<double>& en0,
		sci::matrix<double>& en1,
		double poi_eps = DBL_EPSILON);

	double map_estep_groupNA(
		const sci::vector<double>& alpha,
		const sci::vector<double>& xi,
		const sci::matrix<double>& D0,
		const sci::matrix<double>& D1,
		const sci::matrix<double>& P0,
		const sci::matrix<double>& P1,
		double qv,
		const sci::vector<double>& tdat,
		const sci::vector<int>& gdat,
		const sci::vector<int>& idat,
		sci::vector<double>& eb,
		sci::vector<double>& ez,
		sci::matrix<double>& en0,
		sci::matrix<double>& en1,
		double poi_eps = DBL_EPSILON);

	double hmm_erlang_estep(
		const sci::vector<double>& alpha,
		const sci::vector<double>& xi,
		const sci::matrix<double>& P,
		const sci::vector<int>& shape,
		const sci::vector<double>& rate,
		const sci::vector<double>& tdat,
		sci::vector<double>& eb,
		sci::matrix<double>& en,
		sci::vector<double>& ew0,
		sci::vector<double>& ew1);

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::matrix<double>& en0,
		const sci::matrix<double>& en1,
		sci::vector<double>& alpha,
		sci::matrix<double>& D0,
		sci::matrix<double>& D1);

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::dmatrix<double>& en0,
		const sci::dmatrix<double>& en1,
		sci::vector<double>& alpha,
		sci::dmatrix<double>& D0,
		sci::dmatrix<double>& D1);

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::csrmatrix<double>& en0,
		const sci::csrmatrix<double>& en1,
		sci::vector<double>& alpha,
		sci::csrmatrix<double>& D0,
		sci::csrmatrix<double>& D1);

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::cscmatrix<double>& en0,
		const sci::cscmatrix<double>& en1,
		sci::vector<double>& alpha,
		sci::cscmatrix<double>& D0,
		sci::cscmatrix<double>& D1);

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::coomatrix<double>& en0,
		const sci::coomatrix<double>& en1,
		sci::vector<double>& alpha,
		sci::coomatrix<double>& D0,
		sci::coomatrix<double>& D1);

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::matrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::matrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate);

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::dmatrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::dmatrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate);

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::csrmatrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::csrmatrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate);

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::cscmatrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::cscmatrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate);

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::coomatrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::coomatrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate);

	/// gmmpp

	double gmmpp_estep(
		const sci::vector<double>& alpha,
		const sci::vector<double>& xi,
		const sci::matrix<double>& D0,
		const sci::matrix<double>& D1,
		const sci::vector<double>& tdat,
		const sci::vector<int>& gdat,
		const sci::vector<int>& idat,
		sci::vector<double>& eb,
		sci::vector<double>& ez,
		sci::matrix<double>& en0,
		sci::matrix<double>& en1,
		int inte_divide, double inte_eps);

	double gmmpp_estep(
		const sci::vector<double>& alpha,
		const sci::vector<double>& xi,
		const sci::dmatrix<double>& D0,
		const sci::dmatrix<double>& D1,
		const sci::vector<double>& tdat,
		const sci::vector<int>& gdat,
		const sci::vector<int>& idat,
		sci::vector<double>& eb,
		sci::vector<double>& ez,
		sci::dmatrix<double>& en0,
		sci::dmatrix<double>& en1,
		int inte_divide, double inte_eps);
}
