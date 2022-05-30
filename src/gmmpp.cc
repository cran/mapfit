/*
 *  GMAPApprox.cpp
 *  EMAP
 *
 *  Created by Hiroyuki Okamura on 06/10/08.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>

#include "gamma.h"
#include "erlang.h"

#include "pgamma.h"
#include "gauss_inte.h"

#include "mapfit.h"

 namespace mapfit {

	double xifunc0(int n, double t, double u, double mi, double mj, double ri, double rj) {
		return exp(n * log(ri*u + rj*(t-u)) - mylgamma(n+1.0) - mi*u - mj*(t-u));
	}

	void gauss_init(sci::vector<double>& x, sci::vector<double>& w, double eps) {
		int divide = x.size;
		gauss_inte_w(divide, x.ptr, w.ptr, eps);
	}

	double gam_inte(int n, double t, 
		double mi, double mj, double ri, double rj,
		const sci::vector<double>& x, const sci::vector<double>& w) {

		int divide = x.size;
		sci::vector<double> fx(divide);
		sci::vector<double> fv(divide);

		double s;
		s = gauss_inte_fx(divide, x.ptr, 0.0, t, fx.ptr);
		for (int i=1; i<=divide; i++) {
			fv(i) = xifunc0(n, t, fx(i), mi, mj, ri, rj);
		}
		return gauss_inte_fv(divide, w.ptr, s, fv.ptr);
	}

	double psi_inte(int n, double t, 
		double mi, double mj, double ri, double rj,
		const sci::vector<double>& x, const sci::vector<double>& w) {

		int divide = x.size;
		sci::vector<double> fx(divide);
		sci::vector<double> fv(divide);

		double s;
		s = gauss_inte_fx(divide, x.ptr, 0.0, t, fx.ptr);
		for (int i=1; i<=divide; i++) {
			fv(i) = fx(i) * xifunc0(n, t, fx(i), mi, mj, ri, rj);
		}
		return gauss_inte_fv(divide, w.ptr, s, fv.ptr);
	}

	void makeG(
		int num, double t,
		const sci::dmatrix<double>& D0,
		const sci::dmatrix<double>& D1,
		sci::dmatrix<double>& G,
		const sci::vector<double>& x,
		const sci::vector<double>& w) {

		double tmp, dij, dji;
		int n = D0.nrow;

		for (int j=1; j<=n; j++) {
			for (int i=1; i<=n; i++) {
				if (i == j) {
					dij = 1.0;
					dji = 1.0;
				} else {
					dij = D0(i,j);
					dji = D0(j,i);
				}
				G(i,j) = dij * gam_inte(num, t, -D0(i,i), -D0(j,j), D1(i,i), D1(j,j), x, w);
			}
		}
	}

	void makeGPsi(
		int num, double t,
		const sci::dmatrix<double>& D0,
		const sci::dmatrix<double>& D1,
		sci::dmatrix<double>& G,
		sci::dmatrix<double>& PsiT1,
		sci::dmatrix<double>& PsiT2,
		sci::dmatrix<double>& PsiN1,
		sci::dmatrix<double>& PsiN2,
		const sci::vector<double>& x,
		const sci::vector<double>& w) {

		double tmp, dij, dji;
		int n = D0.nrow;

		for (int j=1; j<=n; j++) {
			for (int i=1; i<=n; i++) {
				if (i == j) {
					dij = 1.0;
					dji = 1.0;
				} else {
					dij = D0(i,j);
					dji = D0(j,i);
				}
				G(i,j) = dij * gam_inte(num, t, -D0(i,i), -D0(j,j), D1(i,i), D1(j,j), x, w);
				tmp = psi_inte(num, t, -D0(i,i), -D0(j,j), D1(i,i), D1(j,j), x, w);
				PsiT1(i,j) = dij * tmp;
				PsiT2(j,i) = dji * tmp;
				if (num != 0) {
					tmp = psi_inte(num-1, t, -D0(i,i), -D0(j,j), D1(i,i), D1(j,j), x, w);
					PsiN1(i,j) = D1(i,i) * dij * tmp;
					PsiN2(j,i) = D1(i,i) * dji * tmp;
				}
			}
		}
	}

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
		int inte_divide, double inte_eps) {

		if (D0.type() != DENSE || D1.type() != DENSE) {
			return 0.0;
		}

		return gmmpp_estep(alpha, xi,
			dynamic_cast<const sci::dmatrix<double>&>(D0),
			dynamic_cast<const sci::dmatrix<double>&>(D1),
			tdat, gdat, idat, eb, ez,
			dynamic_cast<sci::dmatrix<double>&>(en0),
			dynamic_cast<sci::dmatrix<double>&>(en1),
			inte_divide, inte_eps);
	}

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
		int inte_divide, double inte_eps) {

		sci::vector<double> inte_x(inte_divide);
		sci::vector<double> inte_w(inte_divide);
		gauss_init(inte_x, inte_w, inte_eps);

		int n = alpha.size;
		int m = tdat.size;

		double scale;
		double llf = 0.0;
		sci::vector<double> tmpv(n);
		sci::vector<double> tmpb(n);
		sci::vector<double> tmpv2(n);
		sci::dmatrix<double> tmpm(n,n);

		sci::dmatrix<double> G(n,n);
		sci::dmatrix<double> Psi1T(n,n);
		sci::dmatrix<double> Psi2T(n,n);
		sci::dmatrix<double> Psi1N(n,n);
		sci::dmatrix<double> Psi2N(n,n);

		sci::array< sci::vector<double> > vf(m+2, sci::vector<double>(n));
		sci::array< sci::vector<double> > vb(m+2, sci::vector<double>(n));

		eb = 0.0;
		ez = 0.0;
		en0 = 0.0;
		en1 = 0.0;

		vb[m+1] = xi;
		for (int k=m; k>=1; k--) {
			if (idat(k) == 1) {
				sci::dgemv(sci::mat::N, 1.0, D1, vb[k+1], 0.0, tmpv);
			} else {
				tmpv = vb[k+1];
			}
			makeG(gdat(k), tdat(k), D0, D1, G, inte_x, inte_w);
			sci::dgemv(sci::mat::N, 1.0, G, tmpv, 0.0, vb[k]);
			scale = sci::dsum(vb[k]);
			sci::dscal(1.0/scale, vb[k]);
			llf += log(scale);
		}
		eb = vb[1];

		vf[0] = alpha;
		for (int k=1; k<=m; k++) {
			makeGPsi(gdat(k), tdat(k), D0, D1, G, Psi1T, Psi2T, Psi1N, Psi2N, inte_x, inte_w);
			sci::dgemv(sci::mat::T, 1.0, G, vf[k-1], 0.0, tmpv);
			if (idat(k) == 1) {
				sci::dgemv(sci::mat::T, 1.0, D1, tmpv, 0.0, vf[k]);
			} else {
				vf[k] = tmpv;
			}

			// for ez, en0, en1
			if (idat(k) == 1) {
				sci::dgemv(sci::mat::N, 1.0, D1, vb[k+1], 0.0, tmpb);
			} else {
				tmpb = vb[k+1];
			}
			scale = sci::ddot(tmpv, tmpb);

			// en0
			tmpm = 0.0;
			en0 += G *= sci::dger(1.0/scale, vf[k-1], tmpb, tmpm);

			// ez
			sci::dgemv(sci::mat::T, 1.0/scale, Psi2T, vf[k-1], 0.0, tmpv2);
			ez += tmpv2 *= tmpb;
			sci::dgemv(sci::mat::N, 1.0/scale, Psi1T, tmpb, 0.0, tmpv2);
			ez += tmpv2 *= vf[k-1];

			// en1
			if (gdat(k) != 0) {
				sci::dgemv(sci::mat::T, 1.0/scale, Psi2N, vf[k-1], 0.0, tmpv2);
				en1.diag() += tmpv2 *= tmpb;
				sci::dgemv(sci::mat::N, 1.0/scale, Psi1N, tmpb, 0.0, tmpv2);
				en1.diag() += tmpv2 *= vf[k-1];
			}

			scale = sci::dsum(vf[k]);
			sci::dscal(1.0/scale, vf[k]);
		}

		eb *= alpha;
		scale = dsum(eb);
		dscal(1.0/scale, eb);
		llf += log(scale);
		return llf;
	}
}

