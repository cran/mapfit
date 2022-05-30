/*
  Description: estep for PH with weighted time and group/truncated data

   alpha      (in): initial vector
   baralpha   (in): baralpha = alpha (-Q)^-1
   xi         (in): exit vector
   one        (in): one vector
   Q          (in): infinitesimal generator
   P          (in): uniformed generator
   qv         (in): uniformization constant
   tdat       (in): interarrival time
   wdat       (in): weights for interarrivals
   gdat       (in): # of arrivals (-1 means NA)
   gdatlast   (in): # of arrivals in [lasttime, infinity] (-1 means NA)
   idat       (in): indicator whether an arrival occurs at the last instant
   etotal    (out): expected # of arrivals
   eb        (out): expected # of starts
   ey        (out): expected # of exits
   ez        (out): expected sojourn time
   en        (out): expected # of phase transitions
   poi_eps    (in): eps for poisson prob (optional)
   atol       (in): tolerance error in uniformization (optional)

   return value -> llf (log-likelihood)

 */

#include <cmath>
#include <memory>

#include "gamma.h"

#include "poisson.h"
#include "mexp.h"
#include "phfit.h"

namespace mapfit {

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
		double poi_eps) {

		int n = Q.nrow;
		int m = tdat.size;

		double scale;
		double llf = 0.0;
		double tllf = 0.0;
		sci::array< sci::vector<double> > vf(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > vb(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > vc(m+1, sci::vector<double>(n));
		sci::vector<double> blf(m);

		etotal = 0.0;
		eb = 0.0;
		ey = 0.0;
		ez = 0.0;
		sci::dfill(0.0, en);

		double tmax = dmax(tdat);
		double weight;
		sci::range r(0, pois::rightbound(qv*tmax, poi_eps) + 1);
		sci::vector<double> poi(r.size());

		// forward & backward
		vf[0] = alpha;
		vb[0] = xi;
		for (int k=1; k<=m; k++) {
			r.end = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k), r, poi);

			mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, vf[k-1], vf[k], 0.0);
			scale = sci::ddot(vf[k], xi);
			sci::dscal(1.0/scale, vf[k]);
			sci::daxpy(wdat(k), vf[k], ey);
			blf(k) = scale;

			mexp::mexp_unifvec(sci::mat::N, P, qv, r, poi, weight, vb[k-1], vb[k], 0.0);
			scale = sci::ddot(alpha, vb[k]);
			sci::dscal(1.0/scale, vb[k]);
			sci::daxpy(wdat(k), vb[k], eb);

			etotal += wdat(k);
			tllf += log(blf(k));
			llf += wdat(k) * tllf;
		}

		// vc[m] = 0.0;
		daxpy(wdat(m)/blf(m), alpha, vc[m]);
		for (int k=m-1; k>=1; k--) {
			r.end = pois::rightbound(qv*tdat(k+1), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k+1), r, poi);

			mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, vc[k+1], vc[k], 0.0);
			sci::dscal(1.0/blf(k), vc[k]);
			daxpy(wdat(k)/blf(k), alpha, vc[k]);
		}
		for (int k=1; k<=m; k++) {
			r.end = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k), r, poi);

			mexp::mexpc_unif(sci::mat::T, sci::mat::N, P, qv, r, poi, weight,
				vc[k], vb[k-1], vb[k-1], en);
		}

		eb *= alpha;
		ez = sci::diag(en);
		en *= Q;
		ey *= xi;

		return llf;
	}

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
		double poi_eps) {

		int n = Q.nrow;
		int m = tdat.size;

		double tmp, nn, uu;
		double llf = 0.0;
		sci::vector<double> vf(n);
		sci::array< sci::vector<double> > barvf(m+1, sci::vector<double>(n));
		sci::vector<double> tildevf(n);

		sci::array< sci::vector<double> > vb(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > barvb(m+1, sci::vector<double>(n));
		sci::vector<double> tildevb(n);

		sci::vector<double> wg(m+1);
		sci::vector<double> wp(m+1);
		sci::array< sci::vector<double> > vc(m+1, sci::vector<double>(n));

		eb = 0.0;
		ey = 0.0;
		ez = 0.0;
		sci::dfill(0.0, en);

		double tmax = dmax(tdat);
		double weight;
		sci::range r(0, pois::rightbound(qv*tmax, poi_eps) + 1);
		sci::vector<double> poi(r.size());

		barvf[0] = baralpha;
		barvb[0] = one;
		vb[0] = xi;
		nn = 0.0;
		uu = 0.0;

		for (int k=1; k<=m; k++) {
			r.end = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k), r, poi);

			mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, barvf[k-1], barvf[k], 0.0);
			mexp::mexp_unifvec(sci::mat::N, P, qv, r, poi, weight, barvb[k-1], barvb[k], 0.0);
			sci::dgemv(sci::mat::N, -1.0, Q, barvb[k], 0.0, vb[k]);

			tildevf = barvf[k-1] - barvf[k];
			tildevb = barvb[k-1] - barvb[k];

			if (gdat(k) >= 0 && tdat(k) != 0.0) {
				tmp = ddot(alpha, tildevb);
				llf += gdat(k) * log(tmp) - lfact(gdat(k));
				nn += gdat(k);
				uu += tmp;
				sci::daxpy(gdat(k)/tmp, tildevb, eb);
				sci::daxpy(gdat(k)/tmp, tildevf, ey);
				wg(k) = gdat(k) / tmp;
			}
			if (idat(k) == 1) {
				sci::dgemv(sci::mat::T, -1.0, Q, barvf[k], 0.0, vf);
				tmp = sci::ddot(alpha, vb[k]);
				llf += log(tmp);
				nn += 1.0;
				sci::daxpy(1.0/tmp, vb[k], eb);
				sci::daxpy(1.0/tmp, vf, ey);
				wp(k) = 1.0 / tmp;
			}
		}
		// for the interval [t_m, infinity)
		if (gdatlast >= 0) {
			tmp = sci::ddot(alpha, barvb[m]);
			llf += gdatlast * log(tmp) - lfact(gdatlast);
			nn += gdatlast;
			uu += tmp;
			daxpy(gdatlast/tmp, barvb[m], eb);
			daxpy(gdatlast/tmp, barvf[m], ey);
			wg(m+1) = gdatlast / tmp;
		}
		// compupte weights for unobserved periods
		for (int k=1; k<=m; k++) {
			if (gdat(k) == -1) {
				daxpy(nn/uu, barvb[k-1]-barvb[k], eb);
				daxpy(nn/uu, barvf[k-1]-barvf[k], ey);
				wg(k) = nn / uu;
			}
		}
		if (gdatlast == -1) {
			daxpy(nn/uu, barvb[m], eb);
			daxpy(nn/uu, barvf[m], ey);
			wg(m+1) = nn / uu;
		}
		llf += mylgamma(nn + 1.0) - nn * log(uu);
		// compute vectors for convolution
		// vc[m] = 0.0;
		daxpy(wg(m+1)-wg(m), baralpha, vc[m]);
		if (idat(m) == 1) {
			daxpy(wp(m), alpha, vc[m]);
		}
		for (int k=m-1; k>=1; k--) {
			r.end = pois::rightbound(qv*tdat(k+1), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k+1), r, poi);

			mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, vc[k+1], vc[k], 0.0);
			daxpy(wg(k+1)-wg(k), baralpha, vc[k]);
			if (idat(k) == 1) {
				daxpy(wp(k), alpha, vc[k]);
			}
		}
		for (int k=1; k<=m; k++) {
			r.end = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k), r, poi);

			sci::dger(wg(k+1)-wg(k), baralpha, barvb[k], en);
			mexp::mexpc_unif(sci::mat::T, sci::mat::N, P, qv, r, poi, weight,
				vc[k], vb[k-1], vb[k-1], en);
		}
		sci::dger(wg(1), baralpha, barvb[0], en);

		etotal = nn / uu;
		eb *= alpha;
		ez = sci::diag(en);
		en *= Q;
		ey *= xi;

		return llf;
	}

	double phase_estep_group_trunc_poi(
		const sci::vector<double>& alpha,
		const sci::vector<double>& baralpha,
		const sci::vector<double>& xi,
		const sci::vector<double>& one,
		const sci::matrix<double>& Q,
		const sci::matrix<double>& P, double qv,
		double omega, // Poisson mean (in)
		const sci::vector<double>& tdat,
		const sci::vector<int>& gdat,
		int gdatlast,
		const sci::vector<int>& idat,
		double& etotal,
		sci::vector<double>& eb,
		sci::vector<double>& ey,
		sci::vector<double>& ez,
		sci::matrix<double>& en,
		double poi_eps) {

		int n = Q.nrow;
		int m = tdat.size;

		double tmp, nn, uu;
		double llf = 0.0;
		sci::vector<double> vf(n);
		sci::array< sci::vector<double> > barvf(m+1, sci::vector<double>(n));
		sci::vector<double> tildevf(n);

		sci::array< sci::vector<double> > vb(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > barvb(m+1, sci::vector<double>(n));
		sci::vector<double> tildevb(n);

		sci::vector<double> wg(m+1);
		sci::vector<double> wp(m+1);
		sci::array< sci::vector<double> > vc(m+1, sci::vector<double>(n));

		eb = 0.0;
		ey = 0.0;
		ez = 0.0;
		sci::dfill(0.0, en);

		double tmax = dmax(tdat);
		double weight;
		sci::range r(0, pois::rightbound(qv*tmax, poi_eps) + 1);
		sci::vector<double> poi(r.size());

		barvf[0] = baralpha;
		barvb[0] = one;
		vb[0] = xi;
		nn = 0.0;
		uu = 0.0;

		for (int k=1; k<=m; k++) {
			r.end = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k), r, poi);

			mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, barvf[k-1], barvf[k], 0.0);
			mexp::mexp_unifvec(sci::mat::N, P, qv, r, poi, weight, barvb[k-1], barvb[k], 0.0);
			sci::dgemv(sci::mat::N, -1.0, Q, barvb[k], 0.0, vb[k]);

			tildevf = barvf[k-1] - barvf[k];
			tildevb = barvb[k-1] - barvb[k];

			if (gdat(k) >= 0 && tdat(k) != 0.0) {
				tmp = ddot(alpha, tildevb);
				llf += gdat(k) * log(tmp) - lfact(gdat(k));
				nn += gdat(k);
				uu += tmp;
				sci::daxpy(gdat(k)/tmp, tildevb, eb);
				sci::daxpy(gdat(k)/tmp, tildevf, ey);
				wg(k) = gdat(k) / tmp;
			}
			if (idat(k) == 1) {
				sci::dgemv(sci::mat::T, -1.0, Q, barvf[k], 0.0, vf);
				tmp = sci::ddot(alpha, vb[k]);
				llf += log(tmp);
				nn += 1.0;
				sci::daxpy(1.0/tmp, vb[k], eb);
				sci::daxpy(1.0/tmp, vf, ey);
				wp(k) = 1.0 / tmp;
			}
		}
		// for the interval [t_m, infinity)
		if (gdatlast >= 0) {
			tmp = sci::ddot(alpha, barvb[m]);
			llf += gdatlast * log(tmp) - lfact(gdatlast);
			nn += gdatlast;
			uu += tmp;
			daxpy(gdatlast/tmp, barvb[m], eb);
			daxpy(gdatlast/tmp, barvf[m], ey);
			wg(m+1) = gdatlast / tmp;
		}
		// compupte weights for unobserved periods
		for (int k=1; k<=m; k++) {
			if (gdat(k) == -1) {
				daxpy(omega, barvb[k-1]-barvb[k], eb);
				daxpy(omega, barvf[k-1]-barvf[k], ey);
				wg(k) = omega;
			}
		}
		if (gdatlast == -1) {
			daxpy(omega, barvb[m], eb);
			daxpy(omega, barvf[m], ey);
			wg(m+1) = omega;
		}
		llf += nn * log(omega) - omega * uu;
		// compute vectors for convolution
		// vc[m] = 0.0;
		daxpy(wg(m+1)-wg(m), baralpha, vc[m]);
		if (idat(m) == 1) {
			daxpy(wp(m), alpha, vc[m]);
		}
		for (int k=m-1; k>=1; k--) {
			r.end = pois::rightbound(qv*tdat(k+1), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k+1), r, poi);

			mexp::mexp_unifvec(sci::mat::T, P, qv, r, poi, weight, vc[k+1], vc[k], 0.0);
			daxpy(wg(k+1)-wg(k), baralpha, vc[k]);
			if (idat(k) == 1) {
				daxpy(wp(k), alpha, vc[k]);
			}
		}
		for (int k=1; k<=m; k++) {
			r.end = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			weight = pois::pmf(qv*tdat(k), r, poi);

			sci::dger(wg(k+1)-wg(k), baralpha, barvb[k], en);
			mexp::mexpc_unif(sci::mat::T, sci::mat::N, P, qv, r, poi, weight,
				vc[k], vb[k-1], vb[k-1], en);
		}
		sci::dger(wg(1), baralpha, barvb[0], en);

		etotal = nn + omega * (1.0 - uu);
		eb *= alpha;
		ez = sci::diag(en);
		en *= Q;
		ey *= xi;

		return llf;
	}
}
