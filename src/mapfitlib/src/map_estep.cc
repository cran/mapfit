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

#include "poisson.hh"
#include "mexp.hh"
#include "mapfit.hh"

#define MAX(a,b) ((a) > (b) ? (a) : (b))

namespace mapfit {

	sci::matrix<double>* new_matrix(const sci::matrix<double>& m);

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
		double poi_eps) {

		int n = alpha.size;
		int m = tdat.size;

		double scale;
		double llf = 0.0;
		sci::vector<double> tmpv(n);
		sci::vector<double> tmpz(n);
		sci::array< sci::vector<double> > vf(m+2, sci::vector<double>(n));
		sci::array< sci::vector<double> > vb(m+2, sci::vector<double>(n));

		// std::unique_ptr< sci::matrix<double> > hen0(new_matrix(en0));
		// std::unique_ptr< sci::matrix<double> > hen1(new_matrix(en1));
		sci::matrix<double>* hen0 = new_matrix(en0);
		sci::matrix<double>* hen1 = new_matrix(en1);

		eb = 0.0;
		ez = 0.0;
		sci::dfill(0.0, en0);
		sci::dfill(0.0, en1);

		double weight;
		double tmax = dmax(tdat);
		int right = pois::rightbound(qv*tmax, poi_eps) + 1;
		int nmax = gdat(imax(gdat)) + 1;
		sci::range r(0, MAX(right, nmax));
		sci::vector<double> poi(r.size());
//		std::unique_ptr<double[]> work(new double[n*(nmax+1)*(right+2)]);
		double* work = new double[n*(nmax+1)*(right+2)];

		vb[m+1] = xi;
		for (int k=m; k>=1; k--) {
			right = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			r.end = MAX(right, gdat(k)+1);
			weight = pois::pmf(qv*tdat(k), r, poi);

			if (idat(k) == 1) {
				sci::dgemv(sci::mat::N, 1.0, D1, vb[k+1], 0.0, tmpv);
			} else {
				tmpv = vb[k+1];
			}
			mapblas::mexp_unifvec_nbackward(P0, P1, qv, r.end, poi, weight,
				gdat(k), tmpv, vb[k], work);
			scale = sci::dsum(vb[k]);
			sci::dscal(1.0/scale, vb[k]);
			llf += log(scale);
		}
		eb = vb[1];

		vf[0] = alpha;
		for (int k=1; k<=m; k++) {
			right = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			r.end = MAX(right, gdat(k)+1);
			weight = pois::pmf(qv*tdat(k), r, poi);

			mapblas::mexp_unifvec_nforward(P0, P1, qv, r.end, poi, weight,
				gdat(k), vf[k-1], tmpv, work);
			if (idat(k) == 1) {
				sci::dgemv(sci::mat::T, 1.0, D1, tmpv, 0.0, vf[k]);
			} else {
				vf[k] = tmpv;
			}
			scale = sci::dsum(vf[k]);
			sci::dscal(1.0/scale, vf[k]);
		}

		for (int k=1; k<=m; k++) {
			right = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			r.end = MAX(right, gdat(k)+1);
			weight = pois::pmf(qv*tdat(k), r, poi);

			if (idat(k) == 1) {
				sci::dgemv(sci::mat::N, 1.0, D1, vb[k+1], 0.0, tmpv);
			} else {
				tmpv = vb[k+1];
			}
			sci::dfill(0.0, *hen0);
			sci::dfill(0.0, *hen1);
			mapblas::mexpc_unifvec_nforward(P0, P1, qv, r.end, poi, weight,
				gdat(k), vf[k-1], tmpv, tmpz, *hen0, *hen1, work);
			scale = sci::ddot(tmpz, tmpv);
			sci::daxpy(1.0/scale, *hen0, en0);
			sci::daxpy(1.0/scale, *hen1, en1);

			if (idat(k) == 1) {
				sci::dger(1.0/scale, tmpz, vb[k+1], en1);
			}
		}
		eb *= alpha;
		scale = dsum(eb);
		dscal(1.0/scale, eb);
		llf += log(scale);
		ez = sci::diag(en0);
		en0 *= D0;
		en1 *= D1;

		delete hen0;
		delete hen1;
		delete [] work;

		return llf;
	}

	sci::matrix<double>* new_matrix(const sci::matrix<double>& m) {
        switch (m.type()) {
        case (DENSE):
            return new sci::dmatrix<double>(dynamic_cast<const sci::dmatrix<double>&>(m));
        case (CSR):
            return new sci::csrmatrix<double>(dynamic_cast<const sci::csrmatrix<double>&>(m));
        case (CSC):
            return new sci::cscmatrix<double>(dynamic_cast<const sci::cscmatrix<double>&>(m));
        case (COO):
            return new sci::coomatrix<double>(dynamic_cast<const sci::coomatrix<double>&>(m));
        default:
            throw;
        }		
	}

	// experimental
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
		double poi_eps) {

		int n = alpha.size;
		int m = tdat.size;

		double scale;
		double llf = 0.0;
		sci::vector<double> tmpv(n);
		sci::vector<double> tmpz(n);
		sci::array< sci::vector<double> > vf(m+2, sci::vector<double>(n));
		sci::array< sci::vector<double> > vb(m+2, sci::vector<double>(n));

		// std::unique_ptr< sci::matrix<double> > hen0(new_matrix(en0));
		// std::unique_ptr< sci::matrix<double> > hen1(new_matrix(en1));
		sci::matrix<double>* hen0 = new_matrix(en0);
		sci::matrix<double>* hen1 = new_matrix(en1);

		eb = 0.0;
		ez = 0.0;
		sci::dfill(0.0, en0);
		sci::dfill(0.0, en1);

		double weight;
		double tmax = dmax(tdat);
		int right = pois::rightbound(qv*tmax, poi_eps) + 1;
		int nmax = gdat(imax(gdat)) + 1;
		sci::range r(0, MAX(right, nmax));
		sci::vector<double> poi(r.size());
//		std::unique_ptr<double[]> work(new double[n*(nmax+1)*(right+2)]);
		double* work = new double[n*(nmax+1)*(right+2)];

		vb[m+1] = xi;
		for (int k=m; k>=1; k--) {
			right = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			r.end = MAX(right, gdat(k)+1);
			weight = pois::pmf(qv*tdat(k), r, poi);

			if (idat(k) == 1) {
				sci::dgemv(sci::mat::N, 1.0, D1, vb[k+1], 0.0, tmpv);
			} else {
				tmpv = vb[k+1];
			}

			if (gdat(k) >= 0) {
				mapblas::mexp_unifvec_nbackward(P0, P1, qv, r.end, poi, weight,
					gdat(k), tmpv, vb[k], work);
			} else {
				mapblas::mexp_unifvec_NAbackward(P0, P1, qv, r.end, poi, weight,
					tmpv, vb[k]);
			}

			scale = sci::dsum(vb[k]);
			sci::dscal(1.0/scale, vb[k]);
			llf += log(scale);
		}
		eb = vb[1];

		vf[0] = alpha;
		for (int k=1; k<=m; k++) {
			right = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			r.end = MAX(right, gdat(k)+1);
			weight = pois::pmf(qv*tdat(k), r, poi);

			if (gdat(k) >= 0) {
				mapblas::mexp_unifvec_nforward(P0, P1, qv, r.end, poi, weight,
					gdat(k), vf[k-1], tmpv, work);
			} else {
				mapblas::mexp_unifvec_NAforward(P0, P1, qv, r.end, poi, weight,
					vf[k-1], tmpv);
			}

			if (idat(k) == 1) {
				sci::dgemv(sci::mat::T, 1.0, D1, tmpv, 0.0, vf[k]);
			} else {
				vf[k] = tmpv;
			}

			scale = sci::dsum(vf[k]);
			sci::dscal(1.0/scale, vf[k]);
		}

		for (int k=1; k<=m; k++) {
			right = pois::rightbound(qv*tdat(k), poi_eps) + 1;
			r.end = MAX(right, gdat(k)+1);
			weight = pois::pmf(qv*tdat(k), r, poi);

			if (idat(k) == 1) {
				sci::dgemv(sci::mat::N, 1.0, D1, vb[k+1], 0.0, tmpv);
			} else {
				tmpv = vb[k+1];
			}
			sci::dfill(0.0, *hen0);
			sci::dfill(0.0, *hen1);
			if (gdat(k) >= 0) {
				mapblas::mexpc_unifvec_nforward(P0, P1, qv, r.end, poi, weight,
					gdat(k), vf[k-1], tmpv, tmpz, *hen0, *hen1, work);
			} else {
				mapblas::mexpc_unifvec_NAforward(P0, P1, qv, r.end, poi, weight,
					vf[k-1], tmpv, tmpz, *hen0, *hen1, work);
			}
			scale = sci::ddot(tmpz, tmpv);
			sci::daxpy(1.0/scale, *hen0, en0);
			sci::daxpy(1.0/scale, *hen1, en1);

			if (idat(k) == 1) {
				sci::dger(1.0/scale, tmpz, vb[k+1], en1);
			}
		}
		eb *= alpha;
		scale = dsum(eb);
		dscal(1.0/scale, eb);
		llf += log(scale);
		ez = sci::diag(en0);
		en0 *= D0;
		en1 *= D1;

		delete hen0;
		delete hen1;
		delete [] work;

		return llf;
	}
}
