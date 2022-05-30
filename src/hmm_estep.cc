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

#include "gamma.h"
#include "erlang.h"

#include "mapfit.h"

namespace mapfit {

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
		sci::vector<double>& ew1) {

		int n = alpha.size;
		int m = tdat.size;

		double scale;
		double llf = 0.0;
		sci::vector<double> tmpv(n);
		sci::vector<double> hew(n);
		sci::array< sci::vector<double> > vf(m+2, sci::vector<double>(n));
		sci::array< sci::vector<double> > vb(m+2, sci::vector<double>(n));
		sci::array< sci::vector<double> > erl(m+1, sci::vector<double>(n));

		eb = 0.0;
		sci::dfill(0.0, en);
		ew0 = 0.0;
		ew1 = 0.0;

		// set Erlang pdf
		for (int k=1; k<=m; k++) {
			for (int i=1; i<=n; i++) {
				erl[k](i) = erlang_pdf(shape(i), rate(i), tdat(k));
			}
		}

		vf[0] = alpha;
		vb[m+1] = xi;
		vf[1] = vf[0];
		vf[1] *= erl[1];
		scale = dsum(vf[1]);
		dscal(1.0/scale, vf[1]);

		for (int k=2; k<=m; k++) {
			sci::dgemv(sci::mat::T, 1.0, P, vf[k-1], 0.0, vf[k]);
			vf[k] *= erl[k];
			scale = dsum(vf[k]);
			dscal(1.0/scale, vf[k]);
		}

		for (int k=m; k>=1; k--) {
			sci::dgemv(sci::mat::N, 1.0, P, vb[k+1], 0.0, vb[k]);
			vb[k] *= erl[k];
			scale = dsum(vb[k]);
			dscal(1.0/scale, vb[k]);
			llf += log(scale);
		}

		for (int k=1; k<=m; k++) {
			sci::dgemv(sci::mat::N, 1.0, P, vb[k+1], 0.0, tmpv);
			scale = ddot(vf[k], tmpv);
			sci::dger(1.0/scale, vf[k], vb[k+1], en);
			hew = vf[k];
			hew *= tmpv;
			daxpy(1.0/scale, hew, ew0);
			daxpy(tdat(k)/scale, hew, ew1);
		}

 		eb = vb[1];
 		eb *= alpha;
		scale = dsum(eb);
		dscal(1.0/scale, eb);
		llf += log(scale);
		en *= P;
		return llf;
	}

}
