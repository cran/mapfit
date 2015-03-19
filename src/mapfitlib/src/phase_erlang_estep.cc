/*
  Description: estep for Erlang-PH with weighted time and group/truncated data

   alpha    (in): initial vector
   shape    (in): shape parameter vector
   rate     (in): rate parameter vector
   tdat     (in): interarrival time
   wdat     (in): weights for interarrivals
   gdat     (in): # of arrivals (-1 means NA)
   gdatlast (in): # of arrivals in [lasttime, infinity] (-1 means NA)
   idat     (in): indicator whether an arrival occurs at the last instant
   etotal  (out): expected # of arrivals
   eb      (out): expected # of starts
   ew      (out): expected sojourn time?

   return value -> llf (log-likelihood)

 */

#include <cmath>

#include "gamma.h"
#include "erlang.hh"
#include "phfit.hh"

namespace mapfit {

	double phase_erlang_estep_wtime(
		const sci::vector<double>& alpha,
		const sci::vector<int>& shape,
		const sci::vector<double>& rate,
		const sci::vector<double>& tdat,
		const sci::vector<double>& wdat,
		double& etotal,
		sci::vector<double>& eb,
		sci::vector<double>& ew) {

		int n = alpha.size;
		int m = tdat.size;

		double scale, tmp;
		double llf = 0.0;
		sci::array< sci::vector<double> > perl0(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > perl1(m+1, sci::vector<double>(n));

		// set Erlang pdf
		tmp = 0.0;
		for (int k=1; k<=m; k++) {
			tmp += tdat(k);
			for (int i=1; i<=n; i++) {
				perl0[k](i) = erlang_pdf(shape(i), rate(i), tmp);
				perl1[k](i) = tmp * perl0[k](i);
			}
		}

		etotal = 0.0;
		eb = 0.0;
		ew = 0.0;
		for (int k=1; k<=m; k++) {
			scale = sci::ddot(alpha, perl0[k]);
			sci::daxpy(wdat(k)/scale, perl0[k], eb);
			sci::daxpy(wdat(k)/scale, perl1[k], ew);
			llf += wdat(k) * log(scale);
			etotal += wdat(k);
		}
		eb *= alpha;
		ew *= alpha;
		return llf;
	}

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
		sci::vector<double>& ew) {

		int n = alpha.size;
		int m = tdat.size;

		sci::array< sci::vector<double> > perl0(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > perl1(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > cerl0(m+2, sci::vector<double>(n));
		sci::array< sci::vector<double> > cerl1(m+2, sci::vector<double>(n));

		// set Erlang pdf and cdf
		double tmp = 0.0;
		cerl0[0] = 0.0;
		cerl1[0] = 0.0;
		for (int k=1; k<=m; k++) {
			tmp += tdat(k);
			for (int i=1; i<=n; i++) {
				perl0[k](i) = erlang_pdf(shape(i), rate(i), tmp);
				perl1[k](i) = tmp * perl0[k](i);
				cerl0[k](i) = erlang_cdf(shape(i), rate(i), tmp);
				cerl1[k](i) = (shape(i) / rate(i)) * erlang_cdf(shape(i)+1, rate(i), tmp);
			}
		}
		cerl0[m+1] = 1.0;
		for (int i=1; i<=n; i++) {
			cerl1[m+1](i) = shape(i) / rate(i);
		}

		// estep
		sci::vector<double> tmpv0(n);
		sci::vector<double> tmpv1(n);
		double scale;
		double nn = 0.0;
		double uu = 0.0;
		double llf = 0.0;
		etotal = 0.0;
		eb = 0.0;
		ew = 0.0;
		for (int k=1; k<=m; k++) {
			if (gdat(k) >= 0 && tdat(k) != 0.0) {
				tmpv0 = cerl0[k] - cerl0[k-1];
				tmpv1 = cerl1[k] - cerl1[k-1];
				scale = sci::ddot(alpha, tmpv0);
				nn += gdat(k);
				uu += scale;
				sci::daxpy(gdat(k)/scale, tmpv0, eb);
				sci::daxpy(gdat(k)/scale, tmpv1, ew);
				llf += gdat(k) * log(scale) - lfact(gdat(k));
			}
			if (idat(k) == 1) {
				scale = sci::ddot(alpha, perl0[k]);
				nn += 1.0;
				sci::daxpy(1.0/scale, perl0[k], eb);
				sci::daxpy(1.0/scale, perl1[k], ew);
				llf += log(scale);
			}
		}
		if (gdatlast >= 0) {
			tmpv0 = cerl0[m+1] - cerl0[m];
			tmpv1 = cerl1[m+1] - cerl1[m];
			scale = sci::ddot(alpha, tmpv0);
			nn += gdatlast;
			uu += scale;
			sci::daxpy(gdatlast/scale, tmpv0, eb);
			sci::daxpy(gdatlast/scale, tmpv1, ew);
			llf += gdatlast * log(scale) - lfact(gdatlast);
		}

		for (int k=1; k<=m; k++) {
			if (gdat(k) == -1) {
				tmpv0 = cerl0[k] - cerl0[k-1];
				tmpv1 = cerl1[k] - cerl1[k-1];
				sci::daxpy(nn/uu, tmpv0, eb);
				sci::daxpy(nn/uu, tmpv1, ew);
			}
		}
		if (gdatlast == -1) {
			tmpv0 = cerl0[m+1] - cerl0[m];
			tmpv1 = cerl1[m+1] - cerl1[m];
			sci::daxpy(nn/uu, tmpv0, eb);
			sci::daxpy(nn/uu, tmpv1, ew);
		}
		llf += mylgamma(nn + 1.0) - nn * log(uu);

		etotal = nn / uu;
		eb *= alpha;
		ew *= alpha;
		return llf;
	}

	double phase_erlang_estep_group_trunc_poi(
		const sci::vector<double>& alpha,
		const sci::vector<double>& shape,
		const sci::vector<double>& rate,
		double omega,
		const sci::vector<double>& tdat,
		const sci::vector<int>& gdat,
		int gdatlast,
		const sci::vector<int>& idat,
		double& etotal,
		sci::vector<double>& eb,
		sci::vector<double>& ew) {

		int n = alpha.size;
		int m = tdat.size;

		sci::array< sci::vector<double> > perl0(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > perl1(m+1, sci::vector<double>(n));
		sci::array< sci::vector<double> > cerl0(m+2, sci::vector<double>(n));
		sci::array< sci::vector<double> > cerl1(m+2, sci::vector<double>(n));

		// set Erlang pdf and cdf
		double tmp = 0.0;
		cerl0[0] = 0.0;
		cerl1[0] = 0.0;
		for (int k=1; k<=m; k++) {
			tmp += tdat(k);
			for (int i=1; i<=n; i++) {
				perl0[k](i) = erlang_pdf(shape(i), rate(i), tmp);
				perl1[k](i) = tmp * perl0[k](i);
				cerl0[k](i) = erlang_cdf(shape(i), rate(i), tmp);
				cerl1[k](i) = (shape(i) / rate(i)) * erlang_cdf(shape(i)+1, rate(i), tmp);
			}
		}
		cerl0[m+1] = 1.0;
		cerl1[m+1] = shape / rate;

		// estep
		sci::vector<double> tmpv0(n);
		sci::vector<double> tmpv1(n);
		double scale;
		double nn = 0.0;
		double uu = 0.0;
		double llf = 0.0;
		etotal = 0.0;
		eb = 0.0;
		ew = 0.0;
		for (int k=1; k<=m; k++) {
			if (gdat(k) >= 0 && tdat(k) != 0.0) {
				tmpv0 = cerl0[k] - cerl0[k-1];
				tmpv1 = cerl1[k] - cerl1[k-1];
				scale = sci::ddot(alpha, tmpv0);
				nn += gdat(k);
				uu += scale;
				sci::daxpy(gdat(k)/scale, tmpv0, eb);
				sci::daxpy(gdat(k)/scale, tmpv1, ew);
				llf += gdat(k) * log(scale) - lfact(gdat(k));
			}
			if (idat(k) == 1) {
				scale = sci::ddot(alpha, perl0[k]);
				nn += 1.0;
				sci::daxpy(1.0/scale, perl0[k], eb);
				sci::daxpy(1.0/scale, perl1[k], ew);
				llf += log(scale);
			}
		}
		if (gdatlast >= 0) {
			tmpv0 = cerl0[m+1] - cerl0[m];
			tmpv1 = cerl1[m+1] - cerl1[m];
			scale = sci::ddot(alpha, tmpv0);
			nn += gdatlast;
			uu += scale;
			sci::daxpy(gdatlast/scale, tmpv0, eb);
			sci::daxpy(gdatlast/scale, tmpv1, ew);
			llf += gdatlast * log(scale) - lfact(gdatlast);
		}

		for (int k=1; k<=m; k++) {
			if (gdat(k) == -1) {
				tmpv0 = cerl0[k] - cerl0[k-1];
				tmpv1 = cerl1[k] - cerl1[k-1];
				sci::daxpy(omega, tmpv0, eb);
				sci::daxpy(omega, tmpv1, ew);
			}
		}
		if (gdatlast == -1) {
			tmpv0 = cerl0[m+1] - cerl0[m];
			tmpv1 = cerl1[m+1] - cerl1[m];
			sci::daxpy(omega, tmpv0, eb);
			sci::daxpy(omega, tmpv1, ew);
		}
		llf += nn * log(omega) - omega * uu;

		etotal = nn + omega * (1.0 - uu);
		eb *= alpha;
		ew *= alpha;
		return llf;
	}
}
