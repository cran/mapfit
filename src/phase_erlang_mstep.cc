/*
  Description: estep for Erlang-PH with weighted time and group/truncated data

   alpha   (inout): initial vector
   shape   (inout): shape parameter vector
   rate    (inout): rate parameter vector
   etotal  (in): expected # of arrivals
   eb      (in): expected # of starts
   ew      (in): expected sojourn time?

 */

#include <cmath>
#include <memory>

#include "mexp.h"
#include "phfit.h"

namespace mapfit {

	void phase_erlang_mstep(
		const double& etotal,
		const sci::vector<double>& eb,
		const sci::vector<double>& ew,
		sci::vector<double>& alpha,
		sci::vector<int>& shape,
		sci::vector<double>& rate) {

		int n = alpha.size;
		alpha = eb;
		dscal(1.0/etotal, alpha);
		for (int i=1; i<=n; i++) {
			rate(i) = shape(i) * eb(i) / ew(i);
		}
	}

}
