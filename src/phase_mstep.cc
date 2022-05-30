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

#include "mexp.h"
#include "phfit.h"

namespace mapfit {

	void phase_mstep(
		const double& etotal,
		const sci::vector<double>& eb,
		const sci::vector<double>& ey,
		const sci::vector<double>& ez,
		const sci::matrix<double>& en,
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::matrix<double>& Q) {
	    switch (Q.type()) {
	    case (DENSE):
	        return phase_mstep(etotal, eb, ey, ez, dynamic_cast<const sci::dmatrix<double>&>(en),
	        	alpha, xi, dynamic_cast<sci::dmatrix<double>&>(Q));
	    case (CSR):
	        return phase_mstep(etotal, eb, ey, ez, dynamic_cast<const sci::csrmatrix<double>&>(en),
	        	alpha, xi, dynamic_cast<sci::csrmatrix<double>&>(Q));
	    case (CSC):
	        return phase_mstep(etotal, eb, ey, ez, dynamic_cast<const sci::cscmatrix<double>&>(en),
	        	alpha, xi, dynamic_cast<sci::cscmatrix<double>&>(Q));
	    case (COO):
	        return phase_mstep(etotal, eb, ey, ez, dynamic_cast<const sci::coomatrix<double>&>(en),
	        	alpha, xi, dynamic_cast<sci::coomatrix<double>&>(Q));
	    default:
	        throw;
	    }
	}

	void phase_mstep(
		const double& etotal,
		const sci::vector<double>& eb,
		const sci::vector<double>& ey,
		const sci::vector<double>& ez,
		const sci::dmatrix<double>& en,
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::dmatrix<double>& Q) {

		int n = Q.nrow;

		alpha = eb;
		dscal(1.0/etotal, alpha);
		xi = ey;
		xi /= ez;

		sci::vector<double> tmp(xi);
		for (int j=1; j<=n; j++) {
			for (int i=1; i<=n; i++) {
				if (i != j) {
					Q(i,j) = en(i,j) / ez(i);
					tmp(i) += Q(i,j);
				}
			}
		}

		for (int i=1; i<=n; i++) {
			Q(i,i) = -tmp(i);
		}
	}

	void phase_mstep(
		const double& etotal,
		const sci::vector<double>& eb,
		const sci::vector<double>& ey,
		const sci::vector<double>& ez,
		const sci::csrmatrix<double>& en,
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::csrmatrix<double>& Q) {

		int n = Q.nrow;

		alpha = eb;
		dscal(1.0/etotal, alpha);
		xi = ey;
		xi /= ez;

		sci::vector<int> idiag(n);
		sci::vector<double> tmp(xi);

	    sci::array<int>::range rowptr = Q.get_rowptr();
	    sci::array<int>::range colind = Q.get_colind();
	    sci::array<double>::range Q_value = Q.get_value();
	    sci::array<double>::const_range en_value = en.get_value();

		for (int i=1; i<=n; i++) {
			for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
				if (i != colind[z]) {
					Q_value[z] = en_value[z] / ez(i);
					tmp(i) += Q_value[z];
				} else {
					idiag(i) = z;
				}
			}
		}

		for (int i=1; i<=n; i++) {
			Q_value[idiag(i)] = -tmp(i);
		}
	}

	void phase_mstep(
		const double& etotal,
		const sci::vector<double>& eb,
		const sci::vector<double>& ey,
		const sci::vector<double>& ez,
		const sci::cscmatrix<double>& en,
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::cscmatrix<double>& Q) {

		int n = Q.nrow;

		alpha = eb;
		dscal(1.0/etotal, alpha);
		xi = ey;
		xi /= ez;

		sci::vector<int> idiag(n);
		sci::vector<double> tmp(xi);

	    sci::array<int>::range colptr = Q.get_colptr();
	    sci::array<int>::range rowind = Q.get_rowind();
	    sci::array<double>::range Q_value = Q.get_value();
	    sci::array<double>::const_range en_value = en.get_value();

		for (int j=1; j<=n; j++) {
			for (int z=colptr[j]; z<colptr[j+1]; z++) {
				if (j != rowind[z]) {
					Q_value[z] = en_value[z] / ez(rowind[z]);
					tmp(rowind[z]) += Q_value[z];
				} else {
					idiag(rowind[z]) = z;
				}
			}
		}

		for (int i=1; i<=n; i++) {
			Q_value[idiag(i)] = -tmp(i);
		}
	}

	void phase_mstep(
		const double& etotal,
		const sci::vector<double>& eb,
		const sci::vector<double>& ey,
		const sci::vector<double>& ez,
		const sci::coomatrix<double>& en,
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::coomatrix<double>& Q) {

		int n = Q.nrow;

		alpha = eb;
		dscal(1.0/etotal, alpha);
		xi = ey;
		xi /= ez;

		sci::vector<int> idiag(n);
		sci::vector<double> tmp(xi);

	    sci::array<int>::range rowind = Q.get_rowind();
	    sci::array<int>::range colind = Q.get_colind();
	    sci::array<double>::range Q_value = Q.get_value();
	    sci::array<double>::const_range en_value = en.get_value();

		for (size_t z=1; z<=Q.nnz; z++) {
			if (rowind[z] != colind[z]) {
				Q_value[z] = en_value[z] / ez(rowind[z]);
				tmp(rowind[z]) += Q_value[z];
			} else {
				idiag(rowind[z]) = z;
			}
		}

		for (int i=1; i<=n; i++) {
			Q_value[idiag(i)] = -tmp(i);
		}
	}

	///// for cf1

	void cf1_swap(int i, int j, sci::vector<double>& alpha,
		const sci::vector<double*>& rate,
		const sci::vector<double*>& diag) {
		double tmp;
		double w = *rate(j) / *rate(i);
		alpha(i) += (1.0 - w) * alpha(j);
		alpha(j) *= w;
		tmp = *rate(j);
		*rate(j) = *rate(i);
		*rate(i) = tmp;
		*diag(i) = -*rate(i);
		*diag(j) = -*rate(j);
	}

	void cf1_sort(sci::vector<double>& alpha,
		const sci::vector<double*>& rate,
		const sci::vector<double*>& diag) {
		int n = alpha.size;
		for (int i=1; i<=n-1; i++) {
			if (*rate(i) > *rate(i+1)) {
				cf1_swap(i, i+1, alpha, rate, diag);
				for (int j=i; j>=2; j--) {
					if (*rate(j-1) < *rate(j)) {
						break;
					}
					cf1_swap(j-1, j, alpha, rate, diag);
				}
			}
		}
	}

	void phase_bidiag_to_cf1(
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::matrix<double>& Q) {
	    switch (Q.type()) {
	    case (DENSE):
	        return phase_bidiag_to_cf1(alpha, xi, dynamic_cast<sci::dmatrix<double>&>(Q));
	    case (CSR):
	        return phase_bidiag_to_cf1(alpha, xi, dynamic_cast<sci::csrmatrix<double>&>(Q));
	    case (CSC):
	        return phase_bidiag_to_cf1(alpha, xi, dynamic_cast<sci::cscmatrix<double>&>(Q));
	    case (COO):
	        return phase_bidiag_to_cf1(alpha, xi, dynamic_cast<sci::coomatrix<double>&>(Q));
	    default:
	        throw;
	    }
	}

	void phase_bidiag_to_cf1(
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::dmatrix<double>& Q) {
		int n = alpha.size;
		sci::vector<double*> rate(n);
		sci::vector<double*> diag(n);
		for (int i=1; i<=n-1; i++) {
			diag(i) = &Q(i,i);
			rate(i) = &Q(i,i+1);
		}
		diag(n) = &Q(n,n);
		rate(n) = &xi(n);
		cf1_sort(alpha, rate, diag);
	}

	void phase_bidiag_to_cf1(
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::csrmatrix<double>& Q) {
		int n = alpha.size;
		sci::vector<double*> rate(n);
		sci::vector<double*> diag(n);

	    sci::array<int>::range rowptr = Q.get_rowptr();
	    sci::array<int>::range colind = Q.get_colind();
	    sci::array<double>::range value = Q.get_value();

		for (int i=1; i<=n; i++) {
			for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
				if (i == colind[z]) {
					diag(i) = &value[z];
				}
				if (i+1 == colind[z]) {
					rate(i) = &value[z];
				}
			}
		}
		rate(n) = &xi(n);
		cf1_sort(alpha, rate, diag);
	}

	void phase_bidiag_to_cf1(
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::cscmatrix<double>& Q) {
		int n = alpha.size;
		sci::vector<double*> rate(n);
		sci::vector<double*> diag(n);

	    sci::array<int>::range colptr = Q.get_colptr();
	    sci::array<int>::range rowind = Q.get_rowind();
	    sci::array<double>::range value = Q.get_value();

		for (int j=1; j<=n; j++) {
			for (int z=colptr[j]; z<colptr[j+1]; z++) {
				if (j == rowind[z]) {
					diag(j) = &value[z];
				}
				if (j == rowind[z]+1) {
					rate(rowind[z]) = &value[z];
				}
			}
		}
		rate(n) = &xi(n);
		cf1_sort(alpha, rate, diag);
	}

	void phase_bidiag_to_cf1(
		sci::vector<double>& alpha,
		sci::vector<double>& xi,
		sci::coomatrix<double>& Q) {
		int n = alpha.size;
		sci::vector<double*> rate(n);
		sci::vector<double*> diag(n);

	    sci::array<int>::range rowind = Q.get_rowind();
	    sci::array<int>::range colind = Q.get_colind();
	    sci::array<double>::range value = Q.get_value();

		for (size_t z=1; z<=Q.nnz; z++) {
			if (rowind[z] == colind[z]) {
				diag(rowind[z]) = &value[z];
			}
			if (rowind[z]+1 == colind[z]) {
				rate(rowind[z]) = &value[z];
			}
		}
		rate(n) = &xi(n);
		cf1_sort(alpha, rate, diag);
	}
}
