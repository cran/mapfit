

#include <cmath>
#include <memory>

#include "mexp.h"
#include "mapfit.h"

namespace mapfit {

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::matrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::matrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate) {
	    switch (P.type()) {
	    case (DENSE):
	        return hmm_erlang_mstep(eb, dynamic_cast<const sci::dmatrix<double>&>(en), ew0, ew1, 
	        	alpha, dynamic_cast<sci::dmatrix<double>&>(P), shape, rate);
	    case (CSR):
	        return hmm_erlang_mstep(eb, dynamic_cast<const sci::csrmatrix<double>&>(en), ew0, ew1, 
	        	alpha, dynamic_cast<sci::csrmatrix<double>&>(P), shape, rate);
	    case (CSC):
	        return hmm_erlang_mstep(eb, dynamic_cast<const sci::cscmatrix<double>&>(en), ew0, ew1, 
	        	alpha, dynamic_cast<sci::cscmatrix<double>&>(P), shape, rate);
	    case (COO):
	        return hmm_erlang_mstep(eb, dynamic_cast<const sci::coomatrix<double>&>(en), ew0, ew1, 
	        	alpha, dynamic_cast<sci::coomatrix<double>&>(P), shape, rate);
	    default:
	        throw;
	    }
	}

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::dmatrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::dmatrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate) {

		int n = alpha.size;

		sci::vector<double> tmp(n);
		for (int j=1; j<=n; j++) {
			for (int i=1; i<=n; i++) {
				tmp(i) += en(i,j);
			}
		}

		for (int j=1; j<=n; j++) {
			for (int i=1; i<=n; i++) {
				P(i,j) = en(i,j) / tmp(i);
			}
		}

		alpha = eb;

		for (int i=1; i<=n; i++) {
			rate(i) = ew0(i) * shape(i) / ew1(i);
		}
	}

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::csrmatrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::csrmatrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate) {

		int n = alpha.size;

		sci::vector<double> tmp(n);

	    sci::array<int>::range rowptr = P.get_rowptr();
	    sci::array<int>::range colind = P.get_colind();
	    sci::array<double>::range value = P.get_value();
	    sci::array<double>::const_range en_value = en.get_value();

		for (int i=1; i<=n; i++) {
			for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
				tmp(i) += en_value[z];
			}
		}

		for (int i=1; i<=n; i++) {
			for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
				value[z] = en_value[z] / tmp(i);
			}
		}

		alpha = eb;

		for (int i=1; i<=n; i++) {
			rate(i) = ew0(i) * shape(i) / ew1(i);
		}
	}

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::cscmatrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::cscmatrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate) {

		int n = alpha.size;

		sci::vector<double> tmp(n);

	    sci::array<int>::range colptr = P.get_colptr();
	    sci::array<int>::range rowind = P.get_rowind();
	    sci::array<double>::range value = P.get_value();
	    sci::array<double>::const_range en_value = en.get_value();

		for (int j=1; j<=n; j++) {
			for (int z=colptr[j]; z<colptr[j+1]; z++) {
				tmp(rowind[z]) += en_value[z];
			}
		}

		for (int j=1; j<=n; j++) {
			for (int z=colptr[j]; z<colptr[j+1]; z++) {
				value[z] = en_value[z] / tmp(rowind[z]);
			}
		}

		alpha = eb;

		for (int i=1; i<=n; i++) {
			rate(i) = ew0(i) * shape(i) / ew1(i);
		}
	}

	void hmm_erlang_mstep(
		const sci::vector<double>& eb,
		const sci::coomatrix<double>& en,
		const sci::vector<double>& ew0,
		const sci::vector<double>& ew1,
		sci::vector<double>& alpha,
		sci::coomatrix<double>& P,
		sci::vector<int>& shape,
		sci::vector<double>& rate) {

		int n = alpha.size;

		sci::vector<double> tmp(n);

	    sci::array<int>::range rowind = P.get_rowind();
	    sci::array<int>::range colind = P.get_colind();
	    sci::array<double>::range value = P.get_value();
	    sci::array<double>::const_range en_value = en.get_value();

		for (size_t z=1; z<=P.nnz; z++) {
			tmp(rowind[z]) += en_value[z];
		}

		for (size_t z=1; z<=P.nnz; z++) {
			value[z] = en_value[z] / tmp(rowind[z]);
		}

		alpha = eb;

		for (int i=1; i<=n; i++) {
			rate(i) = ew0(i) * shape(i) / ew1(i);
		}
	}

}
