
#include <cmath>
#include <memory>

#include "mexp.h"
#include "mapfit.h"

namespace mapfit {

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::matrix<double>& en0,
		const sci::matrix<double>& en1,
		sci::vector<double>& alpha,
		sci::matrix<double>& D0,
		sci::matrix<double>& D1) {
	    switch (D0.type()) {
	    case (DENSE):
	        return map_mstep(eb, ez,
	        	dynamic_cast<const sci::dmatrix<double>&>(en0),
	        	dynamic_cast<const sci::dmatrix<double>&>(en1),
	        	alpha,
	        	dynamic_cast<sci::dmatrix<double>&>(D0),
	        	dynamic_cast<sci::dmatrix<double>&>(D1));
	    case (CSR):
	        return map_mstep(eb, ez,
	        	dynamic_cast<const sci::csrmatrix<double>&>(en0),
	        	dynamic_cast<const sci::csrmatrix<double>&>(en1),
	        	alpha,
	        	dynamic_cast<sci::csrmatrix<double>&>(D0),
	        	dynamic_cast<sci::csrmatrix<double>&>(D1));
	    case (CSC):
	        return map_mstep(eb, ez,
	        	dynamic_cast<const sci::cscmatrix<double>&>(en0),
	        	dynamic_cast<const sci::cscmatrix<double>&>(en1),
	        	alpha,
	        	dynamic_cast<sci::cscmatrix<double>&>(D0),
	        	dynamic_cast<sci::cscmatrix<double>&>(D1));
	    case (COO):
	        return map_mstep(eb, ez,
	        	dynamic_cast<const sci::coomatrix<double>&>(en0),
	        	dynamic_cast<const sci::coomatrix<double>&>(en1),
	        	alpha,
	        	dynamic_cast<sci::coomatrix<double>&>(D0),
	        	dynamic_cast<sci::coomatrix<double>&>(D1));
	    default:
	        throw;
	    }
	}

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::dmatrix<double>& en0,
		const sci::dmatrix<double>& en1,
		sci::vector<double>& alpha,
		sci::dmatrix<double>& D0,
		sci::dmatrix<double>& D1) {

		int n = alpha.size;

		alpha = eb;

		sci::vector<double> tmp(n);
		for (int j=1; j<=n; j++) {
			for (int i=1; i<=n; i++) {
				if (i != j) {
					D0(i,j) = en0(i,j) / ez(i);
					tmp(i) += D0(i,j);
				}
				D1(i,j) = en1(i,j) / ez(i);
				tmp(i) += D1(i,j);
			}
		}

		for (int i=1; i<=n; i++) {
			D0(i,i) = -tmp(i);
		}
	}

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::csrmatrix<double>& en0,
		const sci::csrmatrix<double>& en1,
		sci::vector<double>& alpha,
		sci::csrmatrix<double>& D0,
		sci::csrmatrix<double>& D1) {

		int n = alpha.size;

		alpha = eb;

		sci::vector<int> idiag(n);
		sci::vector<double> tmp(n);

	    sci::array<int>::range rowptr0 = D0.get_rowptr();
	    sci::array<int>::range colind0 = D0.get_colind();
	    sci::array<int>::range rowptr1 = D1.get_rowptr();
	    sci::array<int>::range colind1 = D1.get_colind();
	    sci::array<double>::range D0_value = D0.get_value();
	    sci::array<double>::range D1_value = D1.get_value();
	    sci::array<double>::const_range en0_value = en0.get_value();
	    sci::array<double>::const_range en1_value = en1.get_value();

		for (int i=1; i<=n; i++) {
			for (int z=rowptr0[i]; z<rowptr0[i+1]; z++) {
				if (i != colind0[z]) {
					D0_value[z] = en0_value[z] / ez(i);
					tmp(i) += D0_value[z];
				} else {
					idiag(i) = z;
				}
			}
			for (int z=rowptr1[i]; z<rowptr1[i+1]; z++) {
				D1_value[z] = en1_value[z] / ez(i);
				tmp(i) += D1_value[z];
			}
		}

		for (int i=1; i<=n; i++) {
			D0_value[idiag(i)] = -tmp(i);
		}
	}

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::cscmatrix<double>& en0,
		const sci::cscmatrix<double>& en1,
		sci::vector<double>& alpha,
		sci::cscmatrix<double>& D0,
		sci::cscmatrix<double>& D1) {

		int n = alpha.size;

		alpha = eb;

		sci::vector<int> idiag(n);
		sci::vector<double> tmp(n);

	    sci::array<int>::range colptr0 = D0.get_colptr();
	    sci::array<int>::range rowind0 = D0.get_rowind();
	    sci::array<int>::range colptr1 = D1.get_colptr();
	    sci::array<int>::range rowind1 = D1.get_rowind();
	    sci::array<double>::range D0_value = D0.get_value();
	    sci::array<double>::range D1_value = D1.get_value();
	    sci::array<double>::const_range en0_value = en0.get_value();
	    sci::array<double>::const_range en1_value = en1.get_value();

		for (int j=1; j<=n; j++) {
			for (int z=colptr0[j]; z<colptr0[j+1]; z++) {
				if (j != rowind0[z]) {
					D0_value[z] = en0_value[z] / ez(rowind0[z]);
					tmp(rowind0[z]) += D0_value[z];
				} else {
					idiag(rowind0[z]) = z;
				}
			}
			for (int z=colptr1[j]; z<colptr1[j+1]; z++) {
				D1_value[z] = en1_value[z] / ez(rowind1[z]);
				tmp(rowind1[z]) += D1_value[z];

			}
		}

		for (int i=1; i<=n; i++) {
 			D0_value[idiag(i)] = -tmp(i);
		}
	}

	void map_mstep(
		const sci::vector<double>& eb,
		const sci::vector<double>& ez,
		const sci::coomatrix<double>& en0,
		const sci::coomatrix<double>& en1,
		sci::vector<double>& alpha,
		sci::coomatrix<double>& D0,
		sci::coomatrix<double>& D1) {

		int n = alpha.size;

		alpha = eb;

		sci::vector<int> idiag(n);
		sci::vector<double> tmp(n);

	    sci::array<int>::range rowind0 = D0.get_rowind();
	    sci::array<int>::range colind0 = D0.get_colind();
	    sci::array<int>::range rowind1 = D1.get_rowind();
	    sci::array<int>::range colind1 = D1.get_colind();
	    sci::array<double>::range D0_value = D0.get_value();
	    sci::array<double>::range D1_value = D1.get_value();
	    sci::array<double>::const_range en0_value = en0.get_value();
	    sci::array<double>::const_range en1_value = en1.get_value();

	    for (size_t z=1; z<=en0.nnz; z++) {
	    	if (rowind0[z] != colind0[z]) {
	    		D0_value[z] = en0_value[z] / ez(rowind0[z]);
	    		tmp(rowind0[z]) += D0_value[z];
	    	} else {
	    		idiag(rowind0[z]) = z;
	    	}
	    }

	    for (size_t z=1; z<=en1.nnz; z++) {
    		D1_value[z] = en1_value[z] / ez(rowind1[z]);
    		tmp(rowind1[z]) += D1_value[z];
	    }

		for (int i=1; i<=n; i++) {
			D0_value[idiag(i)] = -tmp(i);
		}
	}
}
