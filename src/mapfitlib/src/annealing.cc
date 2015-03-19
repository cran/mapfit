/*
  Description: annealing of PH/MAP parameters with temprature K

    Q = K * Q (diagonals)
  	Q = Q^K (except for diagonals)
    xi = xi^K
    alpha = alpha^K

       K (in): temprature
       alpha (inout): initial vector
       xi (inout): exit vector
       Q (inout): infinitesimal generator
       theta (out): normalized constant for Q and xi
       c (out): normalized constant for alpha
 */

#include <cmath>
#include <memory>

#include "phfit.hh"

namespace mapfit {

	void phase_annealing(double K,
		sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::matrix<double>& Q,
		double& theta, double& c) {
        switch (Q.type()) {
        case (DENSE):
            phase_annealing(K, alpha, xi, dynamic_cast<sci::dmatrix<double>&>(Q), theta, c);
            break;
        case (CSR):
            phase_annealing(K, alpha, xi, dynamic_cast<sci::csrmatrix<double>&>(Q), theta, c);
            break;
        case (CSC):
            phase_annealing(K, alpha, xi, dynamic_cast<sci::cscmatrix<double>&>(Q), theta, c);
            break;
        case (COO):
            phase_annealing(K, alpha, xi, dynamic_cast<sci::coomatrix<double>&>(Q), theta, c);
            break;
        default:
            throw;
        }
        return;
    }

	void phase_annealing(double K,
		sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::dmatrix<double>& Q,
		double& theta, double& c) {

		int n = Q.nrow;
		sci::vector<double> tmp(n);
		for (int j=1; j<=n; j++) {
			for (int i=1; i<=n; i++) {
				if (i == j) {
					Q(i,j) = K * Q(i,j);
				} else {
					Q(i,j) = pow(Q(i,j), K);
				}
				tmp(i) += Q(i,j);
			}
		}
		for (int i=1; i<n; i++) {
			alpha(i) = pow(alpha(i), K);
			xi(i) = pow(xi(i), K);
			tmp(i) += xi(i);
		}
		theta = sci::dmax(tmp);
		for (int i=1; i<=n; i++) {
			Q(i,i) -= theta;
		}
		c = dsum(alpha);
		dscal(1.0/c, alpha);
	}

	void phase_annealing(double K,
		sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::csrmatrix<double>& Q,
		double& theta, double& c) {

		int n = Q.nrow;
		sci::vector<double> tmp(n);
		sci::vector<int> diag(n);

		sci::array<int>::range rowptr = Q.get_rowptr();
		sci::array<int>::range colind = Q.get_colind();
		sci::array<double>::range value = Q.get_value();

		for (int i=1; i<=n; i++) {
			for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
				if (i == colind[z]) {
					value[z] = K * value[z];
					diag(i) = z;
				} else {
					value[z] = pow(value[z], K);
				}
				tmp(i) += value[z];
			}
		}
		for (int i=1; i<n; i++) {
			alpha(i) = pow(alpha(i), K);
			xi(i) = pow(xi(i), K);
			tmp(i) += xi(i);
		}
		theta = sci::dmax(tmp);
		for (int i=1; i<=n; i++) {
			value[diag(i)] -= theta;
		}
		c = dsum(alpha);
		dscal(1.0/c, alpha);
	}

	void phase_annealing(double K,
		sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::cscmatrix<double>& Q,
		double& theta, double& c) {

		int n = Q.nrow;
		sci::vector<double> tmp(n);
		sci::vector<int> diag(n);

		sci::array<int>::range colptr = Q.get_colptr();
		sci::array<int>::range rowind = Q.get_rowind();
		sci::array<double>::range value = Q.get_value();

		for (int j=1; j<=n; j++) {
			for (int z=colptr[j]; z<colptr[j+1]; z++) {
				if (j == rowind[z]) {
					value[z] = K * value[z];
					diag(j) = z;
				} else {
					value[z] = pow(value[z], K);
				}
				tmp(rowind[z]) += value[z];
			}
		}
		for (int i=1; i<n; i++) {
			alpha(i) = pow(alpha(i), K);
			xi(i) = pow(xi(i), K);
			tmp(i) += xi(i);
		}
		theta = sci::dmax(tmp);
		for (int i=1; i<=n; i++) {
			value[diag(i)] -= theta;
		}
		c = dsum(alpha);
		dscal(1.0/c, alpha);
	}

	void phase_annealing(double K,
		sci::vector<double>& alpha, sci::vector<double>& xi,
		sci::coomatrix<double>& Q,
		double& theta, double& c) {

		int n = Q.nrow;
		sci::vector<double> tmp(n);
		sci::vector<int> diag(n);

		sci::array<int>::range rowind = Q.get_rowind();
		sci::array<int>::range colind = Q.get_colind();
		sci::array<double>::range value = Q.get_value();

		for (size_t z=1; z<=Q.nnz; z++) {
			if (rowind[z] == colind[z]) {
				value[z] = K * value[z];
				diag(rowind[z]) = z;				
			} else {
				value[z] = pow(value[z], K);
			}
			tmp(rowind[z]) += value[z];
		}
		for (int i=1; i<n; i++) {
			alpha(i) = pow(alpha(i), K);
			xi(i) = pow(xi(i), K);
			tmp(i) += xi(i);
		}
		theta = sci::dmax(tmp);
		for (int i=1; i<=n; i++) {
			value[diag(i)] -= theta;
		}
		c = dsum(alpha);
		dscal(1.0/c, alpha);
	}
}
