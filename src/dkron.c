/* kronecker routines */

#include "blas.h"

/*
  !
  ! Description:
  !   This routine computes the following kronecker product:
  ! 
  !     y = alpha (I_{n1} * tras(A) * I_{n2}) x + beta y
  !
  !   where A is an m-by-n square matrix,
  !         I_{n} is an n-by-n identity matrix, and;
  !         * is kronecker product.
  !
  ! Parameters:
  !   m, n: rows and columns of matrix A
  !   n1: size of identity matrix located on left
  !   n2: size of identity matrix located on right
  !
 */

void blas_dkronmv(
	char trans,
	int m, int n, int n1, int n2,
	double alpha,
	const double *A, int lda,
	const double *x,
	double beta, double *y) {

	int i, j;
	if (trans == 'T' || trans == 't') {
		blas_dscal(n*n1*n2, beta, y, 1);
		for (i=0, j=0; i<m*n1*n2; i+=m*n2, j+=n*n2) {
			// blas_dgemm('N', 'N', n2, n, m,
			// 	alpha, &x[i], n2,
			// 	A, lda, 1.0, &y[j], n2);
			blas_dgemm('T', 'T', n, n2, m,
				alpha, A, lda, &x[i], n2,
				1.0, &y[j], n2);
		}
	} else if (trans == 'N' || trans == 'n') {
		blas_dscal(m*n1*n2, beta, y, 1);
		for (i=0, j=0; i<n*n1*n2; i+=n*n2, j+=m*n2) {
			// blas_dgemm('N', 'T', n2, m, n,
			// 	alpha, &x[i], n2,
			// 	A, lda, 1.0, &y[j], n2);
			blas_dgemm('T', 'N', m, n2, n,
				alpha, A, lda, &x[i], n2,
				1.0, &y[j], n2);
		}
	}
}

/*
  !
  ! Description:
  !   This routine computes the following aggregation on kronecker product:
  ! 
  !     B = B + alpha x^T y
  !       and element-wise aggregation of the m-by-n matrix A; I_{n1} * A * I_{n2}
  !
  !   where A is m-by-n matrix
  !         I_{n} is an n-by-n identity matrix, and;
  !         * is kronecker product.
  !
  ! Parameters:
  !   m, n: row and column sizes of square matrix A
  !   n1: size of identity matrix located on left
  !   n2: size of identity matrix located on right
  !
 */

void blas_dkronr(int m, int n, int n1, int n2,
		double alpha, double *x, double *y,
		double *A, int lda) {
	int i;
	for (i=0; i<n*n1*n2; i+=n*n2) {
		blas_dgemm('T', 'N', m, n, n2,
			alpha, &x[i], n2,
			&y[i], n2, 1.0, A, lda);
	}
}

