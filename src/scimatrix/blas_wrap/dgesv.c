/*
 DGESV
 
 Description:
 
 A X = B
 
 The above linear equation is solved by LU decomposition
 
 Arguments:
 n: the number of rows of A
 nrhs: the number of columns of B
 A: matrix A. The result is the matrix that is decomposed by LU
 lda: leading dimension of A
 ipiv: substitued indexes
 B: matrix B. The result is X.
 ldb: leading dimension of B
 
*/

#include <stdio.h>
#include <stdlib.h>

int mylapack_dgesv(int n, int nrhs, double *A, int lda,
		   int *ipiv, double *B, int ldb) {
//  fprintf(stderr, "Please use f77lapack or MKL lapack\n");
//  exit(1);
	return -1;
}




