/*
 DGGGLM
 
 Description:
 
 min_x || B^-1 (d - A x) ||_2
 
 Solve the above least square method
 
 Arguments:
 n: the number of rows of A and B (N>=0)
 m: the number of columns of A (0<=M<=N)
 p: the number of columns of B (P >= N-M)
 A: matrix A (dimension LDA,M).
 lda: leading dimension of A (LDA >= max(1,N))
 B: matrix B (dimension LDB,P)
 ldb: leading dimension of B (LDB >= max(1,N))
 d: vector d. dimension N
 x: vector m. dimension M
 y: vector p. dimension P
 work
 lwork
 
 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK. LWORK >= max(1,N+M+P).
 *          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB,
 *          where NB is an upper bound for the optimal blocksizes for
 *          DGEQRF, SGERQF, DORMQR and SORMRQ.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued by XERBLA.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit.
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          = 1:  the upper triangular factor R associated with A in the
 *                generalized QR factorization of the pair (A, B) is
 *                singular, so that rank(A) < M; the least squares
 *                solution could not be computed.
 *          = 2:  the bottom (N-M) by (N-M) part of the upper trapezoidal
 *                factor T associated with B in the generalized QR
 *                factorization of the pair (A, B) is singular, so that
 *                rank( A B ) < N; the least squares solution could not
 *                be computed.
 */

#include <stdio.h>
#include <stdlib.h>

int mylapack_dggglm(int n, int m, int p, double *A, int lda,
		    double *B, int ldb, double *d, double *x,
		    double *y, double *work, int lwork) {
//  fprintf(stderr, "Please use f77lapack or MKL lapack\n");
//  exit(1);
	return -1;
}

