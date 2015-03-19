/*
 lapack.h
 wrapper for lapack
 */

#ifndef _WRAP_LAPACK_H
#define _WRAP_LAPACK_H

#ifdef __cplusplus
 extern "C" {
#endif
    
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

   int blas_dgesv(int n, int nrhs, double *A, int lda,
		  int *ipiv, double *B, int ldb);
   
   /*
     DPOTRF
     
     Description:
     Cholesky factorization of a real symmetric positive definete matrix A
     
     A = U^T * U if UPLO = 'U'
     A = L * L^T if UPLO = 'L'
     
     Arguments:
     uplo: 'U' or 'L' (see above)
     n: the number of columns (rows) of A
     A: a square matrix
     lda: LDA of A
     
     return value: interger of infomation
     info=0: success
     info<0: if INFO = -i, the i-th argument had an illegal value
     info>0: if INFO = i, the leading minor of order i is not
     positive definite, and the factorization could not be
     completed.
     
   */
   
   int blas_dpotrf(char uplo, int n, double *A, int lda);
   
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
   
   int blas_dggglm(int n, int m, int p, double *A, int lda,
		   double *B, int ldb, double *d, double *x,
		   double *y, double *work, int lwork);
   
#ifdef __cplusplus
}
#endif

#endif
