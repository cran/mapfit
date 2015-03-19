/* wrapper for lapack */

#ifdef MYLAPACK
#include "dgesv.h"
#include "dpotrf.h"
#include "dggglm.h"
#define __DGESV__ mylapack_dgesv
#define __DGGGLM__ mylapack_dggglm
#define __DPOTRF__ mylapack_dpotrf
#endif

#ifdef F77LAPACK
/*----- prototype -----*/
int dgesv_(int* n, int* nrhs, double* A, int* lda,
	   int* ipiv, double* B, int* ldb, int* info);
int dggglm_(int* n, int* m, int* p,
            double* A, int* lda,
            double* B, int* ldb,
            double* d, double* x, double* y,
            double* work, int* lwork, int* info);
int dpotrf_(char* uplo, int* n, double* A, int* lda, int* info);
/*---------------------*/

#define __DGESV__ dgesv_
#define __DGGGLM__ dggglm_
#define __DPOTRF__ dpotrf_
#endif

#ifdef MKLLAPACK
/*----- prototype -----*/
int dgesv_(int* n, int* nrhs, double* A, int* lda,
	   int* ipiv, double* B, int* ldb, int* info);
int dggglm_(int* n, int* m, int* p,
            double* A, int* lda,
            double* B, int* ldb,
            double* d, double* x, double* y,
            double* work, int* lwork, int* info);
int dpotrf_(char* uplo, int* n, double* A, int* lda, int* info);
/*---------------------*/

#define __DGESV__ dgesv_
#define __DGGGLM__ dggglm_
#define __DPOTRF__ dpotrf_
#endif

/*----- wrappers -----*/

int blas_dgesv(int n, int nrhs, double *A, int lda,
	       int *ipiv, double *B, int ldb) {
  int info;
  __DGESV__(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
  return info;
}

int blas_dpotrf(char uplo, int n, double *A, int lda) {
  int info;
  __DPOTRF__(&uplo, &n, A, &lda, &info);
  return info;
}

int blas_dggglm(int n, int m, int p, double *A, int lda,
		double *B, int ldb, double *d, double *x,
		double *y, double *work, int lwork) {
  int info;
  __DGGGLM__(&n, &m, &p, A, &lda, B, &ldb, d, x, y, work, &lwork, &info);
  return info;
}

