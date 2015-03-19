/*
 blas.c
 wrapper for blas
 */

#include <stdlib.h>
#include <math.h>

#ifdef MYBLAS
#include "myblas.h"
#define __DCOPY__ myblas_dcopy
#define __DSCAL__ myblas_dscal
#define __DAXPY__ myblas_daxpy
#define __DDOT__ myblas_ddot
#define __DNRM2__ myblas_dnrm2
#define __DASUM__ myblas_dasum
#define __DGEMV__ myblas_dgemv
#define __DGER__ myblas_dger
#define __DGEMM__ myblas_dgemm
#endif

#ifdef F77BLAS
/*----- prototype -------------*/
void dcopy_(const int *n, const double *x,
	      const int *incx, double *y, const int *incy);
void dscal_(const int *n,
	      const double *alpha,
	      double *x,
	      const int *incx);
void daxpy_(const int *n,
	      const double *alpha,
	      const double *x,
	      const int *incx,
	      double *y,
	      const int *incy);
double ddot_(const int *n,
	       const double *x,
	       const int *incx,
	       const double *y,
	       const int *incy);
double dnrm2_(const int *n,
		const double *x,
		const int *incx);
double dasum_(const int *n,
		const double *x,
		const int *incx);
void dgemv_(const char *trans,
	      const int *m,
	      const int *n,
	      const double *alpha,
	      const double *A,
	      const int *lda,
	      const double *x,
	      const int *incx,
	      const double *beta,
	      double *y,
	      const int *incy);
void dger_(const int *m,
	     const int *n,
	     const double *alpha,
	     const double *x,
	     const int *incx,
	     const double *y,
	     const int *incy,
	     double *A,
	     const int *lda);
void dgemm_(const char *transa,
	      const char *transb,
	      const int *m,
	      const int *n,
	      const int *k,
	      const double *alpha,
	      const double *A,
	      const int *lda,
	      const double *B,
	      const int *ldb,
	      const double *beta,
	      double *C,
	      const int *ldc);
/*------------------*/

#define __DCOPY__ dcopy_
#define __DSCAL__ dscal_
#define __DAXPY__ daxpy_
#define __DDOT__ ddot_
#define __DNRM2__ dnrm2_
#define __DASUM__ dasum_
#define __DGEMV__ dgemv_
#define __DGER__ dger_
#define __DGEMM__ dgemm_
#endif

#ifdef MKLBLAS
/*----- prototype -------------*/
void dcopy_(const int *n, const double *x,
	      const int *incx, double *y, const int *incy);
void dscal_(const int *n,
	      const double *alpha,
	      double *x,
	      const int *incx);
void daxpy_(const int *n,
	      const double *alpha,
	      const double *x,
	      const int *incx,
	      double *y,
	      const int *incy);
double ddot_(const int *n,
	       const double *x,
	       const int *incx,
	       const double *y,
	       const int *incy);
double dnrm2_(const int *n,
		const double *x,
		const int *incx);
double dasum_(const int *n,
		const double *x,
		const int *incx);
void dgemv_(const char *trans,
	      const int *m,
	      const int *n,
	      const double *alpha,
	      const double *A,
	      const int *lda,
	      const double *x,
	      const int *incx,
	      const double *beta,
	      double *y,
	      const int *incy);
void dger_(const int *m,
	     const int *n,
	     const double *alpha,
	     const double *x,
	     const int *incx,
	     const double *y,
	     const int *incy,
	     double *A,
	     const int *lda);
void dgemm_(const char *transa,
	      const char *transb,
	      const int *m,
	      const int *n,
	      const int *k,
	      const double *alpha,
	      const double *A,
	      const int *lda,
	      const double *B,
	      const int *ldb,
	      const double *beta,
	      double *C,
	      const int *ldc);
/*------------------*/

#define __DCOPY__ dcopy_
#define __DSCAL__ dscal_
#define __DAXPY__ daxpy_
#define __DDOT__ ddot_
#define __DNRM2__ dnrm2_
#define __DASUM__ dasum_
#define __DGEMV__ dgemv_
#define __DGER__ dger_
#define __DGEMM__ dgemm_
#endif

// level 1

void blas_dcopy(int n, const double *x, int incx, double *y, int incy) {
  __DCOPY__(&n, x, &incx, y, &incy);
}

void blas_dfill(int n, double v, double *x, int incx) {
  int i;
  for (i=0; i<n; i++, x+=incx) {
    *x = v;
  }
}

void blas_dscal(int n, double alpha, double *x, int incx) {
  __DSCAL__(&n, &alpha, x, &incx);
}

void blas_daxpy(int n, double alpha, const double *x, int incx,
		double *y, int incy) {
  __DAXPY__(&n, &alpha, x, &incx, y, &incy);
}

double blas_ddot(int n, const double *x, int incx, const double *y, int incy) {
  double tmp;
  tmp = __DDOT__(&n, x, &incx, y, &incy);
  return tmp;
}

double blas_dnrm2(int n, const double *x, int incx) {
  double tmp;
  tmp = __DNRM2__(&n, x, &incx);
  return tmp;
}

double blas_dasum(int n, const double *x, int incx) {
  double tmp;
  tmp = __DASUM__(&n, x, &incx);
  return tmp;
}

double blas_dsum(int n, const double *x, int incx) {
  int i;
  double sum = 0.0;
  for (i=0; i<n; i++, x+=incx) {
    sum += *x;
  }
  return sum;
}

int blas_idmax(int n, const double *x, int incx) {
  int i, id;
  double max = *x;
  id = 0;
  for (i=0; i<n; i++, x+=incx) {
    if (*x > max) {
      max = *x;
      id = i;
    }
  }
  return id;
}

int blas_idamax(int n, const double *x, int incx) {
  int i, id;
  double tmp;
  double max = fabs(*x);
  id = 0;
  for (i=0; i<n; i++, x+=incx) {
    tmp = fabs(*x);
    if (tmp > max) {
      max = tmp;
      id = i;
    }
  }
  return id;
}


// level 2

void blas_dgemv(char trans, int m, int n,
		double alpha, const double *A, int lda,
		const double *x, int incx,
		double beta, double *y, int incy) {
  __DGEMV__(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

void blas_dger(int m, int n,
	       double alpha, const double *x, int incx,
	       const double *y, int incy,
	       double *A, int lda) {
  __DGER__(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}

// level 3

void blas_dgemm(char transa, char transb,
		int m, int n, int k,
		double alpha, const double *A, int lda,
		const double *B, int ldb,
		double beta, double *C, int ldc) {
  __DGEMM__(&transa, &transb, &m, &n, &k, &alpha, A, &lda,
	    B, &ldb, &beta, C, &ldc);
}

