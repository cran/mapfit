/*
  myblas

*/

#ifndef _MYBLAS_H
#define _MYBLAS_H

#ifdef __cplusplus
extern "C" {
#endif

void myblas_dcopy(const int *n, const double *x,
		const int *incx, double *y, const int *incy);

void myblas_dscal(const int *n,
		const double *alpha,
		double *x,
		const int *incx);

void myblas_daxpy(const int *n,
		const double *alpha,
		const double *x,
		const int *incx,
		double *y,
		const int *incy);

double myblas_ddot(const int *n,
		 const double *x,
		 const int *incx,
		 const double *y,
		 const int *incy);

double myblas_dnrm2(const int *n,
		  const double *x,
		  const int *incx);

double myblas_dasum(const int *n,
		  const double *x,
		  const int *incx);

// level 2

void myblas_dgemv(const char *trans,
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

void myblas_dger(const int *m,
	       const int *n,
	       const double *alpha,
	       const double *x,
	       const int *incx,
	       const double *y,
	       const int *incy,
	       double *A,
	       const int *lda);

// level 3

void myblas_dgemm(const char *transa,
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

#ifdef __cplusplus
}
#endif
#endif
