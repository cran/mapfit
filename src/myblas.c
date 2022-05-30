/*
 blas.c
 wrapper for blas
 */

#include <stdlib.h>
#include <math.h>

// level 1

void myblas_dcopy(const int *n, const double *x,
		const int *incx, double *y, const int *incy) {
  int i;
  for (i=0; i<*n; i++, x+=*incx, y+=*incy) {
    *y = *x;
  }
}

void myblas_dscal(const int *n,
		const double *alpha,
		double *x,
		const int *incx) {
  int i;
  for (i=0; i<*n; i++, x+=*incx) {
    *x *= *alpha;
  }
}

void myblas_daxpy(const int *n,
		const double *alpha,
		const double *x,
		const int *incx,
		double *y,
		const int *incy) {
  int i;
  for (i=0; i<*n; i++, x+=*incx, y+=*incy) {
    *y += *alpha * *x;
  }
}

double myblas_ddot(const int *n,
		 const double *x,
		 const int *incx,
		 const double *y,
		 const int *incy) {
  int i;
  double sum = 0.0;
  for (i=0; i<*n; i++, x+=*incx, y+=*incy) {
    sum += *x * *y;
  }
  return sum;
}

double myblas_dnrm2(const int *n,
		  const double *x,
		  const int *incx) {
  int i;
  double sum = 0.0;
  for (i=0; i<*n; i++, x+=*incx) {
    sum += *x * *x;
  }
  return sqrt(sum);
}

double myblas_dasum(const int *n,
		  const double *x,
		  const int *incx) {
  int i;
  double sum = 0.0;
  for (i=0; i<*n; i++, x+=*incx) {
    sum += fabs(*x);
  }
  return sum;
}

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
		const int *incy) {
  int i, j;
  double *py;
  if (*trans == 'N' || *trans == 'n') {
    myblas_dscal(m, beta, y, incy);
    for (i=0; i<*n; i++, x+=*incx) {
      for (j=0, py=y; j<*m; j++, py+=*incy) {
	*py += *alpha * A[j + i*(*lda)] * *x;
      }
    }
  } else if (*trans == 'T' || *trans == 't') {
    myblas_dscal(n, beta, y, incy);
    for (i=0; i<*m; i++, x+=*incx) {
      for (j=0, py=y; j<*n; j++, py+=*incy) {
	*py += *alpha * A[i + j*(*lda)] * *x;
      }
    }
  }
}

void myblas_dger(const int *m,
	       const int *n,
	       const double *alpha,
	       const double *x,
	       const int *incx,
	       const double *y,
	       const int *incy,
	       double *A,
	       const int *lda) {
  int i, j;
  const double *px, *py;
  for (i=0, px=x; i<*m; i++, px+=*incx) {
    for (j=0, py=y; j<*n; j++, py+=*incy) {
      A[i + j*(*lda)] += *alpha * *px * *py;
    }
  }
}

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
		const int *ldc) {
  int i, j, l;
  double tmp;
  if ((*transa == 'N' || *transa == 'n') 
      && (*transb == 'N' || *transb == 'n')) {
    for (i=0; i<*m; i++) {
      for (j=0; j<*n; j++) {
	tmp = 0.0;
	for (l=0; l<*k; l++) {
	  tmp += A[l*(*lda)+i] * B[j*(*ldb)+l];
	}
	C[j*(*ldc)+i] = *alpha * tmp + *beta * C[j*(*ldc)+i];
      }
    }
  } else if ((*transa == 'T' || *transa == 't')
	     && (*transb == 'N' || *transb == 'n')) {
    for (i=0; i<*m; i++) {
      for (j=0; j<*n; j++) {
	tmp = 0.0;
	for (l=0; l<*k; l++) {
	  tmp += A[i*(*lda)+l] * B[j*(*ldb)+l];
	}
	C[j*(*ldc)+i] = *alpha * tmp + *beta * C[j*(*ldc)+i];
      }
    }
  } else if ((*transa == 'N' || *transa == 'n')
	     && (*transb == 'T' || *transb == 't')) {
    for (i=0; i<*m; i++) {
      for (j=0; j<*n; j++) {
	tmp = 0.0;
	for (l=0; l<*k; l++) {
	  tmp += A[l*(*lda)+i] * B[l*(*ldb)+j];
	}
	C[j*(*ldc)+i] = *alpha * tmp + *beta * C[j*(*ldc)+i];
      }
    }
  } else if ((*transa == 'T' || *transa == 't')
	     && (*transb == 'T' || *transb == 't')) {
    for (i=0; i<*m; i++) {
      for (j=0; j<*n; j++) {
	tmp = 0.0;
	for (l=0; l<*k; l++) {
	  tmp += A[i*(*lda)+l] * B[l*(*ldb)+j];
	}
	C[j*(*ldc)+i] = *alpha * tmp + *beta * C[j*(*ldc)+i];
      }
    }
  }
}

