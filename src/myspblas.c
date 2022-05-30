/*
myspblas
*/

#include <stdlib.h>
#include <math.h>

#include "blas.h"
#include "spblas.h"

// level 2

void myspblas_dcsrmv_base(
	const char *transa,
	const int *m,
	const int *n,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *indx,
	const int *pntrb,
	const int *pntre,
	const double *x,
	const int *incx,
	const double *beta,
	double *y,
	const int *incy) {
	int i, j;
#ifdef SPBLAS_CBASE
	const int BASE = 0;
#else
	const int BASE = 1;
#endif
	if (*transa == 'N' || *transa == 'n') {
		blas_dscal(*m, *beta, y, *incy);
		for (i=0; i<*m; i++) {
			for (j=pntrb[i]-BASE; j<pntre[i]-BASE; j++) {
				y[i*(*incy)] += *alpha * val[j] * x[(indx[j]-BASE)*(*incx)];
			}
		}
	} else if(*transa == 'T' || *transa == 't') {
		blas_dscal(*n, *beta, y, *incy);
		for (i=0; i<*m; i++) {
			for (j=pntrb[i]-BASE; j<pntre[i]-BASE; j++) {
				y[(indx[j]-BASE)*(*incy)] += *alpha * val[j] * x[i*(*incx)];
			}
		}
	}
}

void myspblas_dcscmv_base(
	const char *transa,
	const int *m,
	const int *n,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *indx,
	const int *pntrb,
	const int *pntre,
	const double *x,
	const int *incx,
	const double *beta,
	double *y,
	const int *incy) {
	int i, j;
#ifdef SPBLAS_CBASE
	const int BASE = 0;
#else
	const int BASE = 1;
#endif
	if (*transa == 'N' || *transa == 'n') {
		blas_dscal(*m, *beta, y, *incy);
		for (i=0; i<*n; i++) {
			for (j=pntrb[i]-BASE; j<pntre[i]-BASE; j++) {
				y[(indx[j]-BASE)*(*incy)] += *alpha * val[j] * x[i*(*incx)];
			}
		}
	} else if(*transa == 'T' || *transa == 't') {
		blas_dscal(*n, *beta, y, *incy);
		for (i=0; i<*n; i++) {
			for (j=pntrb[i]-BASE; j<pntre[i]-BASE; j++) {
				y[i*(*incy)] += *alpha * val[j] * x[(indx[j]-BASE)*(*incx)];
			}
		}
	}
}

void myspblas_dcoomv_base(
	const char *transa,
	const int *m,
	const int *n,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *rowind,
	const int *colind,
	const int *nnz,
	const double *x,
	const int *incx,
	const double *beta,
	double *y,
	const int *incy) {
	int i;
#ifdef SPBLAS_CBASE
	const int BASE = 0;
#else
	const int BASE = 1;
#endif
	if (*transa == 'N' || *transa == 'n') {
		blas_dscal(*m, *beta, y, *incy);
		for (i=0; i<*nnz; i++) {
			y[(rowind[i]-BASE)*(*incy)] += *alpha * val[i] * x[(colind[i]-BASE)*(*incx)];
		}
	} else if(*transa == 'T' || *transa == 't') {
		blas_dscal(*n, *beta, y, *incy);
		for (i=0; i<*nnz; i++) {
			y[(colind[i]-BASE)*(*incy)] += *alpha * val[i] * x[(rowind[i]-BASE)*(*incx)];
		}
	}
}

/////

void myspblas_dcsrmv(
	const char *transa,
	const int *m,
	const int *n,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *indx,
	const int *pntrb,
	const int *pntre,
	const double *x,
	const double *beta,
	double *y) {
	int incx = 1;
	int incy = 1;
	myspblas_dcsrmv_base(transa, m, n,
		alpha, matdescra, val, indx, pntrb, pntre, x, &incx, beta, y, &incy);
}

void myspblas_dcscmv(
	const char *transa,
	const int *m,
	const int *n,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *indx,
	const int *pntrb,
	const int *pntre,
	const double *x,
	const double *beta,
	double *y) {
	int incx = 1;
	int incy = 1;
	myspblas_dcscmv_base(transa, m, n,
		alpha, matdescra, val, indx, pntrb, pntre, x, &incx, beta, y, &incy);
}

void myspblas_dcoomv(
	const char *transa,
	const int *m,
	const int *n,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *rowind,
	const int *colind,
	const int *nnz,
	const double *x,
	const double *beta,
	double *y) {
	int incx = 1;
	int incy = 1;
	myspblas_dcoomv_base(transa, m, n,
		alpha, matdescra, val, rowind, colind, nnz, x, &incx, beta, y, &incy);
}

void myspblas_dcsrr(
	const int *m,
	const int *n,
	const double *alpha, 
	const double *x,
	const int *incx,
	const double *y,
	const int *incy,
	double *val,
	const int *rowptr,
	const int *colind,
	const int *nnz) {
	int i, j, ip, jp;
#ifdef SPBLAS_CBASE
	const int BASE = 0;
#else
	const int BASE = 1;
#endif
	for (i=0, ip=0; i<*m; i++, ip+=*incx) {
		for (j=rowptr[i]-BASE; j<rowptr[i+1]-BASE; j++) {
			jp = (colind[j]-BASE) * (*incy);
			val[j] += *alpha * x[ip] * y[jp];
		}
	}
}

void myspblas_dcscr(
	const int *m,
	const int *n,
	const double *alpha, 
	const double *x,
	const int *incx,
	const double *y,
	const int *incy,
	double *val,
	const int *colptr,
	const int *rowind,
	const int *nnz) {
	int i, j, ip, jp;
#ifdef SPBLAS_CBASE
	const int BASE = 0;
#else
	const int BASE = 1;
#endif
	for (j=0, jp=0; j<*n; j++, jp+=*incy) {
		for (i=colptr[j]-BASE; i<colptr[j+1]-BASE; i++) {
			ip = (rowind[i]-BASE) * (*incx);
			val[i] += *alpha * x[ip] * y[jp];
		}
	}
}

void myspblas_dcoor(
	const int *m,
	const int *n,
	const double *alpha, 
	const double *x,
	const int *incx,
	const double *y,
	const int *incy,
	double *val,
	const int *rowind,
	const int *colind,
	const int *nnz) {
	int i, j, ip, jp;
#ifdef SPBLAS_CBASE
	const int BASE = 0;
#else
	const int BASE = 1;
#endif
	for (i=0; i<*nnz; i++) {
		val[i] += *alpha * x[(rowind[i]-BASE) * (*incx)] * y[(colind[i]-BASE) * (*incy)];
	}
}

// level 3

void myspblas_dcsrmm_base(const char *transa,
	const char *transb,
	const int *m,
	const int *n,
	const int *k,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *indx,
	const int *pntrb,
	const int *pntre,
	const double *B,
	const int *ldb,
	const double *beta,
	double *C,
	const int *ldc) {
	int i, incx, incy;
	if (*transb == 'N' || *transb == 'n') {
		incx = 1;
		incy = 1;
		if (*transa == 'N' || *transa == 'n') {
			for (i=0; i<*n; i++) {
				myspblas_dcsrmv_base(transa, m, k, alpha, matdescra, val, indx, pntrb, pntre,
					&B[i*(*ldb)], &incx, beta, &C[i*(*ldc)], &incy);
			}
		} else if (*transa == 'T' || *transa == 't') {
			for (i=0; i<*n; i++) {
				myspblas_dcsrmv_base(transa, k, m, alpha, matdescra, val, indx, pntrb, pntre,
					&B[i*(*ldb)], &incx, beta, &C[i*(*ldc)], &incy);
			}
		}
	} else if (*transb == 'T' || *transb == 't') {
		incx = *ldb;
		incy = 1;
		if (*transa == 'N' || *transa == 'n') {		
			for (i=0; i<*n; i++) {
				myspblas_dcsrmv_base(transa, m, k, alpha, matdescra, val, indx, pntrb, pntre,
					&B[i], &incx, beta, &C[i*(*ldc)], &incy);
			}
		} else if (*transa == 'T' || *transa == 't') {
			for (i=0; i<*n; i++) {
				myspblas_dcsrmv_base(transa, k, m, alpha, matdescra, val, indx, pntrb, pntre,
					&B[i], &incx, beta, &C[i*(*ldc)], &incy);
			}
		}
	}
}

void myspblas_dcscmm_base(const char *transa,
	const char *transb,
	const int *m,
	const int *n,
	const int *k,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *indx,
	const int *pntrb,
	const int *pntre,
	const double *B,
	const int *ldb,
	const double *beta,
	double *C,
	const int *ldc) {
	int i, incx, incy;
	if (*transb == 'N' || *transb == 'n') {
		incx = 1;
		incy = 1;
		if (*transa == 'N' || *transa == 'n') {
			for (i=0; i<*n; i++) {
				myspblas_dcscmv_base(transa, m, k, alpha, matdescra, val, indx, pntrb, pntre,
					&B[i*(*ldb)], &incx, beta, &C[i*(*ldc)], &incy);
			}
		} else if (*transa == 'T' || *transa == 't') {
			for (i=0; i<*n; i++) {
				myspblas_dcscmv_base(transa, k, m, alpha, matdescra, val, indx, pntrb, pntre,
					&B[i*(*ldb)], &incx, beta, &C[i*(*ldc)], &incy);
			}
		}
	} else if (*transb == 'T' || *transb == 't') {
		incx = *ldb;
		incy = 1;
		if (*transa == 'N' || *transa == 'n') {		
			for (i=0; i<*n; i++) {
				myspblas_dcscmv_base(transa, m, k, alpha, matdescra, val, indx, pntrb, pntre,
					&B[i], &incx, beta, &C[i*(*ldc)], &incy);
			}
		} else if (*transa == 'T' || *transa == 't') {
			for (i=0; i<*n; i++) {
				myspblas_dcscmv_base(transa, k, m, alpha, matdescra, val, indx, pntrb, pntre,
					&B[i], &incx, beta, &C[i*(*ldc)], &incy);
			}
		}
	}
}

void myspblas_dcoomm_base(const char *transa,
	const char *transb,
	const int *m,
	const int *n,
	const int *k,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *rowind,
	const int *colind,
	const int *nnz,
	const double *B,
	const int *ldb,
	const double *beta,
	double *C,
	const int *ldc) {
	int i, incx, incy;
	if (*transb == 'N' || *transb == 'n') {
		incx = 1;
		incy = 1;
		if (*transa == 'N' || *transa == 'n') {
			for (i=0; i<*n; i++) {
				myspblas_dcoomv_base(transa, m, k, alpha, matdescra, val, rowind, colind, nnz,
					&B[i*(*ldb)], &incx, beta, &C[i*(*ldc)], &incy);
			}
		} else if (*transa == 'T' || *transa == 't') {
			for (i=0; i<*n; i++) {
				myspblas_dcoomv_base(transa, k, m, alpha, matdescra, val, rowind, colind, nnz,
					&B[i*(*ldb)], &incx, beta, &C[i*(*ldc)], &incy);
			}
		}
	} else if (*transb == 'T' || *transb == 't') {
		incx = *ldb;
		incy = 1;
		if (*transa == 'N' || *transa == 'n') {		
			for (i=0; i<*n; i++) {
				myspblas_dcoomv_base(transa, m, k, alpha, matdescra, val, rowind, colind, nnz,
					&B[i], &incx, beta, &C[i*(*ldc)], &incy);
			}
		} else if (*transa == 'T' || *transa == 't') {
			for (i=0; i<*n; i++) {
				myspblas_dcoomv_base(transa, k, m, alpha, matdescra, val, rowind, colind, nnz,
					&B[i], &incx, beta, &C[i*(*ldc)], &incy);
			}
		}
	}
}

void myspblas_dcsrmm(const char *transa,
	const int *m,
	const int *n,
	const int *k,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *indx,
	const int *pntrb,
	const int *pntre,
	const double *B,
	const int *ldb,
	const double *beta,
	double *C,
	const int *ldc) {
	char transb = 'N';
	myspblas_dcsrmm_base(transa, &transb, m, n, k, alpha, matdescra, val,
		indx, pntrb, pntre, B, ldb, beta, C, ldc);
}

void myspblas_dcscmm(const char *transa,
	const int *m,
	const int *n,
	const int *k,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *indx,
	const int *pntrb,
	const int *pntre,
	const double *B,
	const int *ldb,
	const double *beta,
	double *C,
	const int *ldc) {
	char transb = 'N';
	myspblas_dcscmm_base(transa, &transb, m, n, k, alpha, matdescra, val,
		indx, pntrb, pntre, B, ldb, beta, C, ldc);
}

void myspblas_dcoomm(const char *transa,
	const int *m,
	const int *n,
	const int *k,
	const double *alpha,
	const char *matdescra,
	const double *val,
	const int *rowind,
	const int *colind,
	const int *nnz,
	const double *B,
	const int *ldb,
	const double *beta,
	double *C,
	const int *ldc) {
	char transb = 'N';
	myspblas_dcoomm_base(transa, &transb, m, n, k, alpha, matdescra, val,
		rowind, colind, nnz, B, ldb, beta, C, ldc);
}
