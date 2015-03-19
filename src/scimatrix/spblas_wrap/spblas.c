/*
spblas
*/

#include "myspblas.h"

#ifdef MYSPBLAS
#define __DCSRMV__ myspblas_dcsrmv
#define __DCSCMV__ myspblas_dcscmv
#define __DCOOMV__ myspblas_dcoomv
#define __DCSRMM__ myspblas_dcsrmm
#define __DCSCMM__ myspblas_dcscmm
#define __DCOOMM__ myspblas_dcoomm
#endif

#ifdef MKLSPBLAS
/*----- prototype -------------*/
void mkl_dcsrmv_(const char *transa,
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
	double *y);

void mkl_dcscmv_(const char *transa,
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
	double *y);

void mkl_dcoomv_(const char *transa,
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
	double *y);

void mkl_dcsrmm_(const char *transa,
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
	const int *ldc);

void mkl_dcscmm_(const char *transa,
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
	const int *ldc);

void mkl_dcoomm_(const char *transa,
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
	const int *ldc);

/*------------------*/

#define __DCSRMV__ mkl_dcsrmv_
#define __DCSCMV__ mkl_dcscmv_
#define __DCOOMV__ mkl_dcoomv_
#define __DCSRMM__ mkl_dcsrmm_
#define __DCSCMM__ mkl_dcscmm_
#define __DCOOMM__ mkl_dcoomm_
#endif

// level 2
#ifdef SPBLAS_CBASE
const char matdescra[] = {'G', '-', '-', 'C', '-', '-'};
#else
const char matdescra[] = {'G', '-', '-', 'F', '-', '-'};
#endif

void spblas_dcsrmv(char transa, int m, int n,
	double alpha, const double *val,
	const int *rowptr, const int *colind, int nnz,
	const double *x, int incx, double beta, double *y, int incy) {
	if (incx == 1 && incy == 1) {
		__DCSRMV__(&transa, &m, &n, &alpha, matdescra,
			val, colind, &rowptr[0], &rowptr[1], x, &beta, y);
	} else {
		myspblas_dcsrmv_base(&transa, &m, &n, &alpha, matdescra,
			val, colind, &rowptr[0], &rowptr[1], x, &incx, &beta, y, &incy);
	}
}

void spblas_dcscmv(char transa, int m, int n,
	double alpha, const double *val,
	const int *rowptr, const int *rowind, int nnz,
	const double *x, int incx, double beta, double *y, int incy) {
	if (incx == 1 && incy == 1) {
		__DCSCMV__(&transa, &m, &n, &alpha, matdescra,
			val, rowind, &rowptr[0], &rowptr[1], x, &beta, y);
	} else {
		myspblas_dcscmv_base(&transa, &m, &n, &alpha, matdescra,
			val, rowind, &rowptr[0], &rowptr[1], x, &incx, &beta, y, &incy);
	}
}

void spblas_dcoomv(char transa, int m, int n,
	double alpha, const double *val,
	const int *rowind, const int *colind, int nnz,
	const double *x, int incx, double beta, double *y, int incy) {
	if (incx == 1 && incy == 1) {
		__DCOOMV__(&transa, &m, &n, &alpha, matdescra,
			val, rowind, colind, &nnz, x, &beta, y);
	} else {
		myspblas_dcoomv_base(&transa, &m, &n, &alpha, matdescra,
			val, rowind, colind, &nnz, x, &incx, &beta, y, &incy);
	}
}

void spblas_dcsrr(char transa, int m, int n,
	double alpha, 
	const double *x, int incx,
	const double *y, int incy,
	double *val, const int *rowptr, const int *colind, int nnz) {
	if (transa == 'N' || transa == 'n') {
		myspblas_dcsrr(&m, &n, 
			&alpha, x, &incx, y, &incy,
			val, rowptr, colind, &nnz);
	} else if (transa == 'T' || transa == 't') {
		myspblas_dcsrr(&m, &n, 
			&alpha, y, &incy, x, &incx,
			val, rowptr, colind, &nnz);
	}
}

void spblas_dcscr(char transa, int m, int n,
	double alpha, 
	const double *x, int incx,
	const double *y, int incy,
	double *val, const int *colptr, const int *rowind, int nnz) {
	if (transa == 'N' || transa == 'n') {
		myspblas_dcscr(&m, &n, 
			&alpha, x, &incx, y, &incy,
			val, colptr, rowind, &nnz);
	} else if (transa == 'T' || transa == 't') {
		myspblas_dcscr(&m, &n, 
			&alpha, y, &incy, x, &incx,
			val, colptr, rowind, &nnz);
	}
}

void spblas_dcoor(char transa, int m, int n,
	double alpha, 
	const double *x, int incx,
	const double *y, int incy,
	double *val, const int *rowind, const int *colind, int nnz) {
	if (transa == 'N' || transa == 'n') {
		myspblas_dcoor(&m, &n, 
			&alpha, x, &incx, y, &incy,
			val, rowind, colind, &nnz);
	} else if (transa == 'T' || transa == 't') {
		myspblas_dcoor(&m, &n, 
			&alpha, y, &incy, x, &incx,
			val, rowind, colind, &nnz);
	}
}

// level 3

void spblas_dcsrmm(char transa, char transb, int m, int n, int k,
	double alpha, const double *val,
	const int *rowptr, const int *colind, int nnz,
	const double *B, int ldb,
	double beta, double *C, int ldc) {
	if (transb == 'N' || transb == 'n') {
		__DCSRMM__(&transa, &m, &n, &k, &alpha, matdescra,
			val, colind, &rowptr[0], &rowptr[1],
			B, &ldb, &beta, C, &ldc);
	} else if (transb == 'T' || transb == 't') {
		myspblas_dcsrmm_base(&transa, &transb, &m, &n, &k, &alpha, matdescra,
			val, colind, &rowptr[0], &rowptr[1],
			B, &ldb, &beta, C, &ldc);
	}
}

void spblas_dcscmm(char transa, char transb, int m, int n, int k,
	double alpha, const double *val,
	const int *colptr, const int *rowind, int nnz,
	const double *B, int ldb,
	double beta, double *C, int ldc) {
	if (transb == 'N' || transb == 'n') {
		__DCSCMM__(&transa, &m, &n, &k, &alpha, matdescra,
			val, rowind, &colptr[0], &colptr[1],
			B, &ldb, &beta, C, &ldc);
	} else {
		myspblas_dcscmm_base(&transa, &transb, &m, &n, &k, &alpha, matdescra,
			val, rowind, &colptr[0], &colptr[1],
			B, &ldb, &beta, C, &ldc);
	}
}

void spblas_dcoomm(char transa, char transb, int m, int n, int k,
	double alpha, const double *val,
	const int *rowind, const int *colind, int nnz,
	const double *B, int ldb,
	double beta, double *C, int ldc) {
	if (transb == 'N' || transb == 'n') {
		__DCOOMM__(&transa, &m, &n, &k, &alpha, matdescra,
			val, rowind, colind, &nnz,
			B, &ldb, &beta, C, &ldc);
	} else {
		myspblas_dcoomm_base(&transa, &transb, &m, &n, &k, &alpha, matdescra,
			val, rowind, colind, &nnz,
			B, &ldb, &beta, C, &ldc);
	}
}
