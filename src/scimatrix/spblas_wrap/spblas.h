/*
spblas.h
*/

#ifndef _WRAP_SPBLAS_H
#define _WRAP_SPBLAS_H

//#define SPBLAS_CBASE

#ifdef __cplusplus
extern "C" {
#endif

	void spblas_dcsrmv(char transa, int m, int n,
		double alpha, const double *val,
		const int *rowptr, const int *colind, int nnz,
		const double *x, int incx, double beta, double *y, int incy);

	void spblas_dcscmv(char transa, int m, int n,
		double alpha, const double *val,
		const int *colptr, const int *rowind, int nnz,
		const double *x, int incx, double beta, double *y, int incy);

	void spblas_dcoomv(char transa, int m, int n,
		double alpha, const double *val,
		const int *rowind, const int *colind, int nnz,
		const double *x, int incx, double beta, double *y, int incy);

	void spblas_dcsrr(char transa, int m, int n,
		double alpha, 
		const double *x, int incx,
		const double *y, int incy,
		double *val, const int *rowptr, const int *colind, int nnz);

	void spblas_dcscr(char transa, int m, int n,
		double alpha, 
		const double *x, int incx,
		const double *y, int incy,
		double *val, const int *colptr, const int *rowind, int nnz);

	void spblas_dcoor(char transa, int m, int n,
		double alpha, 
		const double *x, int incx,
		const double *y, int incy,
		double *val, const int *rowind, const int *colind, int nnz);

	void spblas_dcsrmm(char transa, char transb, int m, int n, int k,
		double alpha, const double *val,
		const int *rowptr, const int *colind, int nnz,
		const double *B, int ldb,
		double beta, double *C, int ldc);

	void spblas_dcscmm(char transa, char transb, int m, int n, int k,
		double alpha, const double *val,
		const int *colptr, const int *rowind, int nnz,
		const double *B, int ldb,
		double beta, double *C, int ldc);

	void spblas_dcoomm(char transa, char transb, int m, int n, int k,
		double alpha, const double *val,
		const int *rowind, const int *colind, int nnz,
		const double *B, int ldb,
		double beta, double *C, int ldc);

	void spblas_dkronmv(
		char trans,
		int m, int n, int n1, int n2,
		double alpha,
		const double *A, const int *rowptr, const int *colind,
		const double *x,
		double beta, double *y);

#ifdef __cplusplus
}
#endif

#endif
