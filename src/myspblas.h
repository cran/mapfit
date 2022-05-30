/* 
myspblas
*/

#ifndef _MYSPBLAS_H_
#define _MYSPBLAS_H_

#ifdef __cplusplus
extern "C" {
#endif

	void myspblas_dcsrmv_base(const char *transa,
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
		const int *incy);

	void myspblas_dcscmv_base(const char *transa,
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
		const int *incy);

	void myspblas_dcoomv_base(const char *transa,
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
		const int *incy);

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
		const int *ldc);

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
		const int *ldc);

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
		const int *ldc);

	void myspblas_dcsrmv(const char *transa,
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

	void myspblas_dcscmv(const char *transa,
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

	void myspblas_dcoomv(const char *transa,
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
		const int *nnz);

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
		const int *nnz);

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
		const int *nnz);

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
		const int *ldc);

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
		const int *ldc);

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
		const int *ldc);

#ifdef __cplusplus
}
#endif

#endif
