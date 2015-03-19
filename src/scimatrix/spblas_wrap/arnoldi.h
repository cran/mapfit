/*
  Arnoldi process
*/

#ifndef _WRAP_ARNOLDI_H
#define _WRAP_ARNOLDI_H

#ifdef __cplusplus
 extern "C" {
#endif
    
   void spblas_dcsrar(char trans, int n, int m,
		      const double *val,
		      const int *rowptr,
		      const int *colind,
		      int nnz,
		      const double *x, int incx,
		      double *h, int ldh,
		      double *v, int ldv,
		      int ite, double *rnorm,
		      double tol,
		      double *work, int *info);
   
   void spblas_dcscar(char trans, int n, int m,
		      const double *val,
		      const int *colptr,
		      const int *rowind,
		      int nnz,
		      const double *x, int incx,
		      double *h, int ldh,
		      double *v, int ldv,
		      int ite, double *rnorm,
		      double tol,
		      double *work, int *info);
   
   void spblas_dcooar(char trans, int n, int m,
		      const double *val,
		      const int *rowind,
		      const int *colind,
		      int nnz,
		      const double *x, int incx,
		      double *h, int ldh,
		      double *v, int ldv,
		      int ite, double *rnorm,
		      double tol,
		      double *work, int *info);

#ifdef __cplusplus
}
#endif

#endif
