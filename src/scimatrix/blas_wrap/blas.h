/*
 blas.h
 wrapper for blas
 */

#ifndef _WRAP_BLAS_H
#define _WRAP_BLAS_H

#ifdef __cplusplus
 extern "C" {
#endif
    
    // level 1
    
    void blas_dcopy(int n, const double *x, int incx, double *y, int incy);
    void blas_dfill(int n, double v, double *x, int incx);
    void blas_dscal(int n, double alpha, double *x, int incx);
    void blas_daxpy(int n, double alpha, const double *x, int incx,
        double *y, int incy);
    
    double blas_dsum(int n, const double *x, int incx);
    double blas_dasum(int n, const double *x, int incx);
    int blas_idmax(int n, const double *x, int incx);
    int blas_idamax(int n, const double *x, int incx);
    double blas_dnrm2(int n, const double *x, int incx);
    double blas_ddot(int n, const double *x, int incx, const double *y, int incy);
    
    // level 2
    void blas_dgemv(char trans, int m, int n,
        double alpha,
        const double *A, int lda,
        const double *x, int incx,
        double beta, double *y, int incy);
    
    void blas_dger(int m, int n,
     double alpha,
     const double *x, int incx,
     const double *y, int incy,
     double *A, int lda);
    
    // level 3
    
    void blas_dgemm(char transa, char transb,
        int m, int n, int k, double alpha,
        const double *A, int lda, const double *B, int ldb,
        double beta, double *C, int ldc);

    // kron
    
    void blas_dkronmv(
        char trans,
        int m, int n, int n1, int n2,
        double alpha,
        const double *A, int lda,
        const double *x,
        double beta, double *y);

    void blas_dkronr(int m, int n, int n1, int n2,
            double alpha, double *x, double *y,
            double *A, int lda);

#ifdef __cplusplus
}
#endif

#endif
