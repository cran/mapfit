/*
 blas.h
 wrapper for blas
 */

#ifndef _WRAP_MYDGGGLM_H
#define _WRAP_MYDGGGLM_H

#ifdef __cplusplus
extern "C" {
#endif
    
  int mylapack_dggglm(int n, int m, int p, double *A, int lda,
                      double *B, int ldb, double *d, double *x,
                      double *y, double *work, int lwork);
  
#ifdef __cplusplus
}
#endif

#endif
