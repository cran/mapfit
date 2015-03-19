/*
  blas.h
 wrapper for blas
 */

#ifndef _WRAP_MYDPOTRF_H
#define _WRAP_MYDPOTRF_H

#ifdef __cplusplus
extern "C" {
#endif
    
  int mylapack_dpotrf(char uplo, int n, double *A, int lda);
    
#ifdef __cplusplus
}
#endif

#endif
