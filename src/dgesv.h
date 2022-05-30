
#ifndef _WRAP_MYDGESV_H
#define _WRAP_MYDGESV_H

#ifdef __cplusplus
extern "C" {
#endif
  
  int mylapack_dgesv(int n, int nrhs, double *A, int lda,
		     int *ipiv, double *B, int ldb);
    
#ifdef __cplusplus
}
#endif

#endif
