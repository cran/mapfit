/*
  DPOTRF
  
  Description:
  Cholesky factorization of a real symmetric positive definete matrix A
  
  A = U^T * U if UPLO = 'U'
  A = L * L^T if UPLO = 'L'
  
  Arguments:
  uplo: 'U' or 'L' (see above)
  n: the number of columns (rows) of A
  A: a square matrix
  lda: LDA of A
  
  return value: interger of infomation
  info=0: success
  info<0: if INFO = -i, the i-th argument had an illegal value
  info>0: if INFO = i, the leading minor of order i is not
  positive definite, and the factorization could not be
  completed.
  
*/

#include <stdio.h>
#include <stdlib.h>

int mylapack_dpotrf(char uplo, int n, double *A, int lda) {
//  fprintf(stderr, "Please use f77lapack or MKL lapack\n");
//  exit(1);
  return -1;
}

