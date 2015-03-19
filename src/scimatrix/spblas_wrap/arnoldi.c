/*
Arnoldi process
*/

#include <stdio.h>

#include "blas.h"
#include "spblas.h"

/*
  Description:

   Arnoldi process

     The matrix A is approximated by h on the subspace

         h = v^T A v on Krylov subspace  {x, A*x, A^2*x, ... A^(m-1)*x}

  Parameters:

   A (in): a n-by-n square sparse matrix
   x (in): a vector with n elements
   h (out): m-by-m square Hessenburg matrix
   v (out): n-by-m matrix of orthogonal vectors on Krylov subspace {x, A*x, A^2*x, ... A^(m-1)*x}

   rnorm (out): 2-norm of vector x (and is used for checking for happy breakdown)
   tol (in): tolerance error for happy breakdown
   ite (in): the number of times of iteration for Arnoldi process to reduce roundoff errors
             ite=1 is enough in most cases
   info (out): the size of effective dimensions of h when a happy breakdown occurs

   work: working area. The required size is n.

*/

void spblas_dcsrar(char trans, int n, int m,
  const double *A,
  const int *rowptr,
  const int *colind,
  int nnz,
  const double *x, int incx,
  double *h, int ldh,
  double *v, int ldv,
  int ite, double *rnorm,
  double tol,
  double *work, int *info) {

  int i, j, l;
  double beta, r;

  *rnorm = blas_dnrm2(nnz, A, 1);
  beta = blas_dnrm2(n, x, incx);

  for (i=0; i<m; i++) {
    blas_dfill(m, 0.0, &h[i*ldh], 1);
    blas_dfill(n, 0.0, &v[i*ldv], 1);
  }

  blas_daxpy(n, 1.0/beta, x, incx, &v[0], 1);
  for (j=0; j<m; j++) {
    spblas_dcsrmv(trans, n, n, 1.0, A, rowptr, colind, nnz,
      &v[j*ldv], 1, 0.0, work, 1);
    for (i=0; i<=j; i++) {
      for (l=0; l<ite; l++) {
        r = blas_ddot(n, work, 1, &v[i*ldv], 1);
        h[i+j*ldh] += r;
        blas_daxpy(n, -r, &v[i*ldv], 1, work, 1);
      }
    }
    if (j != m-1) {
      h[j+1+j*ldh] = blas_dnrm2(n, work, 1);
      if (h[j+1+j*ldh] < *rnorm * tol) {
        *rnorm = beta;
        *info = j+1;
        return;
      }
      blas_daxpy(n, 1.0/h[j+1+j*ldh], work, 1, &v[(j+1)*ldv], 1);
    }
  }
  *rnorm = beta;
  *info = j;
}

void spblas_dcscar(char trans, int n, int m,
  const double *A,
  const int *colptr,
  const int *rowind,
  int nnz,
  const double *x, int incx,
  double *h, int ldh,
  double *v, int ldv,
  int ite, double *rnorm,
  double tol,
  double *work, int *info) {

  int i, j, l;
  double beta, r;

  *rnorm = blas_dnrm2(nnz, A, 1);
  beta = blas_dnrm2(n, x, incx);

  for (i=0; i<m; i++) {
    blas_dfill(m, 0.0, &h[i*ldh], 1);
    blas_dfill(n, 0.0, &v[i*ldv], 1);
  }

  blas_daxpy(n, 1.0/beta, x, incx, &v[0], 1);
  for (j=0; j<m; j++) {
    spblas_dcscmv(trans, n, n, 1.0, A, colptr, rowind, nnz,
      &v[j*ldv], 1, 0.0, work, 1);
    for (i=0; i<=j; i++) {
      for (l=0; l<ite; l++) {
        r = blas_ddot(n, work, 1, &v[i*ldv], 1);
        h[i+j*ldh] += r;
        blas_daxpy(n, -r, &v[i*ldv], 1, work, 1);
      }
    }
    if (j != m-1) {
      h[j+1+j*ldh] = blas_dnrm2(n, work, 1);
      if (h[j+1+j*ldh] < *rnorm * tol) {
        *rnorm = beta;
        *info = j+1;
        return;
      }
      blas_daxpy(n, 1.0/h[j+1+j*ldh], work, 1, &v[(j+1)*ldv], 1);
    }
  }
  *rnorm = beta;
  *info = j;
}

void spblas_dcooar(char trans, int n, int m,
  const double *A,
  const int *rowind,
  const int *colind,
  int nnz,
  const double *x, int incx,
  double *h, int ldh,
  double *v, int ldv,
  int ite, double *rnorm,
  double tol,
  double *work, int *info) {

  int i, j, l;
  double beta, r;

  *rnorm = blas_dnrm2(nnz, A, 1);
  beta = blas_dnrm2(n, x, incx);

  for (i=0; i<m; i++) {
    blas_dfill(m, 0.0, &h[i*ldh], 1);
    blas_dfill(n, 0.0, &v[i*ldv], 1);
  }

  blas_daxpy(n, 1.0/beta, x, incx, &v[0], 1);
  for (j=0; j<m; j++) {
    spblas_dcoomv(trans, n, n, 1.0, A, rowind, colind, nnz,
      &v[j*ldv], 1, 0.0, work, 1);
    for (i=0; i<=j; i++) {
      for (l=0; l<ite; l++) {
        r = blas_ddot(n, work, 1, &v[i*ldv], 1);
        h[i+j*ldh] += r;
        blas_daxpy(n, -r, &v[i*ldv], 1, work, 1);
      }
    }
    if (j != m-1) {
      h[j+1+j*ldh] = blas_dnrm2(n, work, 1);
      if (h[j+1+j*ldh] < *rnorm * tol) {
        *rnorm = beta;
        *info = j+1;
        return;
      }
      blas_daxpy(n, 1.0/h[j+1+j*ldh], work, 1, &v[(j+1)*ldv], 1);
    }
  }
  *rnorm = beta;
  *info = j;
}

