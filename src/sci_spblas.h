#pragma once

#include <cfloat>

#include "sci_blas.h"
#include "sci_csrmatrix.h"
#include "sci_cscmatrix.h"
#include "sci_coomatrix.h"

namespace sci {

    // csrmatrix

    int darnoldi(mat::trans tr,
        int m, const spmatrix<double>& spA,
        const vector<double>& x, dmatrix<double>& H, dmatrix<double>& V,
        double& beta, int ite = 1, double tol = DBL_EPSILON);

    // level 1

    vector<double> diag(const csrmatrix<double>& m);
    double* dcopy(const csrmatrix<double>& m, double *p);
    csrmatrix<double>& dcopy(const double *p, csrmatrix<double>& m);
    csrmatrix<double>& dfill(double v, csrmatrix<double>& m);
    csrmatrix<double>& dscal(double alpha, csrmatrix<double>& m);
    csrmatrix<double>& daxpy(double alpha, const csrmatrix<double>& x, csrmatrix<double>& y);
    double dnrm2(const csrmatrix<double>& m);
    double damax(const csrmatrix<double>& m);
    
    // level 2
    
    vector<double>& dgemv(mat::trans tr, double alpha, const csrmatrix<double>& A,
     const vector<double>& x, double beta, vector<double>& y);

    csrmatrix<double>& dger(double alpha, const vector<double>& x, const vector<double>& y,
        csrmatrix<double>& A);
    csrmatrix<double>& dger(mat::trans tr, double alpha,
        const vector<double>& x, const vector<double>& y, csrmatrix<double>& A);

    // level 3

    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB,
        double alpha, const csrmatrix<double>& A, const dmatrix<double>& B,
        double beta, dmatrix<double>& C);

    // cscmatrix

    // level 1

    vector<double> diag(const cscmatrix<double>& m);
    double* dcopy(const cscmatrix<double>& m, double *p);
    cscmatrix<double>& dcopy(const double *p, cscmatrix<double>& m);
    cscmatrix<double>& dfill(double v, cscmatrix<double>& m);
    cscmatrix<double>& dscal(double alpha, cscmatrix<double>& m);
    cscmatrix<double>& daxpy(double alpha, const cscmatrix<double>& x, cscmatrix<double>& y);
    double dnrm2(const cscmatrix<double>& m);
    double damax(const cscmatrix<double>& m);
    
    // level 2
    
    vector<double>& dgemv(mat::trans tr, double alpha, const cscmatrix<double>& A,
     const vector<double>& x, double beta, vector<double>& y);

    cscmatrix<double>& dger(double alpha, const vector<double>& x, const vector<double>& y,
        cscmatrix<double>& A);
    cscmatrix<double>& dger(mat::trans tr, double alpha,
        const vector<double>& x, const vector<double>& y, cscmatrix<double>& A);

    // level 3

    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB,
        double alpha, const cscmatrix<double>& A, const dmatrix<double>& B,
        double beta, dmatrix<double>& C);

    // coomatrix

    // level 1

    vector<double> diag(const coomatrix<double>& m);
    double* dcopy(const coomatrix<double>& m, double *p);
    coomatrix<double>& dcopy(const double *p, coomatrix<double>& m);
    coomatrix<double>& dfill(double v, coomatrix<double>& m);
    coomatrix<double>& dscal(double alpha, coomatrix<double>& m);
    coomatrix<double>& daxpy(double alpha, const coomatrix<double>& x, coomatrix<double>& y);
    double dnrm2(const coomatrix<double>& m);
    double damax(const coomatrix<double>& m);
    
    // level 2
    
    vector<double>& dgemv(mat::trans tr, double alpha, const coomatrix<double>& A,
     const vector<double>& x, double beta, vector<double>& y);

    coomatrix<double>& dger(double alpha, const vector<double>& x, const vector<double>& y,
        coomatrix<double>& A);
    coomatrix<double>& dger(mat::trans tr, double alpha,
        const vector<double>& x, const vector<double>& y, coomatrix<double>& A);

    // level 3

    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB,
        double alpha, const coomatrix<double>& A, const dmatrix<double>& B,
        double beta, dmatrix<double>& C);

    //////////////////////////////////////////////////////
    // lapack for sparse
    //////////////////////////////////////////////////////

    vector<double>& dgsstep(mat::trans tr,
        double alpha, const csrmatrix<double>& A,
        const vector<double>& x, vector<double>& y);

    vector<double>& dgsstep(mat::trans tr,
        double alpha, const cscmatrix<double>& A,
        const vector<double>& x, vector<double>& y);

    vector<double>& dgsstep(mat::trans tr,
        double alpha, const coomatrix<double>& A,
        const vector<double>& x, vector<double>& y);

    int dgssolve(mat::trans tr,
        double alpha, const csrmatrix<double>& A,
        const vector<double>& x, vector<double>& y, int maxiter, double eps);

    int dgssolve(mat::trans tr,
        double alpha, const cscmatrix<double>& A,
        const vector<double>& x, vector<double>& y, int maxiter, double eps);

    int dgssolve(mat::trans tr,
        double alpha, const coomatrix<double>& A,
        const vector<double>& x, vector<double>& y, int maxiter, double eps);

    int darnoldi(mat::trans tr,
        int m, const csrmatrix<double>& spA,
        const vector<double>& x, dmatrix<double>& H, dmatrix<double>& V,
        double& beta, int ite = 1, double tol = DBL_EPSILON);

    int darnoldi(mat::trans tr,
        int m, const cscmatrix<double>& spA,
        const vector<double>& x, dmatrix<double>& H, dmatrix<double>& V,
        double& beta, int ite = 1, double tol = DBL_EPSILON);

    int darnoldi(mat::trans tr,
        int m, const coomatrix<double>& spA,
        const vector<double>& x, dmatrix<double>& H, dmatrix<double>& V,
        double& beta, int ite = 1, double tol = DBL_EPSILON);
}

