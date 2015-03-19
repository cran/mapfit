#pragma once

#include "sci_matrix.hh"

namespace sci {

    namespace mat {
        enum trans {
            T = 'T',
            N = 'N'
        };

        trans t(trans x);
    }

    double* dcopy(const vector<double>& y, double *p);
    vector<double>& dcopy(const vector<double>& x, vector<double>& y);
    vector<double>& dcopy(const double *p, vector<double>& y);
    vector<double>& dfill(double v, vector<double>& x);
    vector<double>& dscal(double alpha, vector<double>& x);
    vector<double>& daxpy(double alpha, const vector<double>& x, vector<double>& y);
    double ddot(const vector<double>& x, const vector<double>& y);
    double dsum(const vector<double>& x);
    double dasum(const vector<double>& x);
    double dnrm2(const vector<double>& x);
    int idmax(const vector<double>& x);
    int idamax(const vector<double>& x);
    double dmax(const vector<double>& x);
    double damax(const vector<double>& x);

    vector<double> diag(const matrix<double>& m);
    double* dcopy(const matrix<double>& y, double *p);
    matrix<double>* dnewcopy(const matrix<double>& x);
    matrix<double>* dnewcopy(const matrix<double>& x, double *v);
    matrix<double>& dcopy(const matrix<double>& x, matrix<double>& y);
    matrix<double>& dcopy(const double *p, matrix<double>& y);    
    matrix<double>& dfill(double v, matrix<double>& x);
    matrix<double>& dscal(double alpha, matrix<double>& x);
    matrix<double>& daxpy(double alpha, const matrix<double>& x, matrix<double>& y);
    double dnrm2(const matrix<double>& x);
    double damax(const matrix<double>& x);
    vector<double>& dgemv(mat::trans tr, double alpha, const matrix<double>& A, const vector<double>& x, double beta, vector<double>& y);
    matrix<double>& dger(double alpha, const vector<double>& x, const vector<double>& y, matrix<double>& A);
    matrix<double>& dger(mat::trans tr, double alpha, const vector<double>& x, const vector<double>& y, matrix<double>& A);
    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB,
                   double alpha, const matrix<double>& A,
                   const matrix<double>& B, double beta, dmatrix<double>& C);

    // dmatrix

    // level 1

    vector<double> diag(const dmatrix<double>& m);
    double* dcopy(const dmatrix<double>& m, double *p);
    dmatrix<double>& dcopy(const double *p, dmatrix<double>& m);
    dmatrix<double>& dcopy(mat::trans tr, const double *p, dmatrix<double>& m);
    dmatrix<double>& dcopy(const dmatrix<double>& x, dmatrix<double>& y);
    dmatrix<double>& dcopy(mat::trans tr, const dmatrix<double>& x, dmatrix<double>& y);
    dmatrix<double>& dfill(double v, dmatrix<double>& m);
    dmatrix<double>& dscal(double alpha, dmatrix<double>& m);
    dmatrix<double>& daxpy(double alpha, const dmatrix<double>& x, dmatrix<double>& y);
    double dnrm2(const dmatrix<double>& m);
    double damax(const dmatrix<double>& m);
    
    // level 2
    
    vector<double>& dgemv(mat::trans tr, double alpha, const dmatrix<double>& A,
     const vector<double>& x, double beta, vector<double>& y);

    dmatrix<double>& dger(double alpha, const vector<double>& x,
        const vector<double>& y, dmatrix<double>& A);

    dmatrix<double>& dger(mat::trans tr, double alpha,
        const vector<double>& x, const vector<double>& y, dmatrix<double>& A);

    // level 3
    
    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB,
        double alpha, const dmatrix<double>& A, const dmatrix<double>& B,
        double beta, dmatrix<double>& C);

    vector<double>& dkronmv(mat::trans tr, int n1, int n2,
        double alpha, const dmatrix<double>& A,
        const vector<double>& x, double beta, vector<double>& y);

    // lapack
    int dgesv(mat::trans tr, double alpha, const dmatrix<double>& A,
        const vector<double>& x, vector<double>& y);

    int dgesv(const dmatrix<double>& A, const dmatrix<double>& B, dmatrix<double>& C);

    dmatrix<double>& dpotrf(char uplo, dmatrix<double>& A);

    double dpodet(const dmatrix<double>& A);

    int dggglm(const dmatrix<double>& X, const dmatrix<double>& W,
        const vector<double>& z, vector<double>& beta, vector<double>& residual,
        int lwork, double *work);
}
