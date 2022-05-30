/*
 sci blas
 */

#include <cmath>

#include "sci_blas.h"
#include "sci_spblas.h"

namespace sci {

    namespace mat {
        trans t(trans x) {
            if (x == T) {
                return N;
            } else {
                return T;
            }
        }        
    }    

    // level 1

    vector<double> diag(const matrix<double>& m) {
        switch (m.type()) {
        case (DENSE):
            return diag(dynamic_cast<const dmatrix<double>&>(m));
        case (CSR):
            return diag(dynamic_cast<const csrmatrix<double>&>(m));
        case (CSC):
            return diag(dynamic_cast<const cscmatrix<double>&>(m));
        case (COO):
            return diag(dynamic_cast<const coomatrix<double>&>(m));
        default:
            throw;
        }
    }

    double* dcopy(const matrix<double>& x, double *p) {
        switch (x.type()) {
        case (DENSE):
            return dcopy(dynamic_cast<const dmatrix<double>&>(x), p);
        case (CSR):
            return dcopy(dynamic_cast<const csrmatrix<double>&>(x), p);
        case (CSC):
            return dcopy(dynamic_cast<const cscmatrix<double>&>(x), p);
        case (COO):
            return dcopy(dynamic_cast<const coomatrix<double>&>(x), p);
        default:
            throw;
        }
    }

    matrix<double>* dnewcopy(const matrix<double>& x) {
        switch (x.type()) {
        case (DENSE):
            return new dmatrix<double>(dynamic_cast<const dmatrix<double>&>(x));
        case (CSR):
            return new csrmatrix<double>(dynamic_cast<const csrmatrix<double>&>(x));
        case (CSC):
            return new cscmatrix<double>(dynamic_cast<const cscmatrix<double>&>(x));
        case (COO):
            return new coomatrix<double>(dynamic_cast<const coomatrix<double>&>(x));
        default:
            throw;
        }
    }

    matrix<double>* dnewcopy(const matrix<double>& m, double *v) {
        switch (m.type()) {
        case (DENSE): {
            const sci::dmatrix<double>& mm(dynamic_cast<const sci::dmatrix<double>&>(m));
            sci::dmatrix<double>* retm = new sci::dmatrix<double>(mm.nrow, mm.ncol, v);
            return retm;
        }
        case (CSR): {
            const sci::csrmatrix<double>& mm(dynamic_cast<const sci::csrmatrix<double>&>(m));
            sci::csrmatrix<double>* retm = new sci::csrmatrix<double>(mm.nrow, mm.ncol, mm.nnz, v);
            retm->rowptr = mm.rowptr;
            retm->colind = mm.colind;
            return retm;
        }
        case (CSC): {
            const sci::cscmatrix<double>& mm(dynamic_cast<const sci::cscmatrix<double>&>(m));
            sci::cscmatrix<double>* retm = new sci::cscmatrix<double>(mm.nrow, mm.ncol, mm.nnz, v);
            retm->colptr = mm.colptr;
            retm->rowind = mm.rowind;
            return retm;
        }
        case (COO): {
            const sci::coomatrix<double>& mm(dynamic_cast<const sci::coomatrix<double>&>(m));
            sci::coomatrix<double>* retm = new sci::coomatrix<double>(mm.nrow, mm.ncol, mm.nnz, v);
            retm->rowind = mm.rowind;
            retm->colind = mm.colind;
            return retm;
        }
        default:
            throw;
        }
    }

    matrix<double>& dcopy(const double *p, matrix<double>& y) {
        switch (y.type()) {
        case (DENSE):
            return dcopy(p, dynamic_cast<dmatrix<double>&>(y));
        case (CSR):
            return dcopy(p, dynamic_cast<csrmatrix<double>&>(y));
        case (CSC):
            return dcopy(p, dynamic_cast<cscmatrix<double>&>(y));
        case (COO):
            return dcopy(p, dynamic_cast<coomatrix<double>&>(y));
        default:
            throw;
        }
    }

    matrix<double>& dcopy(const matrix<double>& x, matrix<double>& y) {
        switch (PAIR(x.type(),y.type())) {
        case (PAIR(DENSE,DENSE)):
            return dcopy(dynamic_cast<const dmatrix<double>&>(x), dynamic_cast<dmatrix<double>&>(y));
        case (PAIR(CSR,CSR)):
            return dcopy(dynamic_cast<const csrmatrix<double>&>(x), dynamic_cast<csrmatrix<double>&>(y));
        case (PAIR(CSC,CSC)):
            return dcopy(dynamic_cast<const cscmatrix<double>&>(x), dynamic_cast<cscmatrix<double>&>(y));
        case (PAIR(COO,COO)):
            return dcopy(dynamic_cast<const coomatrix<double>&>(x), dynamic_cast<coomatrix<double>&>(y));
        default:
            throw;
        }
    }

    matrix<double>& dfill(double v, matrix<double>& x) {
        switch (x.type()) {
        case (DENSE):
            return dfill(v, dynamic_cast<dmatrix<double>&>(x));
        case (CSR):
            return dfill(v, dynamic_cast<csrmatrix<double>&>(x));
        case (CSC):
            return dfill(v, dynamic_cast<cscmatrix<double>&>(x));
        case (COO):
            return dfill(v, dynamic_cast<coomatrix<double>&>(x));
        default:
            throw;
        }
    }

    matrix<double>& dscal(double alpha, matrix<double>& x) {
        switch (x.type()) {
        case (DENSE):
            return dscal(alpha, dynamic_cast<dmatrix<double>&>(x));
        case (CSR):
            return dscal(alpha, dynamic_cast<csrmatrix<double>&>(x));
        case (CSC):
            return dscal(alpha, dynamic_cast<cscmatrix<double>&>(x));
        case (COO):
            return dscal(alpha, dynamic_cast<coomatrix<double>&>(x));
        default:
            throw;
        }
    }

    matrix<double>& daxpy(double alpha, const matrix<double>& x, matrix<double>& y) {
        switch (PAIR(x.type(),y.type())) {
        case (PAIR(DENSE,DENSE)):
            return daxpy(alpha, dynamic_cast<const dmatrix<double>&>(x), dynamic_cast<dmatrix<double>&>(y));
        case (PAIR(CSR,CSR)):
            return daxpy(alpha, dynamic_cast<const csrmatrix<double>&>(x), dynamic_cast<csrmatrix<double>&>(y));
        case (PAIR(CSC,CSC)):
            return daxpy(alpha, dynamic_cast<const cscmatrix<double>&>(x), dynamic_cast<cscmatrix<double>&>(y));
        case (PAIR(COO,COO)):
            return daxpy(alpha, dynamic_cast<const coomatrix<double>&>(x), dynamic_cast<coomatrix<double>&>(y));
        default:
            throw;
        }
    }

    double dnrm2(const matrix<double>& m) {
        switch (m.type()) {
        case (DENSE):
            return dnrm2(dynamic_cast<const dmatrix<double>&>(m));
        case (CSR):
            return dnrm2(dynamic_cast<const csrmatrix<double>&>(m));
        case (CSC):
            return dnrm2(dynamic_cast<const cscmatrix<double>&>(m));
        case (COO):
            return dnrm2(dynamic_cast<const coomatrix<double>&>(m));
        default:
            throw;
        }
    }

    double damax(const matrix<double>& m) {
        switch (m.type()) {
        case (DENSE):
            return damax(dynamic_cast<const dmatrix<double>&>(m));
        case (CSR):
            return damax(dynamic_cast<const csrmatrix<double>&>(m));
        case (CSC):
            return damax(dynamic_cast<const cscmatrix<double>&>(m));
        case (COO):
            return damax(dynamic_cast<const coomatrix<double>&>(m));
        default:
            throw;
        }
    }

    vector<double>& dgemv(mat::trans tr, double alpha, const matrix<double>& A,
     const vector<double>& x, double beta, vector<double>& y) {
        switch (A.type()) {
        case (DENSE):
            return dgemv(tr, alpha, dynamic_cast<const dmatrix<double>&>(A), x, beta, y);
        case (CSR):
            return dgemv(tr, alpha, dynamic_cast<const csrmatrix<double>&>(A), x, beta, y);
        case (CSC):
            return dgemv(tr, alpha, dynamic_cast<const cscmatrix<double>&>(A), x, beta, y);
        case (COO):
            return dgemv(tr, alpha, dynamic_cast<const coomatrix<double>&>(A), x, beta, y);
        default:
            throw;            
        }
    }

    matrix<double>& dger(double alpha, const vector<double>& x,
        const vector<double>& y, matrix<double>& A) {
        switch (A.type()) {
        case (DENSE):
            return dger(alpha, x, y, dynamic_cast<dmatrix<double>&>(A));
        case (CSR):
            return dger(alpha, x, y, dynamic_cast<csrmatrix<double>&>(A));
        case (CSC):
            return dger(alpha, x, y, dynamic_cast<cscmatrix<double>&>(A));
        case (COO):
            return dger(alpha, x, y, dynamic_cast<coomatrix<double>&>(A));
        default:
            throw;
        }
    }

    matrix<double>& dger(mat::trans tr, double alpha, const vector<double>& x,
        const vector<double>& y, matrix<double>& A) {
        switch (A.type()) {
        case (DENSE):
            return dger(tr, alpha, x, y, dynamic_cast<dmatrix<double>&>(A));
        case (CSR):
            return dger(tr, alpha, x, y, dynamic_cast<csrmatrix<double>&>(A));
        case (CSC):
            return dger(tr, alpha, x, y, dynamic_cast<cscmatrix<double>&>(A));
        case (COO):
            return dger(tr, alpha, x, y, dynamic_cast<coomatrix<double>&>(A));
        default:
            throw;
        }
    }

    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB, double alpha,
        const matrix<double>& A, const matrix<double>& B, double beta, dmatrix<double>& C) {
        switch (PAIR(A.type(),B.type())) {
        case (PAIR(DENSE,DENSE)):
            return dgemm(trA, trB, alpha, dynamic_cast<const dmatrix<double>&>(A),
                dynamic_cast<const dmatrix<double>&>(B), beta, C);
        case (PAIR(CSR,DENSE)):
            return dgemm(trA, trB, alpha, dynamic_cast<const csrmatrix<double>&>(A),
                dynamic_cast<const dmatrix<double>&>(B), beta, C);
        case (PAIR(CSC,DENSE)):
            return dgemm(trA, trB, alpha, dynamic_cast<const cscmatrix<double>&>(A),
                dynamic_cast<const dmatrix<double>&>(B), beta, C);
        case (PAIR(COO,DENSE)):
            return dgemm(trA, trB, alpha, dynamic_cast<const coomatrix<double>&>(A),
                dynamic_cast<const dmatrix<double>&>(B), beta, C);
        default:
            throw;
        }
    }
   
    int darnoldi(mat::trans tr,
        int m, const spmatrix<double>& A,
        const vector<double>& x, dmatrix<double>& H, dmatrix<double>& V,
        double& beta, int ite, double tol) {
        switch (A.type()) {
        case (CSR):
            return darnoldi(tr, m, dynamic_cast<const csrmatrix<double>&>(A), x, H, V, beta, ite, tol);
        case (CSC):
            return darnoldi(tr, m, dynamic_cast<const cscmatrix<double>&>(A), x, H, V, beta, ite, tol);
        case (COO):
            return darnoldi(tr, m, dynamic_cast<const coomatrix<double>&>(A), x, H, V, beta, ite, tol);
        default:
            throw;
        }
    }

    // vector<double>& dcgnr(mat::trans tr, const matrix<double>& A,
    //     const vector<double>& b, const vector<double>& x, vector<double>& y, double tol) {

    //     double alpha0, alpha1, alpha;
    //     mat::trans ttr = mat::t(tr);
    //     vector<double> r(x.size);
    //     vector<double> p(x.size);
    //     vector<double> z(x.size);
    //     vector<double> tmpv(b.size);

    //     // r = A^T * b
    //     dgemv(ttr, 1.0, A, b, 0.0, r);

    //     // r = A^T * b - A^T * A * x
    //     dgemv(tr, 1.0, A, x, 0.0, tmpv); 
    //     dgemv(ttr, -1.0, A, tmpv, 1.0, r);

    //     p = r;
    //     y = x;

    //     for (int k=1; k<=x.size; k++) {

    //         // z = A^T * A * p
    //         dgemv(tr, 1.0, A, p, 0.0, tmpv);
    //         dgemv(ttr, 1.0, A, tmpv, 0.0, z);

    //         // alpha = (r,r)/(r,A*r)
    //         alpha0 = ddot(r, r);
    //         alpha =  alpha0 / ddot(z, p);

    //         // y = y + alpha * p
    //         daxpy(alpha, p, y);

    //         // r = r - alpha * z
    //         daxpy(-alpha, z, r);

    //         if (dnrm2(r) < tol) {
    //             return y;
    //         }

    //         alpha1 = ddot(r, r);
    //         dscal(alpha1 / alpha0, p);
    //         p += r;
    //         alpha0 = alpha1;
    //     }

    //     return y;
    // }


}

