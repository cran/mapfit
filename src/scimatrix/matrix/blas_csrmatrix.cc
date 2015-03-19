#include <cmath>

#include "blas.h"
#include "spblas.h"
#include "arnoldi.h"

#include "sci_spblas.hh"
#include "sci_csrmatrix.hh"

namespace sci {

    vector<double> diag(const csrmatrix<double>& m) {
        if (m.nrow != m.ncol) {
            throw;
        }
        vector<double> result(m.nrow);
        array<int>::const_range rowptr = m.get_rowptr();
        array<int>::const_range colind = m.get_colind();
        array<double>::const_range value = m.get_value();
        
        for (size_t i=1; i<=m.nrow; i++) {
            for (int j=rowptr[i]; j<rowptr[i+1]; j++) {
                if (i == colind[j]) {
                    result(i) = value[j];
                    break;
                }
                if (i < colind[j]) {
                    break;
                }
            }
        }
        return result;
    }
    
    double* dcopy(const csrmatrix<double>& m, double *p) {
        blas_dcopy(m.nnz, m.ptr, 1, p, 1);
        return p;
    }

    csrmatrix<double>& dcopy(const double *p, csrmatrix<double>& m) {
        blas_dcopy(m.nnz, p, 1, m.ptr, 1);
        return m;
    }

    csrmatrix<double>& dfill(double v, csrmatrix<double>& m) {
        blas_dfill(m.nnz, v, m.ptr, 1);
        return m;
    }

    csrmatrix<double>& dscal(double alpha, csrmatrix<double>& m) {
        blas_dscal(m.nnz, alpha, m.ptr, 1);
        return m;
    }

    csrmatrix<double>& daxpy(double alpha, const csrmatrix<double>& x, csrmatrix<double>& y) {
        if (x.nrow != y.nrow || x.ncol != y.ncol || x.nnz != y.nnz) {
            throw;
        }
        blas_daxpy(x.nnz, alpha, x.ptr, 1, y.ptr, 1);
        return y;
    }

    double dnrm2(const csrmatrix<double>& m) {
        double sum = 0.0;
        for (size_t i=0; i<m.nnz; i++) {
            sum += m.ptr[i] * m.ptr[i];
        }
        return sqrt(sum);
    }
    
    double damax(const csrmatrix<double>& m) {
        double max = 0.0;
        for (size_t i=0; i<m.nnz; i++) {
            if (std::abs(m.ptr[i]) > max) {
                max = std::abs(m.ptr[i]);
            }
        }
        return max;
    }
    
    // level 2
    
    vector<double>& dgemv(mat::trans tr, double alpha, const csrmatrix<double>& A,
     const vector<double>& x, double beta, vector<double>& y) {
        switch (tr) {
        case mat::N:
            if (A.nrow != y.size || A.ncol != x.size) {
                throw;
            }
            spblas_dcsrmv('N', A.nrow, A.ncol, alpha, A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz,
                       x.ptr, x.inc, beta, y.ptr, y.inc);
            break;
        case mat::T:
            if (A.nrow != x.size || A.ncol != y.size) {
                throw;
            }
            spblas_dcsrmv('T', A.nrow, A.ncol, alpha, A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz,
                       x.ptr, x.inc, beta, y.ptr, y.inc);
            break;
        }
        return y;
    }

    csrmatrix<double>& dger(double alpha, const vector<double>& x, const vector<double>& y,
        csrmatrix<double>& A) {
        if (x.size != A.nrow || y.size != A.ncol) {
            throw;
        }
        spblas_dcsrr('N', A.nrow, A.ncol, alpha, x.ptr, x.inc, y.ptr, y.inc,
            A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz);
        return A;
    }

    csrmatrix<double>& dger(mat::trans tr, double alpha,
        const vector<double>& x, const vector<double>& y, csrmatrix<double>& A) {
        switch (tr) {
        case mat::N:
            if (x.size != A.nrow || y.size != A.ncol) {
                throw;
            }
            spblas_dcsrr('N', A.nrow, A.ncol, alpha, x.ptr, x.inc, y.ptr, y.inc,
                A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz);
            break;
        case mat::T:
            if (x.size != A.ncol || y.size != A.nrow) {
                throw;
            }
            spblas_dcsrr('T', A.nrow, A.ncol, alpha, x.ptr, x.inc, y.ptr, y.inc,
                A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz);
            break;
        }
        return A;
    }
    
    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB,
        double alpha, const csrmatrix<double>& A, const dmatrix<double>& B,
        double beta, dmatrix<double>& C) {
        if (trA == mat::N && trB == mat::N) {
            if (A.nrow != C.nrow || B.ncol != C.ncol || A.ncol != B.nrow) {
                throw;
            }
            spblas_dcsrmm('N', 'N', C.nrow, C.ncol, A.ncol,
                alpha, A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz,
                B.ptr, B.ld, beta, C.ptr, C.ld);
        } else if (trA == mat::N && trB == mat::T) {
            if (A.nrow != C.nrow || B.nrow != C.ncol || A.ncol != B.ncol) {
                throw;
            }
            spblas_dcsrmm('N', 'T', C.nrow, C.ncol, A.ncol,
                alpha, A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz,
                B.ptr, B.ld, beta, C.ptr, C.ld);
        } else if (trA == mat::T && trB == mat::N) {
            if (A.ncol != C.nrow || B.ncol != C.ncol || A.nrow != B.nrow) {
                throw;
            }
            spblas_dcsrmm('T', 'N', C.nrow, C.ncol, A.nrow,
                alpha, A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz,
                B.ptr, B.ld, beta, C.ptr, C.ld);
        } else if (trA == mat::T && trB == mat::T) {
            if (A.ncol != C.nrow || B.nrow != C.ncol || A.nrow != B.ncol) {
                throw;
            }
            spblas_dcsrmm('T', 'T', C.nrow, C.ncol, A.nrow,
                alpha, A.ptr, A.rowptr.ptr, A.colind.ptr, A.nnz,
                B.ptr, B.ld, beta, C.ptr, C.ld);
        } else {
            throw;
        }
        return C;
    }

    // vector<double> dkronmv(blas::trans tr, int n1, int n2,
    //     double alpha, const dmatrix_base& A,
    //     const vector<double>& x, double beta, vector<double> y) {
    //     switch (tr) {
    //     case blas::N:
    //         if (A.nrow*n1*n2 != y.size || A.ncol*n1*n2 != x.size || x.inc != 1 || y.inc != 1) {
    //             throw;
    //         }
    //         blas_dkronmv('N', A.nrow, A.ncol, n1, n2,
    //             alpha, A.ptr, A.ld, x.ptr, beta, y.ptr);
    //         break;
    //     case blas::T:
    //         if (A.nrow*n1*n2 != x.size || A.ncol*n1*n2 != y.size || x.inc != 1 || y.inc != 1) {
    //             throw;
    //         }
    //         blas_dkronmv('T', A.nrow, A.ncol, n1, n2, alpha, A.ptr, A.ld,
    //                    x.ptr, beta, y.ptr);
    //         break;
    //     }
    //     return y;
    // }

    // lapack?

    vector<double>& dgsstep(mat::trans tr,
        double alpha, const csrmatrix<double>& A,
        const vector<double>& x, vector<double>& y) {
        /* gauss seidal */
        if (A.nrow != A.ncol) {
            throw;
        }

        int n = A.nrow;
        sci::array<int>::const_range rowptr = A.get_rowptr();
        sci::array<int>::const_range colind = A.get_colind();
        sci::array<double>::const_range value = A.get_value();

        switch (tr) {
        case mat::N: {
            double tmp;
            for (size_t i=1; i<=A.nrow; i++) {
                y(i) = x(i);
                for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
                    if (i != colind[z]) {
                        y(i) -= alpha * value[z] * y(colind[z]);
                    } else {
                        tmp = alpha * value[z];
                    }
                }
                y(i) /= tmp;
            }
            return y;
        }
        case mat::T: {
            vector<int> diag(A.nrow);
            vector<double> tmp(x);
            for (size_t i=1; i<=A.nrow; i++) {
                for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
                    if (i > colind[z]) {
                        tmp(colind[z]) -= alpha * value[z] * y(i);
                    } else if (i == colind[z]) {
                        diag(i) = z;
                        break;
                    } else {
                        throw;
                    }
                }
            }
            for (size_t i=1; i<=A.nrow; i++) {
                y(i) = tmp(i) / (alpha * value[diag(i)]);
                for (int z=diag(i)+1; z<rowptr[i+1]; z++) {
                    tmp(colind[z]) -= alpha * value[z] * y(i);
                }
            }
            return y;
        }
        default:
            throw;
        }
    }

    int dgssolve(mat::trans tr,
        double alpha, const csrmatrix<double>& A,
        const vector<double>& x, vector<double>& y, int maxiter, double eps) {
        int n = x.size;
        sci::vector<double> prevy(n);
        for (int k=1; k<=maxiter; k++) {
            prevy = y;
            sci::dgsstep(tr, alpha, A, x, y);
            if (sci::dnrm2(y-prevy)/sci::dnrm2(y) < eps) {
                return 0;
            }
        }
//        std::cout << "Did not converge Gauss-Seidel algorithm." << std::endl;
        return -1;
    }

    int darnoldi(mat::trans tr, int m,
        const csrmatrix<double>& spA, const vector<double>& x,
        dmatrix<double>& H, dmatrix<double>& V,
        double& beta, int ite, double tol) {

        if (spA.nrow != spA.ncol) {
            throw;
        }

        int info;
        int n = spA.nrow;
        vector<double> w(n);
        spblas_dcsrar(tr, n, m, spA.ptr, spA.rowptr.ptr, spA.colind.ptr, spA.nnz,
            x.ptr, x.inc, H.ptr, H.ld, V.ptr, V.ld, 
            ite, &beta, tol, w.ptr, &info);
        return info;  // return effective dimension
    }

}
