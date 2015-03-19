#include <cmath>

#include "blas.h"
#include "spblas.h"
#include "arnoldi.h"

#include "sci_spblas.hh"
#include "sci_cscmatrix.hh"

namespace sci {

    vector<double> diag(const cscmatrix<double>& m) {
        if (m.nrow != m.ncol) {
            throw;
        }
        vector<double> result(m.nrow);
        array<int>::const_range colptr = m.get_colptr();
        array<int>::const_range rowind = m.get_rowind();
        array<double>::const_range value = m.get_value();
        
        for (size_t i=1; i<=m.ncol; i++) {
            for (int j=colptr[i]; j<colptr[i+1]; j++) {
                if (i == rowind[j]) {
                    result(i) = value[j];
                    break;
                }
                if (i < rowind[j]) {
                    break;
                }
            }
        }
        return result;
    }
    
    double* dcopy(const cscmatrix<double>& m, double *p) {
        blas_dcopy(m.nnz, m.ptr, 1, p, 1);
        return p;
    }

    cscmatrix<double>& dcopy(const double *p, cscmatrix<double>& m) {
        blas_dcopy(m.nnz, p, 1, m.ptr, 1);
        return m;
    }

    cscmatrix<double>& dfill(double v, cscmatrix<double>& m) {
        blas_dfill(m.nnz, v, m.ptr, 1);
        return m;
    }

    cscmatrix<double>& dscal(double alpha, cscmatrix<double>& m) {
        blas_dscal(m.nnz, alpha, m.ptr, 1);
        return m;
    }

    cscmatrix<double>& daxpy(double alpha, const cscmatrix<double>& x, cscmatrix<double>& y) {
        if (x.nrow != y.nrow || x.ncol != y.ncol || x.nnz != y.nnz) {
            throw;
        }
        blas_daxpy(x.nnz, alpha, x.ptr, 1, y.ptr, 1);
        return y;
    }

    double dnrm2(const cscmatrix<double>& m) {
        double sum = 0.0;
        for (size_t i=0; i<m.nnz; i++) {
            sum += m.ptr[i] * m.ptr[i];
        }
        return sqrt(sum);
    }
    
    double damax(const cscmatrix<double>& m) {
        double max = 0.0;
        for (size_t i=0; i<m.nnz; i++) {
            if (std::abs(m.ptr[i]) > max) {
                max = std::abs(m.ptr[i]);
            }
        }
        return max;
    }
    
    // level 2
    
    vector<double>& dgemv(mat::trans tr, double alpha, const cscmatrix<double>& A,
     const vector<double>& x, double beta, vector<double>& y) {
        switch (tr) {
        case mat::N:
            if (A.nrow != y.size || A.ncol != x.size) {
                throw;
            }
            spblas_dcscmv('N', A.nrow, A.ncol, alpha, A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz,
                       x.ptr, x.inc, beta, y.ptr, y.inc);
            break;
        case mat::T:
            if (A.nrow != x.size || A.ncol != y.size) {
                throw;
            }
            spblas_dcscmv('T', A.nrow, A.ncol, alpha, A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz,
                       x.ptr, x.inc, beta, y.ptr, y.inc);
            break;
        }
        return y;
    }

    cscmatrix<double>& dger(double alpha, const vector<double>& x, const vector<double>& y,
        cscmatrix<double>& A) {
        if (x.size != A.nrow || y.size != A.ncol) {
            throw;
        }
        spblas_dcscr('N', A.nrow, A.ncol, alpha, x.ptr, x.inc, y.ptr, y.inc,
            A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz);
        return A;
    }

    cscmatrix<double>& dger(mat::trans tr, double alpha,
        const vector<double>& x, const vector<double>& y, cscmatrix<double>& A) {
        switch (tr) {
        case mat::N:
            if (x.size != A.nrow || y.size != A.ncol) {
                throw;
            }
            spblas_dcscr('N', A.nrow, A.ncol, alpha, x.ptr, x.inc, y.ptr, y.inc,
                A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz);
            break;
        case mat::T:
            if (x.size != A.ncol || y.size != A.nrow) {
                throw;
            }
            spblas_dcscr('T', A.nrow, A.ncol, alpha, x.ptr, x.inc, y.ptr, y.inc,
                A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz);
            break;
        }
        return A;
    }
    
    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB,
        double alpha, const cscmatrix<double>& A, const dmatrix<double>& B,
        double beta, dmatrix<double>& C) {
        if (trA == mat::N && trB == mat::N) {
            if (A.nrow != C.nrow || B.ncol != C.ncol || A.ncol != B.nrow) {
                throw;
            }
            spblas_dcscmm('N', 'N', C.nrow, C.ncol, A.ncol,
                alpha, A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz,
                B.ptr, B.ld, beta, C.ptr, C.ld);
        } else if (trA == mat::N && trB == mat::T) {
            if (A.nrow != C.nrow || B.nrow != C.ncol || A.ncol != B.ncol) {
                throw;
            }
            spblas_dcscmm('N', 'T', C.nrow, C.ncol, A.ncol,
                alpha, A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz,
                B.ptr, B.ld, beta, C.ptr, C.ld);
        } else if (trA == mat::T && trB == mat::N) {
            if (A.ncol != C.nrow || B.ncol != C.ncol || A.nrow != B.nrow) {
                throw;
            }
            spblas_dcscmm('T', 'N', C.nrow, C.ncol, A.nrow,
                alpha, A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz,
                B.ptr, B.ld, beta, C.ptr, C.ld);
        } else if (trA == mat::T && trB == mat::T) {
            if (A.ncol != C.nrow || B.nrow != C.ncol || A.nrow != B.ncol) {
                throw;
            }
            spblas_dcscmm('T', 'T', C.nrow, C.ncol, A.nrow,
                alpha, A.ptr, A.colptr.ptr, A.rowind.ptr, A.nnz,
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
        double alpha, const cscmatrix<double>& A,
        const vector<double>& x, vector<double>& y) {
        /* gauss seidal */
        if (A.nrow != A.ncol) {
            throw;
        }

        int n = A.nrow;
        sci::array<int>::const_range colptr = A.get_colptr();
        sci::array<int>::const_range rowind = A.get_rowind();
        sci::array<double>::const_range value = A.get_value();

        switch (tr) {
        case mat::N: {
            vector<int> diag(A.ncol);
            vector<double> tmp(x);
            for (size_t j=1; j<=A.ncol; j++) {
                for (int z=colptr[j]; z<colptr[j+1]; z++) {
                    if (j > rowind[z]) {
                        tmp(rowind[z]) -= alpha * value[z] * y(j);
                    } else if (j == rowind[z]) {
                        diag(j) = z;
                        break;
                    } else {
                        throw;
                    }
                }
            }
            for (size_t j=1; j<=A.nrow; j++) {
                y(j) = tmp(j) / (alpha * value[diag(j)]);
                for (int z=diag(j)+1; z<colptr[j+1]; z++) {
                    tmp(rowind[z]) -= alpha * value[z] * y(j);
                }
            }
            return y;
        }
        case mat::T: {
            double tmp;
            for (size_t j=1; j<=A.ncol; j++) {
                y(j) = x(j);
                for (int z=colptr[j]; z<colptr[j+1]; z++) {
                    if (j != rowind[z]) {
                        y(j) -= alpha * value[z] * y(rowind[z]);
                    } else {
                        tmp = alpha * value[z];
                    }
                }
                y(j) /= tmp;
            }
            return y;
        }
        default:
            throw;
        }
    }

    int dgssolve(mat::trans tr,
        double alpha, const cscmatrix<double>& A,
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
        const cscmatrix<double>& spA, const vector<double>& x,
        dmatrix<double>& H, dmatrix<double>& V,
        double& beta, int ite, double tol) {

        if (spA.nrow != spA.ncol) {
            throw;
        }

        int info;
        int n = spA.nrow;
        vector<double> w(n);
        spblas_dcscar(tr, n, m, spA.ptr, spA.colptr.ptr, spA.rowind.ptr, spA.nnz,
            x.ptr, x.inc, H.ptr, H.ld, V.ptr, V.ld, 
            ite, &beta, tol, w.ptr, &info);
        return info;  // return effective dimension
    }
}
