#include <cmath>
#include <cfloat>

#include "blas.h"
#include "lapack.h"
#include "sci_blas.h"

namespace sci {

    static double dmatrix_double_zero = 0.0;

    template<>
    bool dmatrix<double>::nonzero_check(const double& v) const {
        if (std::abs(v) < DBL_EPSILON) {
            return false;
        } else {
            return true;
        }
    }

    template<>
    const double& dmatrix<double>::zero() const {
        return dmatrix_double_zero;
    }

    vector<double> diag(const dmatrix<double>& m) {
        if (m.nrow != m.ncol) {
            throw;
        }
        vector<double> result(m.nrow);
        for (size_t i=1; i<=m.nrow; i++) {
            result(i) = m(i,i);
        }
        return result;
    }
    
    double* dcopy(const dmatrix<double>& m, double *p) {
        if (m.nrow != m.ld) {
            for (size_t jx=1; jx<=m.ncol; jx++, p+=m.nrow) {
                blas_dcopy(m.nrow, &m(1, jx), 1, p, 1);
            }
        } else {
            blas_dcopy(m.nrow * m.ncol, &m(1,1), 1, p, 1);
        }
        return p;
    }

    dmatrix<double>& dcopy(const double *p, dmatrix<double>& m) {
        if (m.nrow != m.ld) {
            for (size_t jx=1; jx<=m.ncol; jx++, p+=m.nrow) {
                blas_dcopy(m.nrow, p, 1, &m(1,jx), 1);
            }
        } else {
            blas_dcopy(m.nrow * m.ncol, p, 1, &m(1,1), 1);
        }
        return m;
    }

    dmatrix<double>& dcopy(mat::trans tr, const double *p, dmatrix<double>& m) {
        if (tr == mat::N) {
            return dcopy(p, m);
        } else if (tr == mat::T) {
            for (size_t ix=1; ix<=m.nrow; ix++, p+=m.ncol) {
                blas_dcopy(m.ncol, p, 1, &m(ix,1), m.ld);
            }
            return m;
        } else {
            throw;
        }
    }

    // dcopy4

    dmatrix<double>& dcopy(const dmatrix<double>& x, dmatrix<double>& y) {
        if (x.nrow != y.nrow || x.ncol != y.ncol) {
            throw;
        }
        if (x.ld != x.nrow || y.ld != y.nrow) {
            for (size_t jx=1, jy=1; jx<=x.ncol; jx++, jy++) {
                blas_dcopy(x.nrow, &x(1,jx), 1, &y(1,jy), 1);
            }
        } else {
            blas_dcopy(x.nrow * x.ncol, &x(1,1), 1, &y(1,1), 1);
        }
        return y;
    }

    // dcopy5

    dmatrix<double>& dcopy(mat::trans tr, const dmatrix<double>& x, dmatrix<double>& y) {
        if (tr == mat::N) {
            return dcopy(x, y);
        } else if (tr == mat::T) {
            if (x.nrow != y.ncol || x.ncol != y.nrow) {
                throw;
            }
            for (size_t ix=1; ix<=x.ncol; ix++) {
                blas_dcopy(x.nrow, &x(1,ix), 1, &y(ix,1), y.ld);
            }
            return y;
        } else {
            throw;
        }
    }
    
    // dfill

    dmatrix<double>& dfill(double v, dmatrix<double>& m) {
        if (m.ld != m.nrow) {
            for (size_t jx=1; jx<=m.ncol; jx++) {
                blas_dfill(m.nrow, v, &m(1,jx), 1);
            }
        } else {
            blas_dfill(m.nrow * m.ncol, v, &m(1,1), 1);
        }
        return m;
    }

    // dscal

    dmatrix<double>& dscal(double alpha, dmatrix<double>& m) {
        if (m.ld != m.nrow) {
            for (size_t jx=1; jx<=m.ncol; jx++) {
                blas_dscal(m.nrow, alpha, &m(1,jx), 1);
            }
        } else {
            blas_dscal(m.nrow * m.ncol, alpha, &m(1,1), 1);
        }
        return m;
    }

    // daxpy

    dmatrix<double>& daxpy(double alpha, const dmatrix<double>& x, dmatrix<double>& y) {
        if (x.nrow != y.nrow || x.ncol != y.ncol) {
            throw;
        }
        if (x.ld != x.nrow || y.ld != y.nrow) {
            for (size_t jx=1, jy=1; jx<=x.ncol; jx++, jy++) {
                blas_daxpy(x.nrow, alpha, &x(1,jx), 1, &y(1,jy), 1);
            }
        } else {
            blas_daxpy(x.nrow * x.ncol, alpha, &x(1,1), 1, &y(1,1), 1);
        }
        return y;
    }

    // dnrm2

    double dnrm2(const dmatrix<double>& m) {
        double sum = 0.0;
        for (size_t j=1; j<=m.ncol; j++) {
            for (size_t i=1; i<=m.nrow; i++) {
                sum += m(i,j) * m(i,j);
            }
        }
        return sqrt(sum);
    }
    
    // damax

    double damax(const dmatrix<double>& m) {
        double max = 0.0;
        for (size_t j=1; j<=m.ncol; j++) {
            for (size_t i=1; i<=m.nrow; i++) {
                if (std::abs(m(i,j)) > max) {
                    max = std::abs(m(i,j));
                }
            }
        }
        return max;
    }
    
    // level 2
    
    // dgemv

    vector<double>& dgemv(mat::trans tr, double alpha, const dmatrix<double>& A,
     const vector<double>& x, double beta, vector<double>& y) {
        switch (tr) {
        case mat::N:
            if (A.nrow != y.size || A.ncol != x.size) {
                throw;
            }
            blas_dgemv('N', A.nrow, A.ncol, alpha, A.ptr, A.ld,
                       x.ptr, x.inc, beta, y.ptr, y.inc);
            break;
        case mat::T:
            if (A.nrow != x.size || A.ncol != y.size) {
                throw;
            }
            blas_dgemv('T', A.nrow, A.ncol, alpha, A.ptr, A.ld,
                       x.ptr, x.inc, beta, y.ptr, y.inc);
            break;
        }
        return y;
    }

    // dger
    
    dmatrix<double>& dger(double alpha, const vector<double>& x,
        const vector<double>& y, dmatrix<double>& A) {
        if (x.size != A.nrow || y.size != A.ncol) {
            throw;
        }
        blas_dger(A.nrow, A.ncol, alpha, x.ptr, x.inc, y.ptr, y.inc, A.ptr, A.ld);
        return A;
    }

    // dger2
    
    dmatrix<double>& dger(mat::trans tr, double alpha,
        const vector<double>& x, const vector<double>& y, dmatrix<double>& A) {
        switch (tr) {
        case mat::N:
            return dger(alpha, x, y, A);
            break;
        case mat::T:
            if (x.size != A.ncol || y.size != A.nrow) {
                throw;
            }
            blas_dger(A.nrow, A.ncol, alpha, y.ptr, y.inc, x.ptr, x.inc, A.ptr, A.ld);
            break;
        }
        return A;
    }
    
    dmatrix<double>& dgemm(mat::trans trA, mat::trans trB,
        double alpha, const dmatrix<double>& A, const dmatrix<double>& B,
        double beta, dmatrix<double>& C) {
        if (trA == mat::N && trB == mat::N) {
            if (A.nrow != C.nrow || B.ncol != C.ncol || A.ncol != B.nrow) {
                throw;
            }
            blas_dgemm('N', 'N', C.nrow, C.ncol, A.ncol, alpha,
                       A.ptr, A.ld, B.ptr, B.ld,
                       beta, C.ptr, C.ld);
        } else if (trA == mat::N && trB == mat::T) {
            if (A.nrow != C.nrow || B.nrow != C.ncol || A.ncol != B.ncol) {
                throw;
            }
            blas_dgemm('N', 'T', C.nrow, C.ncol, A.ncol, alpha,
                       A.ptr, A.ld, B.ptr, B.ld,
                       beta, C.ptr, C.ld);
        } else if (trA == mat::T && trB == mat::N) {
            if (A.ncol != C.nrow || B.ncol != C.ncol || A.nrow != B.nrow) {
                throw;
            }
            blas_dgemm('T', 'N', C.nrow, C.ncol, A.nrow, alpha,
                       A.ptr, A.ld, B.ptr, B.ld,
                       beta, C.ptr, C.ld);
        } else if (trA == mat::T && trB == mat::T) {
            if (A.ncol != C.nrow || B.nrow != C.ncol || A.nrow != B.ncol) {
                throw;
            }
            blas_dgemm('T', 'T', C.nrow, C.ncol, A.nrow, alpha,
                       A.ptr, A.ld, B.ptr, B.ld,
                       beta, C.ptr, C.ld);
        } else {
            throw;
        }
        return C;
    }

    vector<double>& dkronmv(mat::trans tr, int n1, int n2,
        double alpha, const dmatrix<double>& A,
        const vector<double>& x, double beta, vector<double>& y) {
        switch (tr) {
        case mat::N:
            if (A.nrow*n1*n2 != y.size || A.ncol*n1*n2 != x.size || x.inc != 1 || y.inc != 1) {
                throw;
            }
            blas_dkronmv('N', A.nrow, A.ncol, n1, n2,
                alpha, A.ptr, A.ld, x.ptr, beta, y.ptr);
            break;
        case mat::T:
            if (A.nrow*n1*n2 != x.size || A.ncol*n1*n2 != y.size || x.inc != 1 || y.inc != 1) {
                throw;
            }
            blas_dkronmv('T', A.nrow, A.ncol, n1, n2, alpha, A.ptr, A.ld,
                       x.ptr, beta, y.ptr);
            break;
        }
        return y;
    }

    ///// lapack

    int dgesv(mat::trans tr, double alpha, const dmatrix<double>& A,
        const vector<double>& B, vector<double>& C) {
        if (A.nrow != A.ncol) {
          throw;
        }

        int info;

        vector<int> ipiv(A.nrow);
        dmatrix<double> TA(A.nrow, A.ncol);
        dcopy(tr, A, TA);
        dscal(alpha, TA);
        C = B;

        info = blas_dgesv(TA.nrow, 1, TA.ptr, TA.ld, ipiv.ptr, C.ptr, C.size);
        return info;
    }


    int dgesv(const dmatrix<double>& A, const dmatrix<double>& B, dmatrix<double>& C) {
        if (A.nrow != A.ncol) {
          throw;
        }

        int info;

        vector<int> ipiv(A.nrow);
        dmatrix<double> TA(A);
        C = B;

        info = blas_dgesv(TA.nrow, C.ncol, TA.ptr, TA.ld,
            ipiv.ptr, C.ptr, C.ld);
        return info;
    }

    dmatrix<double>& dpotrf(char uplo, dmatrix<double>& A) {
        blas_dpotrf(uplo, A.nrow, A.ptr, A.ld);
        return A;
    }

    double dpodet(const dmatrix<double>& A) {
        dmatrix<double> tmpA(A);
        dpotrf('U', tmpA);
        double logdet = 0.0;
        for (size_t i=1; i<=tmpA.nrow; i++) {
            logdet += 2.0 * log(tmpA(i,i));
        }
        return logdet;
    }

    int dggglm(const dmatrix<double>& X, const dmatrix<double>& W,
        const vector<double>& z, vector<double>& beta, vector<double>& residual,
        int lwork, double *work) {
        int info;
        dmatrix<double> m_x(X);
        dmatrix<double> m_w(W);
        vector<double> m_z(z);

        info = blas_dggglm(m_x.nrow, m_x.ncol, m_w.ncol,
            m_x.ptr, m_x.nrow, m_w.ptr, m_w.nrow, 
            m_z.ptr, beta.ptr, residual.ptr, work, lwork);

        return info;
    }

}
