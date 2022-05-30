/*
 sci blas
 */

#include <cmath>

#include "blas.h"
#include "sci_vector.h"

namespace sci {

    // level 1
    
    double* dcopy(const vector<double>& x, double *p) {
        blas_dcopy(x.size, x.ptr, x.inc, p, 1);
        return p;
    }
    
    vector<double>& dcopy(const double *p, vector<double>& y) {
        blas_dcopy(y.size, p, 1, y.ptr, y.inc);
        return y;
    }
    
    vector<double>& dcopy(const vector<double>& x, vector<double>& y) {
        if (x.size != y.size) {
            throw;
        }
        blas_dcopy(x.size, x.ptr, x.inc, y.ptr, y.inc);
        return y;
    }
    
    vector<double>& dfill(double v, vector<double>& x) {
        blas_dfill(x.size, v, x.ptr, x.inc);
        return x;
    }
    
    vector<double>& dscal(double alpha, vector<double>& x) {
        blas_dscal(x.size, alpha, x.ptr, x.inc);
        return x;
    }
    
    vector<double>& daxpy(double alpha, const vector<double>& x, vector<double>& y) {
        if (x.size != y.size) {
            throw;
        }
        blas_daxpy(x.size, alpha, x.ptr, x.inc, y.ptr, y.inc);
        return y;
    }
    
    double dsum(const vector<double>& x) {
        return blas_dsum(x.size, x.ptr, x.inc);
    }
    
    double dasum(const vector<double>& x) {
        return blas_dasum(x.size, x.ptr, x.inc);
    }
    
    double dnrm2(const vector<double>& x) {
        return blas_dnrm2(x.size, x.ptr, x.inc);
    }
    
    int idmax(const vector<double>& x) {
        return blas_idmax(x.size, x.ptr, x.inc) + 1;
    }
    
    int idamax(const vector<double>& x) {
        return blas_idamax(x.size, x.ptr, x.inc) + 1;
    }
    
    double dmax(const vector<double>& x) {
        return x(idmax(x));
    }
    
    double damax(const vector<double>& x) {
        return std::abs(x(idamax(x)));
    }
    
    double ddot(const vector<double>& x, const vector<double>& y) {
        if (x.size != y.size) {
            throw;
        }
        return blas_ddot(x.size, x.ptr, x.inc, y.ptr, y.inc);
    }
   
}

