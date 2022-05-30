#pragma once

#include "sci_matrix.h"

namespace sci {

    template<class T>
    vector<double>& dcgnr(mat::trans tr, const T& A,
        const vector<double>& b, const vector<double>& x,
        vector<double>& y, double tol) {

        double alpha0, alpha1, alpha;
        mat::trans ttr = mat::t(tr);
        vector<double> p(x.size);
        vector<double> z(x.size);
        vector<double> r(b.size);
        vector<double> w(b.size);

        // y = x
        y = x;

        // r = b - A * y
        r = b;
        dgemv(tr, -1.0, A, y, 1.0, r);

        // z = A^T r
        dgemv(ttr, 1.0, A, r, 0.0, z);
        p = z;

        alpha0 = ddot(z, z);

        for (int k=1; k<=x.size; k++) {

            // w = A * p
            dgemv(tr, 1.0, A, p, 0.0, w);

            // alpha = (r,r)/(r,A*r)
            alpha =  alpha0 / ddot(w, w);

            // y = y + alpha * p
            daxpy(alpha, p, y);

            // r = r - alpha * w
            daxpy(-alpha, w, r);

            if (dnrm2(r) < tol) {
                return y;
            }

            // z = A * r
            dgemv(ttr, 1.0, A,  r, 0.0, z);

            alpha1 = ddot(z, z);

            // p = z + beta * p
            dscal(alpha1 / alpha0, p);
            p += z;

            alpha0 = alpha1;
        }

        return y;
    }
}
