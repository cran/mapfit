/*
 Scientific sparse matrix package
 */

#pragma once

#include "sci_spmatrix.h"

 namespace sci {

    /* Sparse Matrix Base */
    template <class T>
    class csrmatrix : public spmatrix<T> {
    public:
        array<int> rowptr;
        array<int> colind;

        csrmatrix(size_t r = 1, size_t c = 1, size_t n = 1)
        : spmatrix<T>(n), rowptr(r+1), colind(n) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }
        
        csrmatrix(size_t r, size_t c, size_t n, int* rp, int* cp, T* v)
        : spmatrix<T>(n,v), rowptr(r+1,rp), colind(n,cp) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }

        csrmatrix(size_t r, size_t c, size_t n, T* v)
        : spmatrix<T>(n,v), rowptr(r+1), colind(n) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }

        csrmatrix(const csrmatrix<T>& m)
        : spmatrix<T>(m), rowptr(m.rowptr), colind(m.colind) {
            matrix<T>::nrow = m.nrow;
            matrix<T>::ncol = m.ncol;
            spmatrix<T>::nnz = m.nnz;            
        }

        csrmatrix(const dmatrix<T>& m)
        : spmatrix<T>(m.nnz()), rowptr(m.nrow+1), colind(m.nnz()) {
            matrix<T>::nrow = m.nrow;
            matrix<T>::ncol = m.ncol;
            spmatrix<T>::nnz = m.nnz();
            operator=(m);
        }

        virtual ~csrmatrix() {
        }

        virtual int type() const {
            return CSR;
        }
        
        array<int>::range get_rowptr() {
            return rowptr.alias(sci::range(1,static_cast<int>(matrix<T>::nrow+1)));
        }
        
        array<int>::range get_colind() {
            return colind.alias(sci::range(1,static_cast<int>(spmatrix<T>::nnz)));
        }
        
        array<int>::const_range get_rowptr() const {
            return rowptr.alias(sci::range(1,static_cast<int>(matrix<T>::nrow+1)));
        }
        
        array<int>::const_range get_colind() const {
            return colind.alias(sci::range(1,static_cast<int>(spmatrix<T>::nnz)));
        }
        
        // equal
        csrmatrix<T>& operator=(const T& c) {
            array<T>::operator=(c);
            return *this;
        }
        
        csrmatrix<T>& operator=(const T* a) {
            array<T>::operator=(a);
            return *this;
        }
        
        csrmatrix<T>& operator=(const csrmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            rowptr = m.rowptr.ptr; // use operator=(const T*)
            colind = m.colind.ptr;
            array<T>::operator=(m.ptr);
            return *this;
        }
        
        csrmatrix<T>& operator=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }

            array<int>::range rowptr = get_rowptr();
            array<int>::range colind = get_colind();
            typename array<T>::range value = spmatrix<T>::get_value();

            int z = 1;
            for (size_t i=1; i<=m.nrow; i++) {
                rowptr[i] = z;
                for (size_t j=1; j<=m.ncol; j++) {
                    if (m.nonzero_check(m(i,j))) {
                        colind[z] = j;
                        value[z] = m(i,j);
                        z++;
                    }
                }
            }
            rowptr[m.nrow+1] = z;
            return *this;
        }

        // operator
        csrmatrix<T> operator+(const csrmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            csrmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result += m;
            return result;
        }
        
        csrmatrix<T> operator-(const csrmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            csrmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result -= m;
            return result;
        }
        
        csrmatrix<T> operator*(const csrmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            csrmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result *= m;
            return result;
        }
        
        csrmatrix<T> operator/(const csrmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            csrmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result /= m;
            return result;
        }
        
        csrmatrix<T> operator%(const csrmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            csrmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result %= m;
            return result;
        }
        
        csrmatrix<T>& operator+=(const csrmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator+=(m);
            return *this;
        }
        
        csrmatrix<T>& operator-=(const csrmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator-=(m);
            return *this;
        }
        
        csrmatrix<T>& operator*=(const csrmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator*=(m);
            return *this;
        }
        
        csrmatrix<T>& operator/=(const csrmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator/=(m);
            return *this;
        }
        
        csrmatrix<T>& operator%=(const csrmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator%=(m);
            return *this;
        }
        
        csrmatrix<T>& operator+=(const T& v) {
            array<T>::operator+=(v);
            return *this;
        }

        csrmatrix<T>& operator-=(const T& v) {
            array<T>::operator-=(v);
            return *this;
        }

        csrmatrix<T>& operator*=(const T& v) {
            array<T>::operator*=(v);
            return *this;
        }

        csrmatrix<T>& operator/=(const T& v) {
            array<T>::operator/=(v);
            return *this;
        }

        csrmatrix<T>& operator%=(const T& v) {
            array<T>::operator%=(v);
            return *this;
        }

        // print
        
        virtual std::ostream& print(std::ostream& os) const {
            array<int>::const_range rowptr = get_rowptr();
            array<int>::const_range colind = get_colind();
            typename array<T>::const_range value = spmatrix<T>::get_value();
            
            for (size_t i=1; i<=matrix<T>::nrow; i++) {
                for (int z=rowptr[i]; z<rowptr[i+1]; z++) {
                    os << "(" << i << "," << colind[z] << ")="
                    << value[z] << std::endl;
                }
            }
            return os;
        }
    };

    template<class T>
    dmatrix<T>::dmatrix(const csrmatrix<T>& m) : matrix<T>(m.nrow * m.ncol) {
        matrix<T>::nrow = m.nrow;
        matrix<T>::ncol = m.ncol;
        ld = m.nrow;
        operator=(m);
    }

    template<class T>
    dmatrix<T>& dmatrix<T>::operator=(const csrmatrix<T>& m) {
        if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
            throw;
        }

        array<int>::const_range rowptr = m.get_rowptr();
        array<int>::const_range colind = m.get_colind();
        typename array<T>::const_range value = m.get_value();

        operator=(zero());

        for (size_t i=1; i<=m.nrow; i++) {
            for (int j=rowptr[i]; j<rowptr[i+1]; j++) {
                operator()(i, colind[j]) = value[j];
            }
        }
        return *this;
    }
}

