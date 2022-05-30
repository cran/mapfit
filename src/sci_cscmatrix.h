/*
 Scientific sparse matrix package
 */

#pragma once

#include "sci_spmatrix.h"

 namespace sci {

    /* Sparse Matrix Base */
    template <class T>
    class cscmatrix : public spmatrix<T> {
    public:
        array<int> colptr;
        array<int> rowind;

        cscmatrix(size_t r = 1, size_t c = 1, size_t n = 1)
        : spmatrix<T>(n), colptr(c+1), rowind(n) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }
        
        cscmatrix(size_t r, size_t c, size_t n, int* cp, int* rp, T* v)
        : spmatrix<T>(n,v), colptr(c+1,cp), rowind(n,rp) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }

        cscmatrix(size_t r, size_t c, size_t n, T* v)
        : spmatrix<T>(n,v), colptr(c+1), rowind(n) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }

        cscmatrix(const cscmatrix<T>& m)
        : spmatrix<T>(m), colptr(m.colptr), rowind(m.rowind) {
            matrix<T>::nrow = m.nrow;
            matrix<T>::ncol = m.ncol;
            spmatrix<T>::nnz = m.nnz;            
        }

        cscmatrix(const dmatrix<T>& m)
        : spmatrix<T>(m.nnz()), colptr(m.ncol+1), rowind(m.nnz()) {
            matrix<T>::nrow = m.nrow;
            matrix<T>::ncol = m.ncol;
            spmatrix<T>::nnz = m.nnz();
            operator=(m);
        }

        virtual ~cscmatrix() {
        }

        virtual int type() const {
            return CSC;
        }
        
        array<int>::range get_colptr() {
            return colptr.alias(sci::range(1,static_cast<int>(matrix<T>::ncol+1)));
        }
        
        array<int>::range get_rowind() {
            return rowind.alias(sci::range(1,static_cast<int>(spmatrix<T>::nnz)));
        }
        
        array<int>::const_range get_colptr() const {
            return colptr.alias(sci::range(1,static_cast<int>(matrix<T>::ncol+1)));
        }
        
        array<int>::const_range get_rowind() const {
            return rowind.alias(sci::range(1,static_cast<int>(spmatrix<T>::nnz)));
        }
        
        // equal
        cscmatrix<T>& operator=(const T& c) {
            array<T>::operator=(c);
            return *this;
        }
        
        cscmatrix<T>& operator=(const T* a) {
            array<T>::operator=(a);
            return *this;
        }
        
        cscmatrix<T>& operator=(const cscmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            colptr = m.colptr.ptr; // use operator=(const T*)
            rowind = m.rowind.ptr;
            array<T>::operator=(m.ptr);
            return *this;
        }
        
        cscmatrix<T>& operator=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }

            array<int>::range colptr = get_colptr();
            array<int>::range rowind = get_rowind();
            typename array<T>::range value = spmatrix<T>::get_value();

            int z = 1;
            for (size_t i=1; i<=m.ncol; i++) {
                colptr[i] = z;
                for (size_t j=1; j<=m.nrow; j++) {
                    if (m.nonzero_check(m(j,i))) {
                        rowind[z] = j;
                        value[z] = m(j,i);
                        z++;
                    }
                }
            }
            colptr[m.ncol+1] = z;
            return *this;
        }

        // operator
        cscmatrix<T> operator+(const cscmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            cscmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result += m;
            return result;
        }
        
        cscmatrix<T> operator-(const cscmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            cscmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result -= m;
            return result;
        }
        
        cscmatrix<T> operator*(const cscmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            cscmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result *= m;
            return result;
        }
        
        cscmatrix<T> operator/(const cscmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            cscmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result /= m;
            return result;
        }
        
        cscmatrix<T> operator%(const cscmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            cscmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result %= m;
            return result;
        }
        
        cscmatrix<T>& operator+=(const cscmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator+=(m);
            return *this;
        }
        
        cscmatrix<T>& operator-=(const cscmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator-=(m);
            return *this;
        }
        
        cscmatrix<T>& operator*=(const cscmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator*=(m);
            return *this;
        }
        
        cscmatrix<T>& operator/=(const cscmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator/=(m);
            return *this;
        }
        
        cscmatrix<T>& operator%=(const cscmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator%=(m);
            return *this;
        }
        
        cscmatrix<T>& operator+=(const T& v) {
            array<T>::operator+=(v);
            return *this;
        }

        cscmatrix<T>& operator-=(const T& v) {
            array<T>::operator-=(v);
            return *this;
        }

        cscmatrix<T>& operator*=(const T& v) {
            array<T>::operator*=(v);
            return *this;
        }

        cscmatrix<T>& operator/=(const T& v) {
            array<T>::operator/=(v);
            return *this;
        }

        cscmatrix<T>& operator%=(const T& v) {
            array<T>::operator%=(v);
            return *this;
        }

        // print
        
        virtual std::ostream& print(std::ostream& os) const {
            array<int>::const_range r_colptr = get_colptr();
            array<int>::const_range r_rowind = get_rowind();
            typename array<T>::const_range r_value = spmatrix<T>::get_value();
            
            for (size_t i=1; i<=matrix<T>::ncol; i++) {
                for (int j=r_colptr[i]; j<r_colptr[i+1]; j++) {
                    os << "(" << r_rowind[j] << "," << i << ")="
                    << r_value[j] << std::endl;
                }
            }
            return os;
        }
    };

    template<class T>
    dmatrix<T>::dmatrix(const cscmatrix<T>& m) : matrix<T>(m.nrow * m.ncol) {
        matrix<T>::nrow = m.nrow;
        matrix<T>::ncol = m.ncol;
        ld = m.nrow;
        operator=(m);
    }

    template<class T>
    dmatrix<T>& dmatrix<T>::operator=(const cscmatrix<T>& m) {
        if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
            throw;
        }

        array<int>::const_range colptr = m.get_colptr();
        array<int>::const_range rowind = m.get_rowind();
        typename array<T>::const_range value = m.get_value();

        operator=(zero());

        for (size_t i=1; i<=m.ncol; i++) {
            for (int j=colptr[i]; j<colptr[i+1]; j++) {
                operator()(rowind[j], i) = value[j];
            }
        }
        return *this;
    }
}

