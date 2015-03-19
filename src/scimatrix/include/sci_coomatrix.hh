/*
 Scientific sparse matrix package
 */

#pragma once

#include "sci_spmatrix.hh"

 namespace sci {

    /* Sparse Matrix Base */
    template <class T>
    class coomatrix : public spmatrix<T> {
    public:
        array<int> rowind;
        array<int> colind;

        coomatrix(size_t r = 1, size_t c = 1, size_t n = 1)
        : spmatrix<T>(n), rowind(n), colind(n) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }
        
        coomatrix(size_t r, size_t c, size_t n, int* rp, int* cp, T* v)
        : spmatrix<T>(n,v), rowind(n,rp), colind(n,cp) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }

        coomatrix(size_t r, size_t c, size_t n, T* v)
        : spmatrix<T>(n,v), rowind(n), colind(n) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            spmatrix<T>::nnz = n;
        }

        coomatrix(const coomatrix<T>& m)
        : spmatrix<T>(m), rowind(m.rowind), colind(m.colind) {
            matrix<T>::nrow = m.nrow;
            matrix<T>::ncol = m.ncol;
            spmatrix<T>::nnz = m.nnz;            
        }

        coomatrix(const dmatrix<T>& m)
        : spmatrix<T>(m.nnz()), rowind(m.nnz()), colind(m.nnz()) {
            matrix<T>::nrow = m.nrow;
            matrix<T>::ncol = m.ncol;
            spmatrix<T>::nnz = m.nnz();
            operator=(m);
        }

        virtual ~coomatrix() {
        }

        virtual int type() const {
            return COO;
        }
        
        array<int>::range get_rowind() {
            return rowind.alias(sci::range(1,static_cast<int>(spmatrix<T>::nnz)));
        }
        
        array<int>::range get_colind() {
            return colind.alias(sci::range(1,static_cast<int>(spmatrix<T>::nnz)));
        }
        
        array<int>::const_range get_rowind() const {
            return rowind.alias(sci::range(1,static_cast<int>(spmatrix<T>::nnz)));
        }
        
        array<int>::const_range get_colind() const {
            return colind.alias(sci::range(1,static_cast<int>(spmatrix<T>::nnz)));
        }
        
        // equal
        coomatrix<T>& operator=(const T& c) {
            array<T>::operator=(c);
            return *this;
        }
        
        coomatrix<T>& operator=(const T* a) {
            array<T>::operator=(a);
            return *this;
        }
        
        coomatrix<T>& operator=(const coomatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            rowind = m.rowind.ptr; // use operator=(const T*)
            colind = m.colind.ptr;
            array<T>::operator=(m.ptr);
            return *this;
        }
        
        coomatrix<T>& operator=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }

            array<int>::range rowind = get_rowind();
            array<int>::range colind = get_colind();
            typename array<T>::range value = spmatrix<T>::get_value();

            int z = 1;
            for (size_t j=1; j<=m.ncol; j++) {
                for (int i=1; i<=m.nrow; i++) {
                    if (m.nonzero_check(m(i,j))) {
                        rowind[z] = i;
                        colind[z] = j;
                        value[z] = m(i,j);
                        z++;
                    }
                }
            }
            return *this;
        }

        // operator
        coomatrix<T> operator+(const coomatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            coomatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result += m;
            return result;
        }
        
        coomatrix<T> operator-(const coomatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            coomatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result -= m;
            return result;
        }
        
        coomatrix<T> operator*(const coomatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            coomatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result *= m;
            return result;
        }
        
        coomatrix<T> operator/(const coomatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            coomatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result /= m;
            return result;
        }
        
        coomatrix<T> operator%(const coomatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            coomatrix<T> result(matrix<T>::nrow, matrix<T>::ncol, spmatrix<T>::nnz);
            result = *this;
            result %= m;
            return result;
        }
        
        coomatrix<T>& operator+=(const coomatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator+=(m);
            return *this;
        }
        
        coomatrix<T>& operator-=(const coomatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator-=(m);
            return *this;
        }
        
        coomatrix<T>& operator*=(const coomatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator*=(m);
            return *this;
        }
        
        coomatrix<T>& operator/=(const coomatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator/=(m);
            return *this;
        }
        
        coomatrix<T>& operator%=(const coomatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol || spmatrix<T>::nnz != m.nnz) {
                throw;
            }
            array<T>::operator%=(m);
            return *this;
        }
        
        coomatrix<T>& operator+=(const T& v) {
            array<T>::operator+=(v);
            return *this;
        }

        coomatrix<T>& operator-=(const T& v) {
            array<T>::operator-=(v);
            return *this;
        }

        coomatrix<T>& operator*=(const T& v) {
            array<T>::operator*=(v);
            return *this;
        }

        coomatrix<T>& operator/=(const T& v) {
            array<T>::operator/=(v);
            return *this;
        }

        coomatrix<T>& operator%=(const T& v) {
            array<T>::operator%=(v);
            return *this;
        }

        // print
        
        virtual std::ostream& print(std::ostream& os) const {
            array<int>::const_range r_rowind = get_rowind();
            array<int>::const_range r_colind = get_colind();
            typename array<T>::const_range r_value = spmatrix<T>::get_value();
            
            for (size_t z=1; z<=spmatrix<T>::nnz; z++) {
                os << "(" << r_rowind[z] << "," << r_colind[z] << ")="
                << r_value[z] << std::endl;
            }
            return os;
        }
    };

    template<class T>
    dmatrix<T>::dmatrix(const coomatrix<T>& m) : matrix<T>(m.nrow * m.ncol) {
        matrix<T>::nrow = m.nrow;
        matrix<T>::ncol = m.ncol;
        ld = m.nrow;
        operator=(m);
    }

    template<class T>
    dmatrix<T>& dmatrix<T>::operator=(const coomatrix<T>& m) {
        if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
            throw;
        }

        array<int>::const_range rowind = m.get_rowind();
        array<int>::const_range colind = m.get_colind();
        typename array<T>::const_range value = m.get_value();

        operator=(zero());

        for (size_t z=1; z<=m.nnz; z++) {
            operator()(rowind[z], colind[z]) = value[z];
        }
        return *this;
    }

}

