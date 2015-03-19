/*
 Scientific vector and matrix package
 */

#pragma once

#include "sci_vector.hh"

// To distinguish a combination of two matrices, we assign prime numbers.
#define BASEP 5

#define DENSE 1
#define CSR 2
#define CSC 3
#define COO 4

#define PAIR(a,b) (BASEP*(a)+(b))

namespace sci {
    
    template<class T> class dmatrix;
    template<class T> class csrmatrix;
    template<class T> class cscmatrix;
    template<class T> class coomatrix;

    template<class T>
    class matrix : public array<T> {
    public:
        size_t nrow;
        size_t ncol;

        matrix(size_t s) : array<T>(s) {}
        matrix(const array<T>& a) : array<T>(a) {}
        matrix(size_t s, T* v) : array<T>(s, v) {}
        matrix(size_t s, const T& v) : array<T>(s, v) {}
        matrix(size_t s, T* p, std::vector<T>* e) : array<T>(s, p, e) {}

        virtual ~matrix() {
        }

        virtual int type() const = 0;
    };

//    template <class T> class submatrix;
    
    template<class T>
    class dmatrix : public matrix<T> {
    public:
        int ld;

        dmatrix(size_t r = 1, size_t c = 1) : matrix<T>(r*c) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            ld = r;
            // std::cout << "matrix contructor 1: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }
        
        dmatrix(size_t r, size_t c, int l) : matrix<T>(l*c) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            ld = l;
            // std::cout << "matrix contructor 2: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }

        dmatrix(const dmatrix<T>& m) : matrix<T>(m) {
            matrix<T>::nrow = m.nrow;
            matrix<T>::ncol = m.ncol;
            ld = m.ld;
            // std::cout << "matrix contructor 3: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }

        dmatrix(const csrmatrix<T>& m);
        dmatrix(const cscmatrix<T>& m);
        dmatrix(const coomatrix<T>& m);

        dmatrix(size_t r, size_t c, T* v) : matrix<T>(r*c, v) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            ld = r;
            // std::cout << "matrix contructor 4: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }

        dmatrix(size_t r, size_t c, T* v, int l) : matrix<T>(l*c, v) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            ld = l;
            // std::cout << "matrix contructor 4: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }

        dmatrix(size_t r, size_t c, const T& v) : matrix<T>(r*c, v) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            ld = r;
            // std::cout << "matrix contructor 4: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }

        dmatrix(size_t r, size_t c, T* p, int l, std::vector<T>* e) : matrix<T>(l*c, p, e) {
            matrix<T>::nrow = r;
            matrix<T>::ncol = c;
            ld = l;
            // std::cout << "matrix contructor 5: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        } 

        virtual ~dmatrix() {
            // std::cout << "matrix destructor: "
            //     << "this=" << this << std::endl;
        }
        
    public:
        
        size_t nnz() const {
            size_t z = 0;
            for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                    if (nonzero_check(array<T>::ptr[ix+jx])){
                        z++;
                    }
                }
            }
            return z;
        }

        virtual int type() const {
            return DENSE; // dense
        }

        virtual bool nonzero_check(const T& x) const;
        virtual const T& zero() const;

        // inner class
        typedef matrix_range< T, dmatrix<T> > range;
        typedef matrix_range< const T, const dmatrix<T> > const_range;
        
        // method
        range alias() {
            return range(*this, 1, matrix<T>::nrow, 1, matrix<T>::ncol);
        }
        
        range alias(sci::range x, sci::range y) {
            return range(*this, x.begin, x.end, y.begin, y.end);
        }
        
        const_range alias() const {
            return const_range(*this, 1, matrix<T>::nrow, 1, matrix<T>::ncol);
        }
        
        const_range alias(sci::range x, sci::range y) const {
            return const_range(*this, x.begin, x.end, y.begin, y.end);
        }
        
        T& operator()(int x, int y) {
            return array<T>::ptr[x-1 + (y-1)*ld];
        }
        
        const T& operator()(int x, int y) const {
            return array<T>::ptr[x-1 + (y-1)*ld];
        }
        
        dmatrix<T> operator()(sci::range x, sci::range y);
        const dmatrix<T> operator()(sci::range x, sci::range y) const;
        
        // to vec
        vector<T> operator()(int x, sci::range y) {
            size_t c = y.end - y.begin + 1;
            return vector<T>(c, &array<T>::ptr[x-1+(y.begin-1)*ld], ld, array<T>::elem);
        }
        
        const vector<T> operator()(int x, sci::range y) const {
            size_t c = y.end - y.begin + 1;
            return vector<T>(c, &array<T>::ptr[x-1+(y.begin-1)*ld], ld, array<T>::elem);
        }
        
        vector<T> operator()(sci::range x, int y) {
            size_t r = x.end - x.begin + 1;
            return vector<T>(r, &array<T>::ptr[x.begin-1+(y-1)*ld], 1, array<T>::elem);
        }
        
        const vector<T> operator()(sci::range x, int y) const {
            size_t r = x.end - x.begin + 1;
            return vector<T>(r, &array<T>::ptr[x.begin-1+(y-1)*ld], 1, array<T>::elem);
        }
        
        vector<T> diag() {
            if (matrix<T>::nrow != matrix<T>::ncol) {
                throw;
            }
            return vector<T>(matrix<T>::nrow,
                &array<T>::ptr[0], ld+1, array<T>::elem);
        }
        
        const vector<T> diag() const {
            if (matrix<T>::nrow != matrix<T>::ncol) {
                throw;
            }
            return vector<T>(matrix<T>::nrow,
                &array<T>::ptr[0], ld+1, array<T>::elem);
        }

        dmatrix<T>& operator=(const T& c) {
            for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                    array<T>::ptr[ix+jx] = c;
                }
            }
            return *this;
        }
        
        dmatrix<T>& operator=(const T* a) {
            size_t z = 0;
            for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                    array<T>::ptr[ix+jx] = a[z];
                    z++;
                }
            }
            return *this;
        }
        
        dmatrix<T>& operator=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            for (size_t j=0, jx=0, jy=0; j<matrix<T>::ncol; j++, jx+=ld, jy+=m.ld) {
                for (size_t i=0, ix=0, iy=0; i<matrix<T>::nrow; i++, ix+=1, iy+=1) {
                    array<T>::ptr[ix+jx] = m.ptr[iy+jy];
                }
            }
            return *this;
        }
        
        dmatrix<T>& operator=(const csrmatrix<T>& m);
        dmatrix<T>& operator=(const cscmatrix<T>& m);
        dmatrix<T>& operator=(const coomatrix<T>& m);

        // void shallow_copy(const matrix<T>& m) {
        //     array<T>::shallow_copy(m);
        //     nrow = m.nrow;
        //     ncol = m.ncol;
        //     ld = m.ld;
        // }

        dmatrix<T> operator+(const dmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            dmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol);
            
            for (size_t j=0, jx=0, jy=0, jz=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld, jz+=result.ld) {
                for (size_t i=0, ix=0, iy=0, iz=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1, iz+=1) {
                    result.ptr[iz+jz] = array<T>::ptr[ix+jx] + m.ptr[iy+jy];
                }
            }
            return result;
        }
        
        dmatrix<T> operator-(const dmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            dmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol);
            
            for (size_t j=0, jx=0, jy=0, jz=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld, jz+=result.ld) {
                for (size_t i=0, ix=0, iy=0, iz=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1, iz+=1) {
                    result.ptr[iz+jz] = array<T>::ptr[ix+jx] - m.ptr[iy+jy];
                }
            }
            return result;
        }
        
        dmatrix<T> operator*(const dmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            dmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol);
            
            for (size_t j=0, jx=0, jy=0, jz=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld, jz+=result.ld) {
                for (size_t i=0, ix=0, iy=0, iz=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1, iz+=1) {
                    result.ptr[iz+jz] = array<T>::ptr[ix+jx] * m.ptr[iy+jy];
                }
            }
            return result;
        }
        
        dmatrix<T> operator/(const dmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            dmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol);
            
            for (size_t j=0, jx=0, jy=0, jz=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld, jz+=result.ld) {
                for (size_t i=0, ix=0, iy=0, iz=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1, iz+=1) {
                    result.ptr[iz+jz] = array<T>::ptr[ix+jx] / m.ptr[iy+jy];
                }
            }
            return result;
        }
        
        dmatrix<T> operator%(const dmatrix<T>& m) const {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            dmatrix<T> result(matrix<T>::nrow, matrix<T>::ncol);
            
            for (size_t j=0, jx=0, jy=0, jz=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld, jz+=result.ld) {
                for (size_t i=0, ix=0, iy=0, iz=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1, iz+=1) {
                    result.ptr[iz+jz] = array<T>::ptr[ix+jx] % m.ptr[iy+jy];
                }
            }
            return result;
        }
        
        ////
        
        dmatrix<T>& operator+=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            
            for (size_t j=0, jx=0, jy=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld) {
                for (size_t i=0, ix=0, iy=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1) {
                    array<T>::ptr[ix+jx] += m.ptr[iy+jy];
                }
            }
            return *this;
        }
        
        dmatrix<T>& operator-=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            
            for (size_t j=0, jx=0, jy=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld) {
                for (size_t i=0, ix=0, iy=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1) {
                    array<T>::ptr[ix+jx] -= m.ptr[iy+jy];
                }
            }
            return *this;
        }
        
        dmatrix<T>& operator*=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            
            for (size_t j=0, jx=0, jy=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld) {
                for (size_t i=0, ix=0, iy=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1) {
                    array<T>::ptr[ix+jx] *= m.ptr[iy+jy];
                }
            }
            return *this;
        }
        
        dmatrix<T>& operator/=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            
            for (size_t j=0, jx=0, jy=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld) {
                for (size_t i=0, ix=0, iy=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1) {
                    array<T>::ptr[ix+jx] /= m.ptr[iy+jy];
                }
            }
            return *this;
        }
        
        dmatrix<T>& operator%=(const dmatrix<T>& m) {
            if (matrix<T>::nrow != m.nrow || matrix<T>::ncol != m.ncol) {
                throw;
            }
            
            for (size_t j=0, jx=0, jy=0;
                 j<matrix<T>::ncol;
                 j++, jx+=ld, jy+=m.ld) {
                for (size_t i=0, ix=0, iy=0;
                     i<matrix<T>::nrow;
                     i++, ix+=1, iy+=1) {
                    array<T>::ptr[ix+jx] %= m.ptr[iy+jy];
                }
            }
            return *this;
        }
        
        dmatrix<T>& operator+=(const T& v) {
            for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                    array<T>::ptr[ix+jx] += v;
                }
            }
            return *this;
        }

        dmatrix<T>& operator-=(const T& v) {
            for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                    array<T>::ptr[ix+jx] -= v;
                }
            }
            return *this;
        }

        dmatrix<T>& operator*=(const T& v) {
            for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                    array<T>::ptr[ix+jx] *= v;
                }
            }
            return *this;
        }

        dmatrix<T>& operator/=(const T& v) {
            for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                    array<T>::ptr[ix+jx] /= v;
                }
            }
            return *this;
        }

        dmatrix<T>& operator%=(const T& v) {
            for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                    array<T>::ptr[ix+jx] %= v;
                }
            }
            return *this;
        }

        virtual std::ostream& print(std::ostream& os) const {
            for (size_t i=0, ix=0; i<matrix<T>::nrow; i++, ix+=1) {
                for (size_t j=0, jx=0; j<matrix<T>::ncol; j++, jx+=ld) {
                    os << array<T>::ptr[ix+jx] << " ";
                }
                os << std::endl;
            }
            return os;
        }
        
        template <class sT>
        friend std::ostream& operator<< (std::ostream& os, const dmatrix<sT>& v);

    };
    
    template<class T>
    dmatrix<T> dmatrix<T>::operator()(sci::range x, sci::range y) {
        size_t r = x.end - x.begin + 1;
        size_t c = y.end - y.begin + 1;
        return dmatrix<T>(r, c, &array<T>::ptr[x.begin-1+(y.begin-1)*ld], ld, array<T>::elem);
    }
    
    template<class T>
    const dmatrix<T> dmatrix<T>::operator()(sci::range x, sci::range y) const {
        size_t r = x.end - x.begin + 1;
        size_t c = y.end - y.begin + 1;
        return dmatrix<T>(r, c, &array<T>::ptr[x.begin-1+(y.begin-1)*ld], ld, array<T>::elem);
    }

    template<class T>
    std::ostream& operator<<(std::ostream& os, const dmatrix<T>& v) {
        return v.print(os);
    }
    
}
