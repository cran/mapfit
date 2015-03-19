/*
 Scientific vector and matrix package
 */

#pragma once

#include "sci_array.hh"

namespace sci {
    
    // template<class T> class subvector;

    ///////////////////////// vector
    
    template<class T>
    class vector : public array<T> {
    public:
        int inc;

        vector(size_t s = 1) : array<T>(s) {
            inc = 1;
            // std::cout << "vector contructor 1: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }
        
        vector(size_t s, int i) : array<T>(s*i) {
            inc = i;
            // std::cout << "vector contructor 2: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }

        vector(const vector<T>& v) : array<T>(v) {
            inc = v.inc;
            // std::cout << "vector contructor 3: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }
        
        vector(size_t s, T* v, int i = 1) : array<T>(s, v) {
            inc = i;
            // std::cout << "vector contructor 4: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }

        vector(size_t s, const T& v) : array<T>(s, v) {
            inc = 1;
            // std::cout << "vector contructor 4: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        }

        vector(size_t s, T* p, int i, std::vector<T>* e) : array<T>(s, p, e) {
            inc = i;
            // std::cout << "vector contructor 6: "
            //     << " this=" << this
            //     << " elem=" << array<T>::elem.get()
            //     << " ptr=" << array<T>::ptr << std::endl;
        } 
        
        virtual ~vector() {
            // std::cout << "vector destructor: "
            //     << "this=" << this << std::endl;
        }

    public:
        
        // inner class
        typedef vector_range< T, vector<T> > range;
        typedef vector_range< const T, const vector<T> > const_range;
        
        // method
        range alias() {
            return range(*this, 1, array<T>::size);
        }
        
        range alias(sci::range x) {
            return range(*this, x.begin, x.end);
        }
        
        const_range alias() const {
            return const_range(*this, 1, array<T>::size);
        }
        
        const_range alias(sci::range x) const {
            return const_range(*this, x.begin, x.end);
        }
        
        // elements: 1 origin
        T& operator()(int x) {
            return array<T>::ptr[(x-1)*inc];
        }
        
        const T& operator()(int x) const {
            return array<T>::ptr[(x-1)*inc];
        }
        
        vector<T> operator()(sci::range x);
        const vector<T> operator()(sci::range x) const;
        
        // equal
        vector<T>& operator=(const T& c) {
            for (size_t i=0, ix=0; i<array<T>::size; i++, ix+=inc) {
                array<T>::ptr[ix] = c;
            }
            return *this;
        }
        
        vector<T>& operator=(const T* a) {
            for (size_t i=0, ix=0; i<array<T>::size; i++, ix+=inc) {
                array<T>::ptr[ix] = a[i];
            }
            return *this;
        }
        
        vector<T>& operator=(const vector<T>& v) {
            if (array<T>::size != v.size) {
                throw;
            }
            for (size_t i=0, ix=0, iy=0; i<array<T>::size;
                 i++, ix+=inc, iy+=v.inc) {
                array<T>::ptr[ix] = v.ptr[iy];
            }
            return *this;
        }

        // void shallow_copy(const vector<T>& v) {
        //     array<T>::shallow_copy(v);
        //     inc = v.inc;
        // }

        // operators
        vector<T> operator+(const vector<T>& v) const {
            if (array<T>::size != v.size) {
                throw;
            }
            vector<T> result(array<T>::size);
            
            for (size_t i=0, ix=0, iy=0, iz=0;
                 i<array<T>::size;
                 i++, ix+=inc, iy+=v.inc, iz+=result.inc) {
                result.ptr[iz] = array<T>::ptr[ix] + v.ptr[iy];
            }
            return result;
        }
        
        vector<T> operator-(const vector<T>& v) const {
            if (array<T>::size != v.size) {
                throw;
            }
            vector<T> result(array<T>::size);
            
            for (size_t i=0, ix=0, iy=0, iz=0;
                 i<array<T>::size;
                 i++, ix+=inc, iy+=v.inc, iz+=result.inc) {
                result.ptr[iz] = array<T>::ptr[ix] - v.ptr[iy];
            }
            return result;
        }
        
        vector<T> operator*(const vector<T>& v) const {
            if (array<T>::size != v.size) {
                throw;
            }
            vector<T> result(array<T>::size);
            
            for (size_t i=0, ix=0, iy=0, iz=0;
                 i<array<T>::size;
                 i++, ix+=inc, iy+=v.inc, iz+=result.inc) {
                result.ptr[iz] = array<T>::ptr[ix] * v.ptr[iy];
            }
            return result;
        }
        
        vector<T> operator/(const vector<T>& v) const {
            if (array<T>::size != v.size) {
                throw;
            }
            vector<T> result(array<T>::size);
            
            for (size_t i=0, ix=0, iy=0, iz=0;
                 i<array<T>::size;
                 i++, ix+=inc, iy+=v.inc, iz+=result.inc) {
                result.ptr[iz] = array<T>::ptr[ix] / v.ptr[iy];
            }
            return result;
        }
        
        vector<T> operator%(const vector<T>& v) const {
            if (array<T>::size != v.size) {
                throw;
            }
            vector<T> result(array<T>::size);
            
            for (size_t i=0, ix=0, iy=0, iz=0;
                 i<array<T>::size;
                 i++, ix+=inc, iy+=v.inc, iz+=result.inc) {
                result.ptr[iz] = array<T>::ptr[ix] % v.ptr[iy];
            }
            return result;
        }
        
        ////
        
        vector<T>& operator+=(const vector<T>& v) {
            if (array<T>::size != v.size) {
                throw;
            }
            for (size_t i=0, ix=0, iy=0;
                 i < array<T>::size; i++, ix+=inc, iy+=v.inc) {
                array<T>::ptr[ix] += v.ptr[iy];
            }
            return *this;
        }
        
        vector<T>& operator-=(const vector<T>& v) {
            if (array<T>::size != v.size) {
                throw;
            }
            for (size_t i=0, ix=0, iy=0;
                 i<array<T>::size; i++, ix+=inc, iy+=v.inc) {
                array<T>::ptr[ix] -= v.ptr[iy];
            }
            return *this;
        }
        
        vector<T>& operator*=(const vector<T>& v) {
            if (array<T>::size != v.size) {
                throw;
            }
            for (size_t i=0, ix=0, iy=0;
                 i<array<T>::size; i++, ix+=inc, iy+=v.inc) {
                array<T>::ptr[ix] *= v.ptr[iy];
            }
            return *this;
        }
        
        vector<T>& operator/=(const vector<T>& v) {
            if (array<T>::size != v.size) {
                throw;
            }
            for (size_t i=0, ix=0, iy=0;
                 i<array<T>::size; i++, ix+=inc, iy+=v.inc) {
                array<T>::ptr[ix] /= v.ptr[iy];
            }
            return *this;
        }
        
        vector<T>& operator%=(const vector<T>& v) {
            if (array<T>::size != v.size) {
                throw;
            }
            for (size_t i=0, ix=0, iy=0;
                 i<array<T>::size; i++, ix+=inc, iy+=v.inc) {
                array<T>::ptr[ix] %= v.ptr[iy];
            }
            return *this;
        }
        
        vector<T>& operator+=(const T& v) {
            for (size_t i=0, ix=0; i < array<T>::size; i++, ix+=inc) {
                array<T>::ptr[ix] += v;
            }
            return *this;
        }

        vector<T>& operator-=(const T& v) {
            for (size_t i=0, ix=0; i < array<T>::size; i++, ix+=inc) {
                array<T>::ptr[ix] -= v;
            }
            return *this;
        }

        vector<T>& operator*=(const T& v) {
            for (size_t i=0, ix=0; i < array<T>::size; i++, ix+=inc) {
                array<T>::ptr[ix] *= v;
            }
            return *this;
        }

        vector<T>& operator/=(const T& v) {
            for (size_t i=0, ix=0; i < array<T>::size; i++, ix+=inc) {
                array<T>::ptr[ix] /= v;
            }
            return *this;
        }

        vector<T>& operator%=(const T& v) {
            for (size_t i=0, ix=0; i < array<T>::size; i++, ix+=inc) {
                array<T>::ptr[ix] %= v;
            }
            return *this;
        }

        ////// print
        virtual std::ostream& print(std::ostream& os) const {
            for (size_t i=0, ix=0; i<array<T>::size; i++, ix+=inc) {
                os << array<T>::ptr[ix] << " ";
            }
            return os;
        }
        
        template <class sT>
        friend std::ostream& operator<< (std::ostream& os,
                                         const vector<sT>& v);
        
    };
    
    template<class T>
    vector<T> vector<T>::operator()(sci::range x) {
        size_t s = x.end - x.begin + 1;
        return vector<T>(s, &array<T>::ptr[(x.begin-1)*inc], inc, array<T>::elem);
    }
    
    template<class T>
    const vector<T> vector<T>::operator()(sci::range x) const {
        size_t s = x.end - x.begin + 1;
        return vector<T>(s, &array<T>::ptr[(x.begin-1)*inc], inc, array<T>::elem);
    }

    template<class T>
    std::ostream& operator<<(std::ostream& os, const vector<T>& v) {
        return v.print(os);
    }

    // generic method
    template<class T>
    int imax(const vector<T>& x) {
        int m = 1;
        for (size_t i=2; i<=x.size; i++) {
            if (x(m) < x(i)) {
                m = i;
            }
        }
        return m;
    }

    template<class T>
    int imin(const vector<T>& x) {
        int m = 1;
        for (size_t i=2; i<=x.size; i++) {
            if (x(m) > x(i)) {
                m = i;
            }
        }
        return m;
    }
}
