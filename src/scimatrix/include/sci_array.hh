/*
 Scientific vector and matrix package
 */

#pragma once

#include <vector>
#include <iostream>
#include <memory>

#include "sci_range.hh"

namespace sci {

    ///////////////////////// array
  
    template<class T>
    class array : public object {
    public:
        bool master;
        T* ptr;
        size_t size;
//        std::shared_ptr< std::vector<T> > elem;
        std::vector<T>* elem;

    public:
        array(size_t s = 1) : elem(new std::vector<T>(s)) {
            master = true;
            size = s;
            ptr = &((*elem)[0]);
            // std::cout << "array contructor 1: "
            //     << "this=" << this
            //     << "elem=" << elem.get()
            //     << " ptr=" << ptr << std::endl;
        }
        
        array(const array<T>& a) : elem(new std::vector<T>(&a.ptr[0], &a.ptr[a.size])) {
            master = true;
            size = a.size;
            ptr = &((*elem)[0]);

            // std::cout << "array contructor 2: "
            //     << "this=" << this
            //     << "elem=" << elem.get()
            //     << " ptr=" << ptr << std::endl;
        }
        
        array(size_t s, T* v) : elem(0) {
            master = false;
            size = s;
            ptr = v;
            // std::cout << "array contructor 3: "
            //     << "this=" << this
            //     << "elem=" << elem.get()
            //     << " ptr=" << ptr << std::endl;
        }

        array(size_t s, const T& t) : elem(new std::vector<T>(s,t)) {
            master = true;
            size = s;
            ptr = &((*elem)[0]);
            // std::cout << "array contructor 3: "
            //     << "this=" << this
            //     << "elem=" << elem.get()
            //     << " ptr=" << ptr << std::endl;
        }

        virtual ~array() {
            // std::cout << "array destructor: "
            //     << "this=" << this << std::endl;
            if (master) {
                delete elem;
            }
        }
        
    protected:
        array(size_t s, T* p, std::vector<T>* e) : elem(e) {
            master = false;
            size = s;
            ptr = p;
            // std::cout << "contructor 4: "
            //     << "this=" << this
            //     << "elem=" << elem.get()
            //     << " ptr=" << ptr << std::endl;
        }

    public:
        // inner
        typedef array_range< T, array<T> > range;
        typedef array_range< const T, const array<T> > const_range;
        
        // method
        range alias() {
            return range(*this, 0, size-1);
        }
        
        range alias(sci::range r) {
            return range(*this, r.begin, r.end);
        }
        
        const_range alias() const {
            return const_range(*this, 0, size-1);
        }
        
        const_range alias(sci::range r) const {
            return const_range(*this, r.begin, r.end);
        }
        
        T& operator[](int x) {
            return ptr[x];
        }
        
        const T& operator[](int x) const {
            return ptr[x];
        }
        
        array<T> operator[](sci::range r);
        const array<T> operator[](sci::range r) const;
                
        // equal
        array<T>& operator=(const T& v) {
            for (size_t i=0; i<size; i++) {
                ptr[i] = v;
            }
            return *this;
        }
        
        array<T>& operator=(const T* a) {
            for (size_t i=0; i<size; i++) {
                ptr[i] = a[i];
            }
            return *this;
        }
        
        array<T>& operator=(const array<T>& a) {
            if (size != a.size) {
                throw;
            }
            for (size_t i=0; i<size; i++) {
                ptr[i] = a.ptr[i];
            }
            return *this;
        }

        // template <class S>
        // array_base<T>& operator=(const array_base<S>& a) {
        //     if (size != a.size) {
        //         throw;
        //     }
        //     for (int i=0; i<size; i++) {
        //         ptr[i] = a.ptr[i];
        //     }
        //     return *this;
        // }

        // void shallow_copy(const array_base<T>& v) {
        //     size = v.size;
        //     ptr = v.ptr;
        // }

        // operators
        array<T> operator+(const array<T>& a) const {
            if (size != a.size) {
                throw;
            }
            array<T> result(size);
            for (size_t i=0; i<size; i++) {
                result.ptr[i] = ptr[i] + a.ptr[i];
            }
            return result;
        }
        
        array<T> operator-(const array<T>& a) const {
            if (size != a.size) {
                throw;
            }
            array<T> result(size);
            for (size_t i=0; i<size; i++) {
                result.ptr[i] = ptr[i] - a.ptr[i];
            }
            return result;
        }
        
        array<T> operator*(const array<T>& a) const {
            if (size != a.size) {
                throw;
            }
            array<T> result(size);
            for (size_t i=0; i<size; i++) {
                result.ptr[i] = ptr[i] * a.ptr[i];
            }
            return result;
        }
        
        array<T> operator/(const array<T>& a) const {
            if (size != a.size) {
                throw;
            }
            array<T> result(size);
            for (size_t i=0; i<size; i++) {
                result.ptr[i] = ptr[i] / a.ptr[i];
            }
            return result;
        }
        
        array<T> operator%(const array<T>& a) const {
            if (size != a.size) {
                throw;
            }
            array<T> result(size);
            for (size_t i=0; i<size; i++) {
                result.ptr[i] = ptr[i] % a.ptr[i];
            }
            return result;
        }
        
        //
        array<T>& operator+=(const array<T>& a) {
            if (size != a.size) {
                throw;
            }
            for (size_t i=0; i<size; i++) {
                ptr[i] += a.ptr[i];
            }
            return *this;
        }
        
        array<T>& operator-=(const array<T>& a) {
            if (size != a.size) {
                throw;
            }
            for (size_t i=0; i<size; i++) {
                ptr[i] -= a.ptr[i];
            }
            return *this;
        }
        
        array<T>& operator*=(const array<T>& a) {
            if (size != a.size) {
                throw;
            }
            for (size_t i=0; i<size; i++) {
                ptr[i] *= a.ptr[i];
            }
            return *this;
        }
        
        array<T>& operator/=(const array<T>& a) {
            if (size != a.size) {
                throw;
            }
            for (size_t i=0; i<size; i++) {
                ptr[i] /= a.ptr[i];
            }
            return *this;
        }
        
        array<T>& operator%=(const array<T>& a) {
            if (size != a.size) {
                throw;
            }
            for (size_t i=0; i<size; i++) {
                ptr[i] %= a.ptr[i];
            }
            return *this;
        }
        
        array<T>& operator+=(const T& v) {
            for (size_t i=0; i<size; i++) {
                ptr[i] += v;
            }
            return *this;
        }
        
        array<T>& operator-=(const T& v) {
            for (size_t i=0; i<size; i++) {
                ptr[i] -= v;
            }
            return *this;
        }

        array<T>& operator*=(const T& v) {
            for (size_t i=0; i<size; i++) {
                ptr[i] *= v;
            }
            return *this;
        }

        array<T>& operator/=(const T& v) {
            for (size_t i=0; i<size; i++) {
                ptr[i] /= v;
            }
            return *this;
        }

        array<T>& operator%=(const T& v) {
            for (size_t i=0; i<size; i++) {
                ptr[i] %= v;
            }
            return *this;
        }

        // print
        virtual std::ostream& print(std::ostream& os) const {
            for (size_t i=0; i<size; i++) {
                os << ptr[i] << " ";
            }
            return os;
        }
        
        template <class sT>
        friend std::ostream& operator<< (std::ostream& os, const array<sT>& v);

        void resize(size_t s) {
            size = s;
            (*elem).resize(s);
            ptr = &((*elem)[0]);
        }
    };

    template<class T>
    array<T> array<T>::operator[](sci::range r) {
        size_t s = r.end - r.begin + 1;
        return array<T>(s, &ptr[r.begin], elem);
    }
    
    template<class T>
    const array<T> array<T>::operator[](sci::range r) const {
        size_t s = r.end - r.begin + 1;
        return array<T>(s, &ptr[r.begin], elem);
    }

    template<class T>
    std::ostream& operator<<(std::ostream& os, const array<T>& v) {
        return v.print(os);
    }
    
}
