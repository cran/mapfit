/*
 Scientific sparse matrix package
 */

#pragma once

#include "sci_matrix.h"

 namespace sci {

    /* Sparse Matrix Base */
    template <class T>
    class spmatrix : public matrix<T> {
    public:
        size_t nnz;

        spmatrix(size_t s) : matrix<T>(s) {}
        spmatrix(const array<T>& a) : matrix<T>(a) {}
        spmatrix(size_t s, T* v) : matrix<T>(s, v) {}
        spmatrix(size_t s, const T& v) : matrix<T>(s, v) {}
        spmatrix(size_t s, T* p, std::vector<T>* e) : matrix<T>(s, p, e) {}

        virtual ~spmatrix() {
        }
        
        typename array<T>::range get_value() {
            return array<T>::alias(sci::range(1,static_cast<int>(nnz)));
        }
        
        typename array<T>::const_range get_value() const {
            return array<T>::alias(sci::range(1,static_cast<int>(nnz)));
        }
        
        // print
        
        virtual std::ostream& print(std::ostream& os) const = 0;
        
        template <class sT>
        friend std::ostream& operator<<(std::ostream& os, const spmatrix<sT>& m);

    };
    
    template<class T>
    std::ostream& operator<<(std::ostream& os, const spmatrix<T>& m) {
        return m.print(os);
    }
}

