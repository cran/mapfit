////////
///  range
///////

#pragma once

#include <iostream>

namespace sci {
    
    /////////////////////////////////////////
    // prototype
    /////////////////////////////////////////
    
    // template <class T> class array;
    // template <class T> class vector;
    // template <class T> class matrix;
    
    //  typedef size_t int;
    class object {
    public:
        virtual ~object() { }
    };
    
    //////////////////////////////////////
    
    class range : public object {
    public:
        range(int b, int e) {
            begin = b;
            end = e;
        }
        
        ~range() { }
        
        int begin;
        int end;

        int size() {
            return end - begin + 1;
        }
    };
    
    ///////////////
    
    template <class T, class Ary>
    class array_range : public object {
    public:
        // constractor
        array_range(Ary& a, int lb, int ub) : ary(a) {
            //      std::cout << "constructor: range (array) " << this << std::endl;
            albound = lb;
            aubound = ub;
        }
        
        virtual ~array_range() {
            //      std::cout << "destructor: range (array) " << this << std::endl;
        }
        
        // method
        T& operator[](int i) {
            if (i < albound || aubound < i) {
                throw;
            }
            return ary.ptr[i-albound];
        }
        
        Ary operator[](sci::range r) {
            if (r.begin < albound || aubound < r.end) {
                throw;
            }
            return Ary(r.end-r.begin+1, &ary.ptr[r.begin-albound]);
        }

    private:
        int albound;
        int aubound;
        Ary& ary;
    };
    
    /////////////////////////
    
    template <class T, class Vec>
    class vector_range : public object {
    public:
        // constractor
        vector_range(Vec& v, int xlb, int xub) : vec(v) {
            //      std::cout << "constructor: range (vec) " << this << std::endl;
            xlbound = xlb;
            xubound = xub;
        }
        
        virtual ~vector_range() {
            //      std::cout << "destructor: range (vec) " << this << std::endl;
        }
        
        // method
        T& operator()(int x) {
            if (x < xlbound || xubound < x) {
                throw;
            }
            return vec.ptr[(x-xlbound)*vec.inc];
        }
        
        Vec operator()(sci::range x) {
            if (x.begin < xlbound || xubound < x.end) {
                throw;
            }
            return Vec(x.end-x.begin+1,
                       &vec.ptr[(x.begin-xlbound)*vec.inc],
                       vec.inc);
        }
        
    private:
        int xlbound;
        int xubound;
        Vec& vec;
    };
    
    /////////////////////////
    
    template <class T, class Mat>
    class matrix_range : public object {
    public:
        // constractor
        matrix_range(Mat& m, int xlb, int xub,
                     int ylb, int yub)
        : mat(m) {
            //      std::cout << "constructor: range (mat) " << this << std::endl;
            xlbound = xlb;
            xubound = xub;
            ylbound = ylb;
            yubound = yub;
        }
        
        virtual ~matrix_range() {
            //      std::cout << "destructor: range (mat) " << this << std::endl;
        }
        
        // method
        T& operator()(int x, int y) {
            if (x < xlbound || xubound < x || y < ylbound || yubound < y) {
                throw;
            }
            return mat.ptr[x-xlbound + (y-ylbound)*mat.ld];
        }
        
        Mat operator()(sci::range x, sci::range y) {
            if (x.begin < xlbound || xubound < x.end
                || y.begin < ylbound || ylbound < y.end) {
                throw;
            }
            return Mat(x.end-x.begin+1, y.end-y.begin+1,
                       &mat.ptr[x.begin-xlbound
                                + (y.begin-ylbound)*mat.ld],
                       mat.ld);
        }
        
    private:
        int xlbound;
        int xubound;
        int ylbound;
        int yubound;
        Mat& mat;
    };
    
}

