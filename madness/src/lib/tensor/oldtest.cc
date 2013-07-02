/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: test.cc 1105 2009-03-28 20:58:10Z rjharrison $
*/

//  !!! This file is deprecated ... contents are being slowly
//  migrated into test.cc (using Google test environment)


#include <tensor/tensor.h>

#include <iostream>
#include <cstdio>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <ctime>


/// \file tensor/oldtest.cc
/// \brief Test code for Tensor, TensorIterator, SliceTensor, etc.

using namespace madness;

typedef std::complex<float> float_complex;
typedef std::complex<double> double_complex;

void error(const char *msg, int code) {
    std::cerr << msg << " " << code << std::endl;
    std::exit(1);
}

template <typename T> double mynorm(T x) {
    return ((double) x)*x;
}

template <> double mynorm<float_complex>(float_complex x) {
    return (double) std::norm(x);
}

template <> double mynorm<double_complex>(double_complex x) {
    return (double) std::norm(x);
}

template <typename T, typename Q>
inline
bool
check(const T& t, const Q& q, double tol=1e-7) {
    double err = std::abs(t - (T)q)/std::max<double>(std::abs(t),1.0);
    bool ok = (err <= tol);
    if (!ok) {
        std::cout.setf(std::ios::scientific);
        std::cout << "check failed " << t << " " << q << " " << err << " " << tol << std::endl;
    }
    return ok;
}


template <typename T, typename Q> void Test1() {
    //
    // Test basic declarations & type conversions.
    //

    try {

        Tensor<T> t;
        std::cout << "t\n";
        std::cout << t.ndim() << std::endl;
        std::cout << t.size() << std::endl;

        Tensor<Q> q1(1);
        std::cout << "q1\n";
        std::cout << q1.ndim() << std::endl;
        std::cout << q1.size() << std::endl;
        std::cout << q1.dim(0) << std::endl;
        std::cout << q1.stride(0) << std::endl;
        Tensor<Q> q2(1,2);
        std::cout << q1.ndim() << std::endl;
        Tensor<Q> q3(1,2,3);
        std::cout << q1.ndim() << std::endl;
        Tensor<Q> q4(1,2,3,4);
        std::cout << q1.ndim() << std::endl;
        Tensor<Q> q5(1,2,3,4,5);
        std::cout << q1.ndim() << std::endl;
        Tensor<Q> q6(1,2,3,4,5,6);
        std::cout << q1.ndim() << std::endl;

        long d[] = {4,3,2,1};
        Tensor<Q> qg(4L,d);

        // should be initialized to zero
        t = Tensor<T>(q1);
        std::cout << t;
        ITERATOR1(t,if (q1(IND1) != 0) error("test1: failed",10));
        t = Tensor<T>(q2);
        ITERATOR2(t,if (q2(IND2) != 0) error("test1: failed",20));
        t = Tensor<T>(q3);
        ITERATOR3(t,if (q3(IND3) != 0) error("test1: failed",30));
        t = Tensor<T>(q4);
        ITERATOR4(t,if (q4(IND4) != 0) error("test1: failed",40));
        t = Tensor<T>(q5);
        ITERATOR5(t,if (q5(IND5) != 0) error("test1: failed",50));
        t = Tensor<T>(q6);
        ITERATOR6(t,if (q6(IND6) != 0) error("test1: failed",60));

        long count;
        q1.fillindex();
        count = 0;
        ITERATOR1(q1,if (q1(IND1) != count++) error("test1: failed",110));
        q2.fillindex();
        count = 0;
        ITERATOR2(q2,if (q2(IND2) != count++) error("test1: failed",210));
        q3.fillindex();
        count = 0;
        ITERATOR3(q3,if (q3(IND3) != count++) error("test1: failed",310));
        q4.fillindex();
        count = 0;
        ITERATOR4(q4,if (q4(IND4) != count++) error("test1: failed",410));
        q5.fillindex();
        count = 0;
        ITERATOR5(q5,if (q5(IND5) != count++) error("test1: failed",510));
        q6.fillindex();
        count = 0;
        ITERATOR6(q6,if (q6(IND6) != count++) error("test1: failed",610));

        q1.fillrandom();
        q2.fillrandom();
        q3.fillrandom();
        q4.fillrandom();
        q5.fillrandom();
        q6.fillrandom();
        if (q1.id() == TensorTypeData<long>::id) {
            q1.fillindex();
            q2.fillindex();
            q3.fillindex();
            q4.fillindex();
            q5.fillindex();
            q6.fillindex();
        }


        // if T != Q, this is type conversion with deep copy performed
        // by first converting Q to a new tensor of type T, and then
        // doing shallow assignment of same type.

        // if T == Q, this is actually a shallow copy ...
        t = Tensor<T>(q1);
        ITERATOR1(t,if (!check(t(IND1),q1(IND1))) error("test1: failed",1));
        t = Tensor<T>(q2);
        std::cout << "q2\n" << q2 << std::endl;
        std::cout << "t\n" << t << std::endl;
        ITERATOR2(t,if (!check(t(IND2),q2(IND2))) error("test1: failed",2));
        t = Tensor<T>(q3);
        ITERATOR3(t,if (!check(t(IND3),q3(IND3))) error("test1: failed",3));
        t = Tensor<T>(q4);
        ITERATOR4(t,if (!check(t(IND4),q4(IND4))) error("test1: failed",4));
        t = Tensor<T>(q5);
        ITERATOR5(t,if (!check(t(IND5),q5(IND5))) error("test1: failed",5));
        t = Tensor<T>(q6);
        ITERATOR6(t,if (!check(t(IND6),q6(IND6))) error("test1: failed",6));

    }
    catch (TensorException e) {
        std::cout << "This exception is unexpected\n";
        std::cout << e;
        std::exit(1);
    }

    std::cout << "Test1<" << tensor_type_names[TensorTypeData<T>::id] << "," <<
              tensor_type_names[TensorTypeData<Q>::id] << "> OK" << std::endl;
}


template <typename T> void Test3() {
    // printing

    Tensor<T> t(3,4,5);
    ITERATOR3(t,t(_i,_j,_k) = _i*10000 + _j*100 + _k);

    std::cout << std::endl << t ;

    std::cout << "Test3<" << tensor_type_names[TensorTypeData<T>::id] <<
              "> ... verify elements are 10000*i + 100*j + k \n";
}

template <class T> void Test5() {
    // arithmetic
    Tensor<T> a(2,7,5,2,11), b(2,7,5,2,11), c(2,7,5,2,11), d(2,7,5,2,11), e, f;
    a.fillrandom();
    b.fillrandom();
    c.fillrandom();
    d.fillrandom();

    if (a.id() == TensorTypeData<long>::id) {
        // Random number generates really big numbers that screw up the
        // testing logic and cause overflow ... replace with index
        a.fillindex();
        b.fillindex();
        c.fillindex();
        d.fillindex();
    }

    e = a*0.125;
    ITERATOR5(e,if (e(IND5) != (T)(a(IND5)*0.125)) error("test5: failed",0));

    f = 0.125*a;
    ITERATOR5(e,if (f(IND5) != (T)(0.125*a(IND5))) error("test5: failed",1));
    ITERATOR5(e,if (e(IND5) != f(IND5)) {
        std::cout << e(IND5) <<" " << f(IND5) << std::endl;
        error("test5: failed",2);
    }
             );

    e = a + b;
    std::cout.setf(std::ios::scientific);
    ITERATOR5(e,if (!check(e(IND5),(T)(a(IND5)+b(IND5)))) {
        std::cout << e(IND5) << " " << a(IND5)+b(IND5) << " " << e(IND5)-(a(IND5)+b(IND5)) << std::endl;
        error("test5: failed",200);
    }
             );

    e = a - b;
    ITERATOR5(e,if (!check(e(IND5),a(IND5)-b(IND5))) error("test5: failed",3));

    e = b/32.0;
    ITERATOR5(e,if (!check(e(IND5), b(IND5)/32.0)) error("test5: failed",4));

    e = a*3.14159 + b;
    ITERATOR5(e,if (std::abs(e(IND5) - (T)(a(IND5)*3.14159 + b(IND5))) > 1e-7) { // 1e-7 for float
        std::cout << _i << " " << _j << " " << _k << " " << _l << " " << _m << " " << e(IND5) << " " << (a(IND5)*3.14159 + b(IND5));
        error("test5: failed",5);
    }
             );

    if (a.id() == TensorTypeData<long>::id) {
        // Tensor expression does type conversion to integer for each
        // intermediate ... elementwise converts only at the end
        std::printf("skipping test6 for long\n");
    }
    else {
        e = (a*3.14159 + b/37.0 - 81.0*c + d*99.0)/2.0;
        ITERATOR5(e,if (std::abs(e(IND5) - (T)((a(IND5)*3.14159 + b(IND5)/37.0 - 81.0*c(IND5) + d(IND5)*99.0)/2.0)) > 1e-12)
                  error("test5: failed",6));
    }

    // unary-, inplace operations

    f = -a;
    ITERATOR5(e,if (f(IND5) != -a(IND5)) error("test5: failed",7));

    f = copy(e);
    f += a;
    ITERATOR5(e,if (std::abs(f(IND5) - e(IND5) - a(IND5)) > 1e-12)
              error("test5: failed",8));

    f = copy(e);
    f -= a;
    ITERATOR5(e,if (std::abs(f(IND5) - e(IND5) + a(IND5)) > 1e-12)
              error("test5: failed",9));

    f = copy(e);
    f *= (T) 2;
    ITERATOR5(e,if (std::abs(f(IND5) - ((T) 2)*e(IND5)) > 1e-12)
              error("test5: failed",10));

    f = copy(e);
    f *= (T) 2;
    ITERATOR5(e,if (std::abs(f(IND5) - (T) 2*e(IND5)) > 1e-12)
              error("test5: failed",11));

    f = copy(e);
    f.emul(a);
    ITERATOR5(e,if (std::abs(f(IND5) - e(IND5)*a(IND5)) > 1e-12)
              error("test5: failed",12));

    f = copy(e);
    f.gaxpy((T) 3.14159, a, (T) 2.71828);
    ITERATOR5(e,if (!check(f(IND5),(e(IND5)*((T) 3.14159) + a(IND5)*((T) 2.71828)))) {
        std::cout << "F\n" << f;
        std::cout << "X\n" << (e*((T)3.14159)+ a*((T) 2.71828)) << std::endl;
        error("test5: failed",13);
    }
             );

    std::cout << "Test5<" << tensor_type_names[TensorTypeData<T>::id] << "> OK\n";
}


template <class T> void Test6() {
    // sum, min, max, absmin, absmax, normf, trace, and all that
    // sum already tested above

    Tensor<T> f(9,5,3,8), g(5,9,8,3);
    long ind[TENSOR_MAXDIM];

    f.fill(0);
    f(3,2,1,5) = (T) 1;
    if (f.id() != TensorTypeData<float_complex>::id &&
            f.id() != TensorTypeData<double_complex>::id) {
        if (f.max() != (T) 1) error("test6: failed",0);
        if (f.max(ind) != (T) 1) error("test6: failed",1);
        if (ind[0]!=3 || ind[1]!=2 || ind[2]!=1 || ind[3]!=5)
            error("test6: failed",2);
    }

    f.fill(0);
    f(3,2,1,5) = (T) -1;
    if (long(f.absmax()+0.0001) != 1) error("test6: failed",3);
    if (long(f.absmax(ind)+0.0001) != 1) error("test6: failed",4);
    if (ind[0]!=3 || ind[1]!=2 || ind[2]!=1 || ind[3]!=5)
        error("test6: failed",5);

    f.fill(1);
    f(3,2,1,5) = (T) -2;
    if (f.id() != TensorTypeData<float_complex>::id &&
            f.id() != TensorTypeData<double_complex>::id) {
        if (f.min() != (T) -2) error("test6: failed",6);
        if (f.min(ind) != (T) -2) error("test6: failed",7);
        if (ind[0]!=3 || ind[1]!=2 || ind[2]!=1 || ind[3]!=5)
            error("test6: failed",8);
    }

    f.fill(-12);
    f(3,2,1,5) = (T) -2;
    if (long(f.absmin()+1e-5) != 2) error("test6: failed",9);
    if (long(f.absmin(ind)+1e-5) != 2) error("test6: failed",10);
    if (ind[0]!=3 || ind[1]!=2 || ind[2]!=1 || ind[3]!=5)
        error("test6: failed",11);

    f.fillindex();
    double x = f.normf(), s=0.0;
    ITERATOR4(f,s+=mynorm(f(IND4)));
    std::cout << x << " " << sqrt(s) << "\n";
    if (std::abs(std::sqrt(s) - x) > x*1e-6) error("test6: failed",12);

    g = g.fillindex().swapdim(0,1).swapdim(-1,-2);
    T y = f.trace(g);
    T z = g.trace(f);
    std::cout << y << " " << z << std::endl;
    if (std::abs(y-z) > 2e-6*std::abs(y)) error("test6: failed",13);
    T q = 0;
    ITERATOR4(f,q+=f(IND4)*g(IND4));
    if (std::abs(y-q) > 1e-6*std::abs(y)) error("test6: failed",14);


//   double_complex sum1=0, sum2=t4.trace(t5), sum3=t5.trace(t4);

//   ITERATOR(t4,sum1+=t4(_i,_j,_k,_l,_m,_n)*t5(_i,_j,_k,_l,_m,_n));
//   std::printf("trace test: %f %f     %f %f    %f %f\n",
// 	      std::real(sum1),std::imag(sum1),
// 	      std::real(sum2),std::imag(sum2),
// 	      std::real(sum3),std::imag(sum3));


    std::cout << "Test6<" << tensor_type_names[TensorTypeData<T>::id] << "> OK\n";

}

template <class T> void Test7() {
    // inner, outer, transform ...
    Tensor<T> fred(12,13),dave(14,13);

    fred.fillrandom();
    dave.fillrandom();

    Tensor<T> george = outer(fred,dave);

    ITERATOR4(george,
              if (std::abs(george(_i,_j,_k,_l) - fred(_i,_j)*dave(_k,_l)) > 1e-7)
              error("test7: failed",1));

    Tensor<T> mary = inner(fred,dave,-1,-1);

    ITERATOR2(mary,
              T sum = 0;
              for (int k=0; k<fred.dim(1); ++k) sum += fred(_i,k)*dave(_j,k);
              if (std::abs(sum-mary(_i,_j))> std::abs(sum)*1e-6) error("test7: failed",2));
    dave = dave.swapdim(0,1);
    mary = inner(fred,dave);
    ITERATOR2(mary,
              T sum = 0;
              for (int k=0; k<fred.dim(1); ++k) sum += fred(_i,k)*dave(k,_j);
              if (std::abs(sum-mary(_i,_j))> std::abs(sum)*1e-6) error("test7: failed",3));

    Tensor<T> alfred(3,4,5,6),samantha(7,5,8,2);
    alfred.fillrandom();
    samantha.fillrandom();
    Tensor<T> john = inner(alfred,samantha,2,1);
    ITERATOR6(john,
              T sum = 0;
              for (int k=0; k<alfred.dim(2); ++k) sum += alfred(_i,_j,k,_k)*samantha(_l,k,_m,_n);
              if (std::abs(sum-john(_i,_j,_k,_l,_m,_n))> std::abs(sum)*1e-6) error("test7: failed",4));


    alfred   = Tensor<T>(6,5,4,6);
    samantha = Tensor<T>(6,7,8,6);
    alfred.fillrandom();
    samantha.fillrandom();

    john = inner(alfred,samantha,0,0);
    ITERATOR6(john,
              T sum = 0;
              for (int k=0; k<alfred.dim(0); ++k) sum += alfred(k,_i,_j,_k)*samantha(k,_l,_m,_n);
    if (std::abs(sum-john(_i,_j,_k,_l,_m,_n))> std::abs(sum)*1e-6) {
        std::cout << sum << " " << john(_i,_j,_k,_l,_m,_n) << " " << sum - john(_i,_j,_k,_l,_m,_n) << std::endl;
            error("test7: failed",410);
        }
             );

    john = inner(alfred,samantha,0,-1);
    ITERATOR6(john,
              T sum = 0;
              for (int k=0; k<alfred.dim(0); ++k) sum += alfred(k,_i,_j,_k)*samantha(_l,_m,_n,k);
              if (std::abs(sum-john(_i,_j,_k,_l,_m,_n))> std::abs(sum)*1e-6) error("test7: failed",411));

    john = inner(alfred,samantha,-1,0);
    ITERATOR6(john,
              T sum = 0;
              for (int k=0; k<alfred.dim(0); ++k) sum += alfred(_i,_j,_k,k)*samantha(k,_l,_m,_n);
              if (std::abs(sum-john(_i,_j,_k,_l,_m,_n))> std::abs(sum)*1e-6) error("test7: failed",412));

    john = inner(alfred,samantha,-1,-1);
    ITERATOR6(john,
              T sum = 0;
              for (int k=0; k<alfred.dim(0); ++k) sum += alfred(_i,_j,_k,k)*samantha(_l,_m,_n,k);
              if (std::abs(sum-john(_i,_j,_k,_l,_m,_n))> std::abs(sum)*1e-6) error("test7: failed",413));


    const long n = 8; // Was 5 but changed to even value to test fast_transform
    Tensor<T> x(n,n,n,n);
    Tensor<T> r(n,n,n,n);
    Tensor<T> c(n,n);
    x.fillrandom();
    c.fillrandom();
    x -= (T) 0.5;
    c -= (T) 0.5;

    std::cout << c;
    Tensor<T> q = inner(x,c,0,0);
    Tensor<T> qq(n,n,n,n);
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            for (int k=0; k<n; ++k) {
                for (int l=0; l<n; ++l) {
                    T sum = 0;
                    for (int ii=0; ii<n; ++ii) {
                        sum += x(ii,j,k,l)*c(ii,i);
                    }
                    if (!check(sum,q(j,k,l,i),1e-5)) {
                        std::cout << sum << " " << q(j,k,l,i) << " " << sum-q(j,k,l,i) << std::endl;
                        error("test7: failed",41);
                    }
                }
            }
        }
    }

    Tensor<T> y = transform(x,c);
    r.fill(0);
    // stupid n^8 loop (rather than n^5 optimal) to ensure correctness
    for (int i=0; i<n; ++i) {
        for (int ii=0; ii<n; ++ii) {
            for (int j=0; j<n; ++j) {
                for (int jj=0; jj<n; ++jj) {
                    for (int k=0; k<n; ++k) {

                        for (int kk=0; kk<n; ++kk) {
                            for (int l=0; l<n; ++l) {
                                for (int ll=0; ll<n; ++ll) {
                                    r(i,j,k,l) += x(ii,jj,kk,ll)*c(ii,i)*c(jj,j)*c(kk,k)*c(ll,l);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if ((r-y).normf() > 1e-5*r.normf()) error("test7: failed",5);

    r = copy(x);
    Tensor<T> workspace = copy(x);
    fast_transform(x,c,r,workspace);
    if ((r-y).normf() > 1e-6*r.normf()) {
        std::cout << "R\n";
        std::cout << r;
        std::cout << "Y\n";
        std::cout << y;
        error("test7: failed",555);
    }

    x = Tensor<T>(n,n,n,n,n);	// odd # dimensions
    x.fillrandom();
    x -= (T) 0.5;
    y = transform(x,c);
    workspace = copy(x);
    r = copy(x);
    fast_transform(x,c,r,workspace);
    if ((r-y).normf() > 1e-6*r.normf()) error("test7: failed",666);


    r = Tensor<T>(7,9);
    x = Tensor<T>(9);
    y = Tensor<T>(7);
    r.fillrandom();
    x.fillrandom();

    double tol = r.normf()*x.normf()*1e-6;

    if ((inner(r,x)-inner(x,r.swapdim(0,1))).normf() > tol)
        error("test7: failed",6);

    c = inner(r,x);
    for (int i=0; i<7; ++i) {
        T sum = 0;
        for (int j=0; j<9; ++j) sum += r(i,j)*x(j);
        //std::cout << i << " -> " << sum << " " << c(i) << std::endl;
        if (std::abs(sum-c(i)) > tol) error("test7: failed",7);
    }


    long NNN=24;
    r = Tensor<T>(NNN,NNN,NNN); // result
    y = Tensor<T>(NNN,NNN,NNN); // workspace
    x = Tensor<T>(NNN,NNN,NNN); // input
    c = Tensor<T>(NNN,NNN);
    x.fillrandom();
    c.fillrandom();
    double start = std::clock();
    for (long i=0; i<10000; ++i) {
        fast_transform(x,c,r,y);
    }
    double used = (std::clock()-start)/CLOCKS_PER_SEC;
    double mops = 10000.0*2*3*1e-6*double(NNN)*double(NNN)*double(NNN)*NNN/used;
    std::cout << "TRANSFORM MOPS=" << mops << "   " << used << std::endl;

    std::cout << "Test7<" << tensor_type_names[TensorTypeData<T>::id] << "> OK\n";
}

template <class T> void Test8() {
    // test complex only operations and type restrictions
}

int main() {

    std::cout << "bounds checking enabled = " << Tensor<double>::bounds_checking()
              << std::endl;

    Tensor<double> a(2);
    std::cout << "made A\n";
    try {
        std::cout << "testing access\n";
        std::cout << a(2) << std::endl;
        if (a.bounds_checking()) error("bounds check failed to detect error",1);
    }
    catch (TensorException e) {
        std::cout << "If bounds checking, this exception was expected" << std::endl;
        std::cout << e;
        if (!a.bounds_checking())
            error("bounds checking was not enabled",0);
    }
    try {
        std::cout << a(-1) << std::endl;
        if (a.bounds_checking()) error("bounds check failed to detect error",2);
    }
    catch (TensorException e) {
        std::cout << "If bounds checking, this exception was expected" << std::endl;
        std::cout << e;
        if (!a.bounds_checking())
            error("bounds checking was not enabled",0);
    }

    Test1<double,double>();
    Test1<double,long>();
    Test1<long,double>();
    Test1<double,float>();
    Test1<float,double>();
    Test1<double_complex,double>();
    Test1<double_complex,long>();
    Test1<float_complex,float>();
    //Test1<double,double_complex>(); // Correctly fails to compile

    std::cout << "\n after test1 count=" <<
              Tensor<long>().get_instance_count() << std::endl;

    std::cout << std::endl;
    Test3<double>();
    Test3<long>();
    Test3<float>();
    Test3<float_complex>();
    Test3<double_complex>();
    std::cout << std::endl;
    std::cout << "\n after test3 count=" <<
              Tensor<long>().get_instance_count() << std::endl;

    std::cout << std::endl;

    Test5<double>();
    Test5<long>();
    //Test5<float>();
    //Test5<float_complex>(); // float * double -> mess
    std::cout << "skipping Test5<float/float_complex>" << std::endl;
    Test5<double_complex>();
    std::cout << std::endl;
    std::cout << "\n after test5 count=" <<
              Tensor<long>().get_instance_count() << std::endl;

    std::cout << std::endl;
    Test6<double>();
    Test6<long>();
    Test6<float>();
    Test6<float_complex>();
    Test6<double_complex>();
    std::cout << std::endl;
    std::cout << "\n after test6 count=" <<
              Tensor<long>().get_instance_count() << std::endl;

    std::cout << std::endl;
    Test7<double>();
    Test7<float>();
    Test7<float_complex>();
    Test7<double_complex>();
    std::cout << std::endl;

    std::cout << "\n after tests count=" <<
              Tensor<long>().get_instance_count() << std::endl;

    Tensor<double_complex> f(1);
    Tensor<double> g = real(f);
    // real(g);  // will not link ... since we did not instantiate invalid operations

    std::cout << "\n\n Tensor seems vaguely OK ... fingers crossed\n";

    return 0;
}
