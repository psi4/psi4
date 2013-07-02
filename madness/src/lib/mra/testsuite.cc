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


  $Id: test.cc 257 2007-06-25 19:09:38Z HartmanBaker $
*/

/// \file testsuite.cc
/// \brief The QA/test suite for Function

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <unistd.h>
#include <cstdio>
#include <constants.h>
#include <mra/qmprop.h>

#include <misc/ran.h>

const double PI = 3.1415926535897932384;

using namespace madness;

template <typename T, std::size_t NDIM>
struct lbcost {
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        return 1.0;
    }
};

template <typename T>
T complexify(T c) {
    return c;
}

template <> double_complex complexify<double_complex>(double_complex c) {
    return c*double_complex(0.5,-sqrt(3.0)*0.5);
}

template <> float_complex complexify<float_complex>(float_complex c) {
    return c*float_complex(0.5,-sqrt(3.0)*0.5);
}


template <typename T, std::size_t NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;

    Gaussian(const coordT& center, double exponent, T coefficient)
            : center(center), exponent(exponent), coefficient(complexify(coefficient)) {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    };
};

template <typename T, std::size_t NDIM>
class DerivativeGaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;
    const int axis;

    DerivativeGaussian(const coordT& center, double exponent, T coefficient, int axis)
            : center(center), exponent(exponent), coefficient(complexify(coefficient)), axis(axis) {};

    DerivativeGaussian(const Gaussian<T,NDIM>& g, int axis)
            : center(g.center), exponent(g.exponent), coefficient(g.coefficient), axis(axis) {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return -2.0*exponent*(x[axis]-center[axis])*coefficient*exp(-exponent*sum);
    };
};


template <typename T, typename L, typename R>
inline T product(L l, R r) {
    return T(l*r);
}

template <typename T, typename L, typename R>
inline T sum(L l, R r) {
    return T(l+r);
}


/// Makes a square-normalized Gaussian with random origin and exponent
template <typename T, std::size_t NDIM>
Gaussian<T,NDIM>*
RandomGaussian(const Tensor<double> cell, double expntmax=1e5) {
    typedef Vector<double,NDIM> coordT;
    coordT origin;
    for (std::size_t i=0; i<NDIM; ++i) {
        origin[i] = RandomValue<double>()*(cell(i,1)-cell(i,0)) + cell(i,0);
    }
    double lo = log(0.1);
    double hi = log(expntmax);
    double expnt = exp(RandomValue<double>()*(hi-lo) + lo);
    T coeff = pow(2.0*expnt/PI,0.25*NDIM);
    //print("RandomGaussian: origin", origin, "expnt", expnt, "coeff", coeff);
    return new Gaussian<T,NDIM>(origin,expnt,coeff);
}

/// Returns a new functor combining two functors via operation op(left,right)
template <typename resultT, typename L, typename R, typename opT, std::size_t NDIM>
class BinaryOp : public FunctionFunctorInterface<resultT,NDIM> {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<L,NDIM> > functorL;
    typedef std::shared_ptr< FunctionFunctorInterface<R,NDIM> > functorR;

    functorL left;
    functorR right;
    opT op;

public:
    BinaryOp(const functorL& left, const functorR& right, opT& op)
            : left(left), right(right), op(op) {};

    resultT operator()(const coordT& x) const {
        return op((*left)(x),(*right)(x));
    };
};

#define CHECK(value, threshold, message)        \
        do { \
             if (world.rank() == 0) { \
                 bool status = std::abs(value) < threshold;             \
                const char* msgs[2] = {"FAILED","OK"};                  \
                std::printf("%20.20s :%5d :%30.30s : %10.2e  < %10.2e : %s\n", (__FUNCTION__),(__LINE__),(message),(std::abs(value)),(threshold), (msgs[int(status)])); \
                if (!status) ok = false; \
             } \
        } while (0)

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt)


template <typename T, std::size_t NDIM>
void test_basic(World& world) {
    bool ok = true;
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0)
        print("Test compression of a normalized gaussian at origin, type =",
              archive::get_type_name<T>(),", ndim =",NDIM);

    Tensor<double> cell(NDIM,2);
    for (std::size_t i=0; i<NDIM; ++i) {
        cell(i,0) = -11.0-2*i;  // Deliberately asymmetric bounding box
        cell(i,1) =  10.0+i;
    }
    FunctionDefaults<NDIM>::set_cell(cell);
    FunctionDefaults<NDIM>::set_k(7);
    FunctionDefaults<NDIM>::set_thresh(1e-5);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);

    const coordT origin(0.0);
    coordT point;
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);

    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    for (std::size_t i=0; i<NDIM; ++i) point[i] = 0.1*i;

    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);

    double norm = f.norm2();
    double err = f.err(*functor);
    T val = f(point);
    CHECK(std::abs(norm-1.0), 1.0e-10, "norm");
    CHECK(err, 3e-7, "err");
    CHECK(val-(*functor)(point), 1e-8, "error at a point");

    f.compress();
    double new_norm = f.norm2();
    CHECK(new_norm-norm, 1e-14, "new_norm");

    f.reconstruct();
    new_norm = f.norm2();
    double new_err = f.err(*functor);
    CHECK(new_norm-norm, 1e-14, "new_norm");
    CHECK(new_err-err, 1e-14, "new_err");

    f.compress();
    new_norm = f.norm2();
    CHECK(new_norm-norm, 1e-14, "new_norm");

    f.truncate();
    new_norm = f.norm2();
    new_err = f.err(*functor);
    CHECK(new_norm-norm, 1e-9, "new_norm");
    CHECK(new_err, 3e-5, "new_err");

    world.gop.fence();
    if (world.rank() == 0) print("projection, compression, reconstruction, truncation OK\n\n");
}

template <typename T, std::size_t NDIM>
void test_conv(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("Test convergence - log(err)/(n*k) should be roughly const, a least for each value of k");
        print("                 - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }
    const coordT origin(0.0);
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    FunctionDefaults<NDIM>::set_cubic_cell(-10,10);

    for (int k=2; k<=30; k+=2) {
        if (world.rank() == 0) printf("k=%d\n", k);
        int ntop = 5;
        if (NDIM > 2 && k>5) ntop = 4;
        for (int n=1; n<=ntop; ++n) {
            Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor).norefine().initial_level(n).k(k);
            double err2 = f.err(*functor);
            std::size_t size = f.size();
            if (world.rank() == 0)
                printf("   n=%d err=%.2e #coeff=%.2e log(err)/(n*k)=%.2e\n",
                       n, err2, double(size), std::abs(log(err2)/n/k));
        }
    }

    world.gop.fence();
    if (world.rank() == 0) print("test conv OK\n\n");
}

template <typename T, std::size_t NDIM>
struct myunaryop {
    typedef T resultT;
    Tensor<T> operator()(const Key<NDIM>& key, const Tensor<T>& t) const {
        return -t;
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};

template <typename T, std::size_t NDIM>
struct myunaryop_square {
    typedef T resultT;
    Tensor<T> operator()(const Key<NDIM>& key, const Tensor<T>& t) const {
        Tensor<T> result = copy(t);
        T* r = result.ptr();
        for (int i = 0; i < result.size(); ++i) {
            r[i] = r[i]*r[i];
        }
        return result;
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};

template <typename T, std::size_t NDIM>
void test_math(World& world) {
    bool ok = true;
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("Test basic math operations - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }

    FunctionDefaults<NDIM>::set_k(9);
    FunctionDefaults<NDIM>::set_thresh(1e-9);
    FunctionDefaults<NDIM>::set_truncate_mode(0);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_autorefine(false);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_cubic_cell(-10,10);

    const coordT origin(0.0);
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));
    T(*p)(T,T) = &product<T,T,T>;
    functorT functsq(new BinaryOp<T,T,T,T(*)(T,T),NDIM>(functor,functor,p));
    //functorT functsq(new Gaussian<T,NDIM>(origin, 2.0*expnt, coeff*coeff)); // only correct if real

    // First make sure out of place squaring works
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);

    // Out-of-place squaring without autoref
    f.set_autorefine(false);
    world.gop.fence();

    double err = f.err(*functor);
    CHECK(err, 1e-10, "err in f before squaring");
    Function<T,NDIM> fsq = square(f);
    double new_err = f.err(*functor);
    CHECK(new_err-err,1e-14*err,"err in f after squaring");
    double errsq = fsq.err(*functsq);
    CHECK(errsq, 1e-10, "err in fsq");

    // Test same with autorefine
    fsq = unary_op(f, myunaryop<T,NDIM>());
    errsq = (f + fsq).norm2();
    CHECK(errsq, 1e-10, "err in unaryp_op negate");
    fsq.clear();
    f.reconstruct();

    // Test same with autorefine
    fsq = unary_op(f, myunaryop_square<T,NDIM>());
    errsq = (fsq - square(f)).norm2();
    CHECK(errsq, 1e-10, "err in unaryp_op square");
    fsq.clear();
    f.reconstruct();

//     // Test same with autorefine
//     f.set_autorefine(true); world.gop.fence();
//     fsq = square(f);
//     errsq = fsq.err(*functsq);
//     CHECK(errsq, 1e-10, "err in fsq with autoref");

//     // Repeat after agressive truncating to see if autorefine really works
//     f.set_autorefine(false); world.gop.fence();
//     f.truncate(1e-5);
//     f.verify_tree();
//     err = f.err(*functor);
//     CHECK(err, 1e-5, "error in f after truncating");

//     fsq = square(f);
//     errsq = fsq.err(*functsq);
//     CHECK(errsq, 1e-5, "error in fsq after truncating");

//     f.set_autorefine(true); world.gop.fence();
//     fsq = square(f);
//     errsq = fsq.err(*functsq);
//     CHECK(errsq, 1e-5, "error in fsq truncate+autoref");

//     // Finally inplace squaring
//     f.square();
//     double new_errsq = f.err(*functsq);
//     CHECK(new_errsq - errsq, 1e-14*errsq, "err in fsq trunc+auto+inplace");

    fsq.clear();

    // Test adding a constant in scaling function and wavelet bases
    T val = f(origin);
    f.reconstruct();
    f.add_scalar(3.0);
    T val2 = f(origin);

    f.compress();
    f.add_scalar(5.0);
    f.reconstruct();
    val2 = f(origin);
    CHECK(val2-(val+8.0),1e-12,"add scalar in place compressed");

    // Test in-place scaling by a constant in scaling function and wavelet bases
    f.reconstruct();
    f.scale(3.0);
    val2 = f(origin);
    CHECK(val2-3.0*(val+8.0),1e-12,"in-place scaling reconstructed");

    f.compress();
    f.scale(4.0);
    f.reconstruct();
    val2 = f(origin);
    CHECK(val2-12.0*(val+8.0),1e-12,"in-place scaling compressed");

    // Same but using operator notation
    f.reconstruct();
    f *= 7.0;
    val2 = f(origin);
    CHECK(val2-7.0*12.0*(val+8.0),1e-11,"in-place scaling (op) recon");

    f.compress();
    f *= 7.0;
    f.reconstruct();
    val2 = f(origin);
    CHECK(val2-7.0*7.0*12.0*(val+8.0),1e-10,"in-place scaling (op) comp");

    // Test squaring a function by multiplication
    f = Function<T,NDIM>(FunctionFactory<T,NDIM>(world).functor(functor));
    fsq = f*f;
    errsq = fsq.err(*functsq);
    CHECK(errsq, 1e-8, "err in fsq by multiplication");

    // Test norm tree operation
    f.reconstruct();
    f.norm_tree();
    double nnn = f.norm2();
    if (world.rank() == 0) print("nnn", nnn);

    // Test composing operations using general expression(s)
    f.compress();
    err = f.err(*functor);
    f.compress();
    Function<T,NDIM> f6 = f*3.0 + 4.0*f - f;
    new_err = f.err(*functor);
    CHECK(new_err-err,1e-14,"general op unchanged input");
    new_err = (f6 - f.scale(6.0)).norm2();
    CHECK(new_err,1e-13,"general op output");

    if (world.rank() == 0) print("\nTest multiplying random functions");
    default_random_generator.setstate(314159);  // Ensure all processes have the same sequence (for exponents)

    FunctionDefaults<NDIM>::set_autorefine(false);

    int nfunc = 100;
    if (NDIM >= 3) nfunc = 20;
    for (int i=0; i<nfunc; ++i) {
        functorT f1(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),100.0));
        functorT f2(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),100.0));
        T(*p)(T,T) = &product<T,T,T>;
        functorT f3(new BinaryOp<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
        Function<T,NDIM> a = FunctionFactory<T,NDIM>(world).functor(f1);
        a.verify_tree();
        Function<T,NDIM> b = FunctionFactory<T,NDIM>(world).functor(f2);
        b.verify_tree();
        //print("NORMS", a.norm2(), b.norm2());
        //std::cout.flush();
        Function<T,NDIM> c = a*b;
        double err1 = a.err(*f1);
        double err2 = b.err(*f2);
        double err3 = c.err(*f3);
        if (world.rank() == 0) print("  test ",i);
        CHECK(err1,1e-8,"err1");
        CHECK(err2,1e-8,"err2");
        CHECK(err3,1e-8,"err3");

//         double bnorm = b.norm2();
//         if (world.rank() == 0) print("bnorm", bnorm);
//         b.norm_tree();
//         print("++++++++++++++++++++++++++++");
//         Function<T,NDIM> cs = mul_sparse(a,b,1e-4);
//         print("----------------------------");
//         cs.verify_tree();
//         if (world.rank() == 0) print("cs - c", (cs-c).norm2());

    }

    if (world.rank() == 0) print("\nTest multiplying a vector of random functions");
    {
        functorT f1(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),1000.0));
        Function<T,NDIM> left = FunctionFactory<T,NDIM>(world).functor(f1);

        const int nvfunc = 10;
        std::vector< Function<T,NDIM> > vin(nvfunc);
        std::vector<functorT> funcres(nvfunc);
        for (int i=0; i<nvfunc; ++i) {
            functorT f2(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),1000.0));
            T(*p)(T,T) = &product<T,T,T>;
            funcres[i] = functorT(new BinaryOp<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
            vin[i] = FunctionFactory<T,NDIM>(world).functor(f2);
        }
        std::vector< Function<T,NDIM> > vres = mul(world, left, vin);
        for (int i=0; i<nvfunc; ++i) {
            double err = vres[i].err(*funcres[i]);
            CHECK(err, 1e-8, "err");
            vres[i].verify_tree();
        }
    }



    if (world.rank() == 0) print("\nTest adding random functions out of place");
    for (int i=0; i<10; ++i) {
        functorT f1(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),100.0));
        functorT f2(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),100.0));
        T(*p)(T,T) = &sum<T,T,T>;
        functorT f3(new BinaryOp<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
        Function<T,NDIM> a = FunctionFactory<T,NDIM>(world).functor(f1);
        Function<T,NDIM> b = FunctionFactory<T,NDIM>(world).functor(f2);
        Function<T,NDIM> c = a+b;
        a.verify_tree();
        b.verify_tree();
        c.verify_tree();
        double err1 = a.err(*f1);
        double err2 = b.err(*f2);
        double err3 = c.err(*f3);
        if (world.rank() == 0) print("  test ",i);
        CHECK(err1,1e-8,"err1");
        CHECK(err2,1e-8,"err2");
        CHECK(err3,1e-8,"err3");
    }

    if (world.rank() == 0) print("\nTest adding random functions in place");
    for (int i=0; i<10; ++i) {
        functorT f1(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),100.0));
        functorT f2(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),100.0));
        T(*p)(T,T) = &sum<T,T,T>;
        functorT f3(new BinaryOp<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
        Function<T,NDIM> a = FunctionFactory<T,NDIM>(world).functor(f1);
        Function<T,NDIM> b = FunctionFactory<T,NDIM>(world).functor(f2);
        a.verify_tree();
        b.verify_tree();
        a += b;
        double err1 = a.err(*f3);
        double err2 = b.err(*f2);
        if (world.rank() == 0) print("  test ",i);
        CHECK(err1,1e-8,"err1");
        CHECK(err2,1e-8,"err2");
    }

    world.gop.fence();
}


template <typename T, std::size_t NDIM>
void test_diff(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("\nTest differentiation - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }
    const coordT origin(0.0);
    //for (int i=0; i<NDIM; ++i) origin[i] = i/31.4;
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    FunctionDefaults<NDIM>::set_k(10);
    FunctionDefaults<NDIM>::set_thresh(1e-10);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);
    FunctionDefaults<NDIM>::set_truncate_mode(1);
    FunctionDefaults<NDIM>::set_cubic_cell(-10,10);

    START_TIMER;
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);
    END_TIMER("project");

    //f.print_info();  <--------- This is not scalable and might crash the XT

    START_TIMER;
    f.compress();
    END_TIMER("compress");

    START_TIMER;
    f.truncate();
    END_TIMER("truncate");

    START_TIMER;
    f.reconstruct();
    END_TIMER("reconstruct");


    for (std::size_t axis=0; axis<NDIM; ++axis) {
        if (world.rank() == 0) print("doing axis", axis);
        Derivative<T,NDIM> D(world, axis);

        DerivativeGaussian<T,NDIM> df(origin,expnt,coeff,axis);

        START_TIMER;
        Function<T,NDIM> dfdx = D(f);
        END_TIMER("diff");

//         coordT p(0.0);
//         if (world.rank() == 0) {
//             for (int i=0; i<=40; ++i) {
//                 p[axis] = (i-20.0)*0.1;
//                 print("     x, analytic, err",p[axis],df(p), dfdx(p)-df(p));
//             }
//         }
//         world.gop.fence();

        START_TIMER;
        double err = dfdx.err(df);
        END_TIMER("err");

        if (world.rank() == 0) print("    error", err);
    }
    world.gop.fence();
}



namespace madness {
    extern bool test_rnlp();
}

template <typename T, std::size_t NDIM>
void test_op(World& world) {
    test_rnlp();

    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("\nTest separated operators - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }
    const coordT origin(0.5);
    const double expnt = 1.0*100;
    const double coeff = pow(2.0*expnt/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    FunctionDefaults<NDIM>::set_k(10);
    FunctionDefaults<NDIM>::set_thresh(1e-12);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);
    FunctionDefaults<NDIM>::set_truncate_mode(1);
    FunctionDefaults<NDIM>::set_cubic_cell(-10,10);

    START_TIMER;
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);
    END_TIMER("project");

    //f.print_info();  <--------- This is not scalable and might crash the XT


    f.reconstruct();
    double n2 = f.norm2();
    double e2 = f.err(*functor);
    if (world.rank() == 0) {
        print("         f norm is", n2);
        print("     f total error", e2);
    }

//     f.reconstruct();
//     Function<T,NDIM> fff = copy(f);
//     for (int i=0; i<10; ++i) {
//         fff.compress().reconstruct();
//     }
//     f.compress();
//     double ecr = (fff-f).norm2();
//     if (world.rank() == 0) print("error after 10 compress-reconstruct",ecr);

//     fff.reconstruct();
//     for (int i=0; i<10; ++i) {
//         fff.nonstandard(false,true);
//         fff.standard();
//         fff.reconstruct();
//     }
//     fff.compress();
//     ecr = (fff-f).norm2();
//     if (world.rank() == 0) print("error after 10 non-standard compress-reconstruct",ecr);


    // Convolution exp(-a*x^2) with exp(-b*x^2) is
    // exp(-x^2*a*b/(a+b))* (Pi/(a+b))^(NDIM/2)

    Tensor<double> coeffs(1), exponents(1);
    exponents(0L) = 10.0;
    coeffs(0L) = pow(exponents(0L)/PI, 0.5*NDIM);
    SeparatedConvolution<T,NDIM> op(world, coeffs, exponents);
    START_TIMER;
    Function<T,NDIM> r = apply(op,f);
    END_TIMER("apply");
    r.verify_tree();
    f.verify_tree();

    double newexpnt = expnt*exponents(0L)/(expnt+exponents(0L));
    double newcoeff = pow(PI/(expnt+exponents(0L)),0.5*NDIM)*coeff*coeffs(0L);
    functorT fexact(new Gaussian<T,NDIM>(origin, newexpnt, newcoeff));

    T ro = r(origin);
    T eo = (*fexact)(origin);
    double rn = r.norm2();
    double re = r.err(*fexact);
    if (world.rank() == 0) {
        print(" numeric at origin", ro);
        print("analytic at origin", eo);
        print("      op*f norm is", rn);
        print("  op*f total error", re);
    }
//     for (int i=0; i<=100; ++i) {
//         coordT c(-10.0+20.0*i/100.0);
//         print("           ",i,c[0],r(c),r(c)-(*fexact)(c));
//     }
}

/// Computes the electrostatic potential due to a Gaussian charge distribution
class GaussianPotential : public FunctionFunctorInterface<double,3> {
public:
    typedef Vector<double,3> coordT;
    const coordT center;
    const double exponent;
    const double coefficient;

    GaussianPotential(const coordT& center, double expnt, double coefficient)
            : center(center)
            , exponent(sqrt(expnt))
            , coefficient(coefficient*pow(PI/exponent,1.5)*pow(expnt,-0.75)) {}

    double operator()(const coordT& x) const {
        double sum = 00;
        for (int i=0; i<3; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        double r = sqrt(sum);
        return coefficient*erf(exponent*r)/r;
    }
};

void test_coulomb(World& world) {
    typedef Vector<double,3> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

    if (world.rank() == 0) {
        print("\nTest Coulomb operator - type =", archive::get_type_name<double>(),", ndim = 3 (only)\n");
    }

    // Normalized Gaussian exponent a produces potential erf(sqrt(a)*r)/r
    const coordT origin(0.5);
    const double expnt = 100.0;
    const double coeff = pow(1.0/PI*expnt,0.5*3);
    functorT functor(new Gaussian<double,3>(origin, expnt, coeff));


    int k = 10;
    double thresh = 1e-8;

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_truncate_mode(0);
    FunctionDefaults<3>::set_cubic_cell(-10,10);

    START_TIMER;
    Function<double,3> f = FunctionFactory<double,3>(world).functor(functor).thresh(thresh*0.1).initial_level(4);
    END_TIMER("project");

    //f.print_info();  <--------- This is not scalable and might crash the XT

    f.reconstruct();
    double norm = f.norm2(), err = f.err(*functor);
    if (world.rank() == 0) {
        print("         f norm is", norm);
        print("     f total error", err);
        //print(" truncating");
    }

//     START_TIMER;
//     f.truncate();
//     END_TIMER("truncate");
//     START_TIMER;
//     f.reconstruct();
//     END_TIMER("reconstruct");
//     norm = f.norm2();
//     err = f.err(*functor);
//     if (world.rank() == 0) {
//         print("         f norm is", norm);
//         print("     f total error", err);
//     }

    f.reconstruct();
    START_TIMER;
    f.nonstandard(false,true);
    END_TIMER("nonstandard");

    if (world.rank() == 0) {
        print("\nbefore operator");
        //world.am.print_stats();
        print("");
    }
    f.set_thresh(thresh);
    SeparatedConvolution<double,3> op = CoulombOperator(world, 1e-5, thresh);

    FunctionDefaults<3>::set_apply_randomize(true);

    START_TIMER;
    Function<double,3> r = apply_only(op,f) ;
    END_TIMER("apply");


    START_TIMER;
    r.reconstruct();
    END_TIMER("reconstruct result");
    r.verify_tree();

    functorT fexact(new GaussianPotential(origin, expnt, coeff));

    double numeric=r(origin);
    double analytic=(*fexact)(origin);
    double rnorm = r.norm2();
    double rerr = r.err(*fexact);
    if (world.rank() == 0) {
        print(" numeric at origin", numeric);
        print("analytic at origin", analytic);
        print("      op*f norm is", rnorm);
        print("  op*f total error", rerr);
//         for (int i=0; i<=100; ++i) {
//             coordT c(i*0.01);
//             print("           ",i,r(c),(*fexact)(c));
//         }
    }
}

class QMtest : public FunctionFunctorInterface<double_complex,1> {
public:
    typedef Vector<double,1> coordT;
    const double a;
    const double v;
    const double t;

    double_complex operator()(const coordT& coords) const {
        const double x = coords[0];
        double_complex denom(1.0,2.0*a*t);
        const double_complex arg(a*x*x,0.5*v*v*t-x*v);
        return pow(2.0*a/PI,0.25)*sqrt(1.0/denom)*exp(-arg/denom);
    }

    QMtest(double a, double v, double t)
            : a(a), v(v), t(t) {}
};

    struct refop {
        bool operator()(FunctionImpl<double_complex,1>* impl, const Key<1>& key, const Tensor<double_complex>& t) const {
            double tol = impl->truncate_tol(impl->get_thresh(), key);
            double lo, hi;
            impl->tnorm(t, &lo, &hi);
            return hi > tol;;
        }
        template <typename Archive> void serialize(Archive& ar) {}
    };


void test_qm(World& world) {
    /*

      This is the exact kernel of the free-particle propagator

      g(x,t) = exp(I*x^2/(2*t))/sqrt((2*Pi*I)*t)

      This is a square normalized Gaussian (in 1D) with velocity v.

      f(x,a,v) = (2*a/Pi)^(1/4)*exp(-a*x^2+I*x*v)

      This is f(x,a,v) evolved to time t

      f(x,a,v,t) = (2*a/Pi)^(1/4)*sqrt(1/(1+(2*I)*a*t))*exp(-(a*x^2-I*x*v+I*v^2*t*1/2)/(1+(2*I)*a*t))

      The fourier transform of f(x,a,v) decays as exp(-k^2/(4*a))
      and is neglible when k > v + 10*sqrt(a)

      Picking a=v=1 gives an effective bandlimit of c=11.  If we wish this
      to be accurately propagated for a long time we must request a filtered
      bandlimit of c = 11*1.8 = 20.

      The center of the packet is at <x> = v*t.

      The width of the packet as measured by sigma^2 = <x^2 - <x>^2> = (1/4)*(1+4*a^2*t^2)/a.

      So the wave packet is actually spreading faster than it is moving (for our choice
      of parameters).

      The initial support of the Gaussian is about [-5,5] and after 100 units of
      time it has spread to about [-400,600] so for safety we use a range [-600,800].

      The critical time step is 2*pi/c^2= 0.0157 and we shall attempt to propagate at
      10x this which is 0.157.

      The final wave packet is horrible, oscillating thru all space due to the
      complex phase ... this could be reduced with a contact (?) transformation but
      since the point here is test MADNESS to some extent it is better to make things hard.

    */

//     QMtest f(1,1,0.1);

//     for (int i=0; i<10; ++i) {
//         double x = i*0.1;
//         print(x,f(x));
//     }

    typedef std::shared_ptr< FunctionFunctorInterface<double_complex,1> > functorT;
    typedef Vector<double,1> coordT;
    typedef Function<double_complex,1> functionT;
    typedef FunctionFactory<double_complex,1> factoryT;

    //int k = 16;
    //double thresh = 1e-12;

    int k = 16;
    double thresh = 1e-12;
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_refine(true);
    FunctionDefaults<1>::set_initial_level(8);
    FunctionDefaults<1>::set_cubic_cell(-600,800);
    FunctionDefaults<1>::set_truncate_mode(0);
    double width = FunctionDefaults<1>::get_cell_width()(0L);

    double a = 1.0;
    double v = 1.0;
    double ctarget = v + 10.0*sqrt(a);
    double c = 1.86*ctarget; //1.86*ctarget;
    double tcrit = 2*PI/(c*c);
    double tstep = 3*tcrit;

    int nstep = int(100.0/tstep);
    tstep = 100.0/nstep; // so we finish exactly at 100.0

    // For the purpose of testing there is no need to propagate 100 time units.
    // Just 100 steps.
    nstep = 100;

    if (world.rank() == 0) {
        print("\n Testing evolution of a quantum wave packet in",1,"dimensions");
        print("expnt",a,"velocity",v,"bandw",ctarget,"effbandw",c);
        print("tcrit",tcrit,"tstep",tstep,"nstep",nstep,"width",width);
    }

    functorT f(new QMtest(a,v,0.0));
    //SeparatedConvolution<double_complex,1> G = qm_free_particle_propagator<1>(world, k, c, tstep);
    //G.doleaves = true;

    functionT psi = factoryT(world).functor(f).initial_level(12);
    psi.truncate();

    if (world.rank() == 0) {
        print("  step    time      norm      error");
        print(" ------  ------- ---------- ----------");
    }


    Convolution1D<double_complex>* q1d = qm_1d_free_particle_propagator(k, c, tstep, 1400.0);

    for (int i=0; i<nstep; ++i) {
        world.gop.fence();

        psi.reconstruct();
        //psi.refine_general(refop());
        psi.broaden();
        psi.broaden();

        world.gop.fence();
        double norm = psi.norm2();
        double err = psi.err(QMtest(a,v,tstep*i));
        if (world.rank() == 0)
            printf("%6d  %7.3f  %10.8f  %9.1e\n",i, i*tstep, norm, err);

        //         print("psi");
        //         psi.print_tree();

        functionT pp = apply_1d_realspace_push(*q1d, psi, 0);

        //         print("pp before sum down");
        //         pp.print_tree();

        pp.sum_down();

        //         print("pp after sum down");
        //         pp.print_tree();

        //psi.truncate(thresh);
        //psi = apply(G,psi);
        //psi.reconstruct();

        //         print("new psi");
        //         psi.print_tree();

        //double pperr = (pp - psi).norm2();
        //print("ERROR", pperr, pp.norm2());

//         if (pperr > 1e-4) {
//             for (int i=0; i<1001; ++i) {
//                 double x = (i-500)*0.01;
//                 print(x, pp(x), psi(x));
//             }
//             exit(0);
//         }

        psi = pp;

        world.gop.fence();

        psi.truncate();
    }

    // Test program does not need to plot!
//     psi.reconstruct();

//     ofstream plot;

//     if (world.rank() == 0) plot.open("plot.dat",ios::trunc);
//     int npt = 10001;
//     double lo = FunctionDefaults<1>::get_cell()(0,0);
//     double hi = FunctionDefaults<1>::get_cell()(0,1);
//     double h = (hi-lo)/(npt-1);
//     for (int i=0; i<npt; ++i) {
//         double x = lo + i*h;
//         double_complex numeric = psi(x);
//         double_complex exact = QMtest(a,v,tstep*nstep)(x);
//         if (world.rank() == 0) plot << x << " " << numeric.real() << " " << numeric.imag() << " " << std::abs(numeric) << " " << std::abs(numeric-exact) << endl;
//     }
//     if (world.rank() == 0) plot.close();

    return;
}

template <typename T, std::size_t NDIM>
void test_plot(World& world) {
    bool ok = true;
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;
    if (world.rank() == 0) {
        print("\nTest plot cube - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }
    const double L = 4.0;
    FunctionDefaults<NDIM>::set_cubic_cell(-L,L);
    FunctionDefaults<NDIM>::set_k(7);
    FunctionDefaults<NDIM>::set_thresh(1e-7);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);

    const coordT origin(0.6666666);
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);

    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);

    //vector<long> npt(NDIM,21); // recommend this if testing in dimension > 3
    std::vector<long> npt(NDIM,101);
    world.gop.fence();
    Tensor<T> r = f.eval_cube(FunctionDefaults<NDIM>::get_cell(), npt);
    world.gop.fence();
    std::size_t maxlevel = f.max_local_depth();
    if (world.rank() == 0) {
        const double h = (2.0*L - 12e-13)/(npt[0]-1.0);
        for (int i=0; i<npt[0]; ++i) {
            double x = -L + i*h + 2e-13;
            T fplot = r(std::vector<long>(NDIM,i));
            T fnum  = f.eval(coordT(x)).get();
            std::pair<bool,T> fnum2 = f.eval_local_only(coordT(x),maxlevel);

            if (world.size() == 1 && !fnum2.first) print("eval_local_only: non-local but nproc=1!");

            if (fnum2.first) CHECK(fnum-fnum2.second,1e-12,"eval_local_only");
            CHECK(fplot-fnum,1e-12,"plot-eval");
            if (world.rank() == 0 && std::abs(fplot-fnum) > 1e-12) {
                print("bad", i, coordT(x), fplot, fnum, (*functor)(coordT(x)));
            }
        }
    }
    world.gop.fence();

    r = Tensor<T>();
    plotdx(f, "testplot", FunctionDefaults<NDIM>::get_cell(), npt);

    plot_line("testline1", 101, coordT(-L), coordT(L), f);
    plot_line("testline2", 101, coordT(-L), coordT(L), f, f*f);
    plot_line("testline3", 101, coordT(-L), coordT(L), f, f*f, 2.0*f);

    if (world.rank() == 0) print("evaluation of cube/slice for plotting OK");
}

template <typename T, std::size_t NDIM>
void test_io(World& world) {
    if (world.rank() == 0) {
        print("\nTest IO - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }

    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    FunctionDefaults<NDIM>::set_k(5);
    FunctionDefaults<NDIM>::set_thresh(1e-10); // We want lots of boxes
    FunctionDefaults<NDIM>::set_truncate_mode(0);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_cubic_cell(-10,10);

    const coordT origin(0.0);
    const double expnt = 10.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);

    int nio = (world.size()-1)/20 + 1;
    archive::ParallelOutputArchive out(world, "mary", nio);
    out & f;
    out.close();

    Function<T,NDIM> g;

    archive::ParallelInputArchive in(world, "mary", nio);
    in & g;
    in.close();
    in.remove();

    double err = (g-f).norm2();

    if (world.rank() == 0) print("err = ", err);

    //    MADNESS_ASSERT(err == 0.0);

    if (world.rank() == 0) print("test_io OK");
    world.gop.fence();
}

template <typename T, std::size_t NDIM>
void test_apply_push_1d(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0)
        print("Test push1d, type =",archive::get_type_name<T>(),", ndim =",NDIM);

    Tensor<double> cell(NDIM,2);
    const double L = 10.0;
    FunctionDefaults<NDIM>::set_cubic_cell(-L,L);
    FunctionDefaults<NDIM>::set_k(6);
    FunctionDefaults<NDIM>::set_thresh(1e-6);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);

    const coordT origin(0.0);
    const double expnt = 1.0;
    const double coeff = pow(1.0/PI,0.5*NDIM);

    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);

//     f.compress();
//     f.truncate();
//     f.reconstruct();

    double trace = f.trace();
    if (world.rank() == 0)
        print("Trace of f", trace);

    coordT lo(-L), hi(L);
    plot_line("fplot.dat", 201, lo, hi, f);

    GaussianConvolution1D<double> op(6, coeff*2.0*L, expnt*L*L*4.0, 0, false);
    Function<T,NDIM> opf = apply_1d_realspace_push(op, f, 0);

    opf.sum_down();
    trace = opf.trace();
    if (world.rank() == 0)
        print("Trace of opf", trace);
    plot_line("opfplot.dat", 201, lo, hi, opf);

    if (world.rank() == 0) print("result", opf.eval(origin).get());
    world.gop.fence();
}


#define TO_STRING(s) TO_STRING2(s)
#define TO_STRING2(s) #s

int main(int argc, char**argv) {
    initialize(argc, argv);
    try {
        World world(MPI::COMM_WORLD);
        if (world.rank() == 0) {
            print("");
            print("--------------------------------------------");
            print("   MADNESS",PACKAGE_VERSION, "multiresolution testsuite");
            print("--------------------------------------------");
            print("");
            print("   number of processors ...", world.size());
            print("    processor frequency ...", cpu_frequency());
            print("            host system ...", HOST_SYSTEM);
            print("          configured by ...", MADNESS_CONFIGURATION_USER);
            print("          configured on ...", MADNESS_CONFIGURATION_HOST);
            print("          configured at ...", MADNESS_CONFIGURATION_DATE);
            print("                    CXX ...", MADNESS_CONFIGURATION_CXX);
            print("               CXXFLAGS ...", MADNESS_CONFIGURATION_CXXFLAGS);
#ifdef WORLD_WATCHDOG
            print("               watchdog ...", WATCHDOG_BARK_INTERVAL, WATCHDOG_TIMEOUT);
#endif
#ifdef OPTERON_TUNE
            print("             tuning for ...", "opteron");
#elif defined(CORE_DUO_TUNE)
            print("             tuning for ...", "core duo");
#else
            print("             tuning for ...", "core2");
#endif
#ifdef BOUNDS_CHECKING
            print(" tensor bounds checking ...", "enabled");
#endif
#ifdef TENSOR_INSTANCE_COUNT
            print("  tensor instance count ...", "enabled");
#endif
            //         print(" ");
            //         IndexIterator::test();
        }

        startup(world,argc,argv);
        if (world.rank() == 0) print("Initial tensor instance count", BaseTensor::get_instance_count());
        PROFILE_BLOCK(testsuite);

        std::cout.precision(8);

        test_basic<double,1>(world);
        test_conv<double,1>(world);
        test_math<double,1>(world);
        test_diff<double,1>(world);
        test_op<double,1>(world);
        test_plot<double,1>(world);
        test_apply_push_1d<double,1>(world);

        test_io<double,1>(world);

        // stupid location for this test
        GenericConvolution1D<double,GaussianGenericFunctor<double> > gen(10,GaussianGenericFunctor<double>(100.0,100.0),0);
        GaussianConvolution1D<double> gau(10, 100.0, 100.0, 0, false);
        Tensor<double> gg = gen.rnlp(4,0);
        Tensor<double> hh = gau.rnlp(4,0);
        MADNESS_ASSERT((gg-hh).normf() < 1e-13);
        if (world.rank() == 0) print(" generic and gaussian operator kernels agree\n");

        test_qm(world);

        test_basic<double_complex,1>(world);
        test_conv<double_complex,1>(world);
        test_math<double_complex,1>(world);
        test_diff<double_complex,1>(world);
        test_op<double_complex,1>(world);
        test_plot<double_complex,1>(world);
        test_io<double_complex,1>(world);

        //TaskInterface::debug = true;
        test_basic<double,2>(world);
        test_conv<double,2>(world);
        test_math<double,2>(world);
        test_diff<double,2>(world);
        test_op<double,2>(world);
        test_plot<double,2>(world);
        test_io<double,2>(world);

        test_basic<double,3>(world);
        test_conv<double,3>(world);
        test_math<double,3>(world);
        test_diff<double,3>(world);
        test_op<double,3>(world);
        test_coulomb(world);
        test_plot<double,3>(world);
        test_io<double,3>(world);

        //test_plot<double,4>(world); // slow unless reduce npt in test_plot

        if (world.rank() == 0) print("entering final fence");
        world.gop.fence();
        if (world.rank() == 0) {
            print("done with final fence");
            print(" ");
            print("Final tensor instance count", BaseTensor::get_instance_count());
        }

        ThreadPool::end();
        print_stats(world);

    }
    catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    finalize();
    return 0;
}

