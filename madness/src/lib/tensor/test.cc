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


  $Id: test.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/

/// \file tensor/test.cc
/// \brief New test code for Tensor class using Google unit test

#include <tensor/tensor.h>
#include <world/print.h>

#ifdef HAVE_GOOGLE_TEST

// The test code deliberately uses only the dumb ITERATOR macros
// in order to test the optimized iterators used by the implementation.

using madness::_;
using madness::___;
using madness::_reverse;

#include <iostream>
#include <gtest/gtest.h>

namespace {

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

    template <typename T>
    class TensorTest : public ::testing::Test {
    public:
        TensorTest() {}

        virtual ~TensorTest() {}

        virtual void SetUp() {}

        virtual void TearDown() {}

        static T* move_test_ptr;

        static madness::Tensor<T> f() { // For testing move semantics
            madness::Tensor<T> result(3);
            move_test_ptr = result.ptr();
            return result;
        }
    };

    template <typename T> T* TensorTest<T>::move_test_ptr;

    typedef ::testing::Types<int, long, float, double, float_complex, double_complex> TensorTestTypes;
    TYPED_TEST_CASE(TensorTest, TensorTestTypes);

    TYPED_TEST(TensorTest, Basic) {
        for (int ndim=0; ndim<=TENSOR_MAXDIM; ++ndim) {
            try {
                std::vector<long> dim(ndim);
                for (int pass=0; pass<10; ++pass) {
                    long nelem = 1;
                    for (int i=0; i<ndim; ++i) {
                        dim[i] = (madness::RandomValue<long>()&0x5) + 1;
                        nelem *= dim[i];
                    }

                    // Use ASSERT not EXPECT otherwise we get a huge volume of output
                    // from failing tests.  Also, subsequent tests usually rely on
                    // succcess of previous ones.

                    madness::Tensor<TypeParam> empty; // Verify default constructor
                    ASSERT_EQ(empty.size(),0);
                    ASSERT_EQ(empty.ptr(),(TypeParam*)0);
                    ASSERT_EQ(empty.ndim(),-1);
                    ASSERT_EQ(empty.id(),madness::TensorTypeData<TypeParam>::id);

                    madness::Tensor<TypeParam> d(dim);
                    ASSERT_EQ(d.id(),madness::TensorTypeData<TypeParam>::id);
                    ASSERT_EQ(d.ndim(),ndim);
                    ASSERT_EQ(d.size(),nelem);
                    ASSERT_NE(d.ptr(),(TypeParam*)0);
                    for (int i=0; i<ndim; ++i) {
                        ASSERT_EQ(d.dim(i),dim[i]);
                    }
                    // Should be initialized to zero
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(0)));

                    d.fillindex();   // Check fillindex and indexing
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index)));

                    d.fill(TypeParam(33)); // Check fill
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(33)));

                    d = TypeParam(21);     // Check fill
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(21)));

                    madness::Tensor<TypeParam> q(d); // Copy constructor
                    ASSERT_EQ(q.id(),madness::TensorTypeData<TypeParam>::id);
                    ASSERT_EQ(q.ptr(),d.ptr()); // Was it shallow?
                    ASSERT_EQ(q.id(),madness::TensorTypeData<TypeParam>::id);
                    ASSERT_EQ(q.ndim(),ndim);
                    ASSERT_EQ(q.size(),nelem);
                    ITERATOR(d, ASSERT_EQ(q(IND),d(IND)));

                    q.clear();  // Check that clear is working
                    ASSERT_EQ(q.size(),0);
                    ASSERT_EQ(q.ptr(),(TypeParam*)0);
                    ASSERT_EQ(q.ndim(),-1);
                    ASSERT_EQ(q.id(),madness::TensorTypeData<TypeParam>::id);

                    q = d;    // Assignment
                    ASSERT_EQ(q.ptr(),d.ptr()); // Was it shallow?
                    ITERATOR(d, ASSERT_EQ(q(IND),d(IND)));

                    q = copy(d); // Deep copy
                    ASSERT_NE(q.ptr(),d.ptr()); // Was it deep?
                    ITERATOR(d, ASSERT_EQ(q(IND),d(IND)));

                    madness::Tensor<TypeParam> r = TensorTest<TypeParam>::f(); // Should invoke move semantics
                    ASSERT_EQ(r.ptr(),TensorTest<TypeParam>::move_test_ptr);

                    r.clear();  // Does clear work?
                    ASSERT_EQ(r.size(),0);
                    ASSERT_EQ(r.ptr(),(TypeParam*)0);
                    ASSERT_EQ(r.ndim(),-1);
                    ASSERT_EQ(r.id(),madness::TensorTypeData<TypeParam>::id);

                    r = TensorTest<TypeParam>::f(); // Should invoke move semantics
                    ASSERT_EQ(r.ptr(),TensorTest<TypeParam>::move_test_ptr);

                    q.fillindex();

                    d.fillindex();

                    q *= TypeParam(2); // Check inplace scaling
                    d.scale(TypeParam(3));
                    ITERATOR(q, ASSERT_EQ(q(IND),TypeParam(_index*2)));
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index*3)));

                    d += TypeParam(2); // Check inplace scalar addition
                    ITERATOR(q, ASSERT_EQ(q(IND),TypeParam(_index*2)));
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index*3+2)));

                    d -= q;     // Check inplace tensor subtraction
                    ITERATOR(q, ASSERT_EQ(q(IND),TypeParam(_index*2)));
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index+2)));

                    d += q;     // Check inplace tensor addition
                    ITERATOR(q, ASSERT_EQ(q(IND),TypeParam(_index*2)));
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index*3+2)));
                }
            }
            catch (const madness::TensorException& e) {
                if (ndim != 0) std::cout << e;
                EXPECT_EQ(ndim,0);
            }
            catch(...) {
                std::cout << "Caught unknown exception" << std::endl;
                EXPECT_EQ(1,0);
            }
        }
    }

    TYPED_TEST(TensorTest, Slicing1d) {
        madness::Tensor<TypeParam> t(5L);
        t(_) = 3;
        ITERATOR(t,EXPECT_EQ(t(IND),TypeParam(3)));
    }


    TYPED_TEST(TensorTest, SliceAndDice) {
        madness::Tensor<TypeParam> a(5,6,7), b;

        a.fillindex();
        ITERATOR3(a,ASSERT_EQ(a(_,_,_)(IND3),a(IND3)));

        b = a(___);
        ITERATOR3(b,ASSERT_EQ(b(IND3),a(IND3)));
        b.fillrandom();      // verify it is indeed a view of the same thing
        ITERATOR3(b,ASSERT_EQ(b(IND3),a(IND3)));

        b = a(_reverse,_,_);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(4-_i,_j,_k)));
        b.fillrandom();      // verify it is indeed a view of the same thing
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(4-_i,_j,_k)));

        b = a(_,_reverse,_);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k),a(_i,5-_j,_k)));
        b.fillrandom();      // verify it is indeed a view of the same thing
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k),a(_i,5-_j,_k)));

        b = a(_,_,_reverse);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k),a(_i,_j,6-_k)));
        b.fillrandom();      // verify it is indeed a view of the same thing
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k),a(_i,_j,6-_k)));

        // Middle of each dimension
        b = a(madness::Slice(1,-2,1),madness::Slice(3,-2,1),madness::Slice(2,-2,1));
        ASSERT_EQ(b.dim(0),(a.dim(0)-2));
        ASSERT_EQ(b.dim(1),(a.dim(1)-4));
        ASSERT_EQ(b.dim(2),(a.dim(2)-3));
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k),a(1+_i,_j+3,2+_k)));

        // subpatch assignment, chaining, defaulting args on slice ...
        // 5 6 7 view property
        a.fill(0);
        b = a(madness::Slice(1,-1),madness::Slice(1,-2),madness::Slice(2,-2)) = TypeParam(1);
        ASSERT_EQ(long(std::abs(a.sum())+0.1), (a.dim(0)-1)*(a.dim(1)-2)*(a.dim(2)-3));
        ITERATOR3(b,ASSERT_EQ(b(IND3), TypeParam(1)));
        ITERATOR3(b,ASSERT_EQ(a(_i+1,_j+1,_k+2), TypeParam(1)));

        // elimination of dimensions
        a.fill(0);
        b = a(2,_,_) = TypeParam(1);
        ASSERT_EQ(b.ndim(),2);
        ASSERT_EQ(b.dim(0),a.dim(1));
        ASSERT_EQ(b.dim(1),a.dim(2));
        ITERATOR2(b,ASSERT_EQ(b(IND2), a(2,_i,_j)));

        // interaction with reordering, etc
        b = a.swapdim(0,2)(_,-1,_);//  a(i,j,k)->a(k,j,i)->b(k,i) with j=last
        ASSERT_TRUE(a.iscontiguous());
        ASSERT_FALSE(b.iscontiguous());
        ASSERT_EQ(b.ndim(),2);
        ASSERT_EQ(b.dim(0),a.dim(2));
        ASSERT_EQ(b.dim(1),a.dim(0));
        ITERATOR2(b,ASSERT_EQ(b(IND2),a(_j,5,_i)));

        // from my slice to your slice, ...
        b = madness::Tensor<TypeParam>(5,6,7);
        ASSERT_EQ(b.normf(),TypeParam(0));
        b.fillindex();

        a = TypeParam(0);
        a(___) = b(___);

        ASSERT_NE(a.ptr(),b.ptr()); // Ensure it was deep
        ITERATOR3(a,ASSERT_EQ(a(IND3),b(IND3)));

        madness::Tensor<long> g(5,6,7);
        a.fill(0);
        g.fillindex();
        a(___) = g(___);		// test mixed type slice assignment
        ITERATOR3(a,ASSERT_EQ(a(IND3),TypeParam(g(IND3))));

        // test filling a sub patch
        b = madness::Tensor<TypeParam>(5,6,7);		// Fresh ones to be clear about what we have
        a = madness::Tensor<TypeParam>(5,6,7);
        b.fillrandom();
        std::vector<madness::Slice> patch = madness::vector_factory(madness::Slice(1,-2),madness::Slice(1,-2),madness::Slice(1,-2));

        a.fill(0);
        madness::Tensor<TypeParam> c = a(patch) = b(patch);
        ITERATOR3(c,ASSERT_EQ(c(_i,_j,_k), a(_i+1,_j+1,_k+1)));
        ITERATOR3(c,ASSERT_EQ(c(_i,_j,_k), b(_i+1,_j+1,_k+1)));

        // check did not go out of bounds
        EXPECT_NEAR(c.normf(), a.normf(), 2e-7);
        EXPECT_NEAR(b(patch).normf(), a.normf(), 2e-7);

        // patch assignment interaction with reordering of source and target

        a = madness::Tensor<TypeParam>(7,6,5);
        b = madness::Tensor<TypeParam>(5,7,6);

        a = a.swapdim(0,2);
        a(patch).fillindex();
        b = b.swapdim(1,2);
        b(patch) = a(patch);		// Assign patch of reordered a to patch of b.

        long n = a(patch).size();
        ASSERT_EQ(long(std::abs(a.sum())), n*(n-1)/2);
        ASSERT_EQ(long(std::abs(b.sum())),n*(n-1)/2);
        ITERATOR3(c,ASSERT_EQ(a(_i,_j,_k), b(_i,_j,_k)));
    }


    TYPED_TEST(TensorTest, Transpose) {
        for (long n=1; n<=21; ++n) {
            madness::Tensor<TypeParam> a(n,n+3);
            ASSERT_TRUE(a.iscontiguous());
            ASSERT_EQ(a.size(),n*(n+3));
            a.fillrandom();
            ITERATOR(a, ASSERT_NE(a(IND),TypeParam(0)));

            madness::Tensor<TypeParam> aT = transpose(a);
            ASSERT_TRUE(aT.iscontiguous());
            ASSERT_EQ(aT.size(),n*(n+3));
            ASSERT_NE(aT.ptr(),a.ptr()); // Was it deep?
            ITERATOR2(aT, ASSERT_EQ(aT(_i,_j),a(_j,_i)));

            a.fillrandom();
            aT.clear();
            aT = transpose(a);
            ASSERT_EQ(aT.size(),n*(n+3));
            ASSERT_TRUE(a.iscontiguous());
            ASSERT_TRUE(aT.iscontiguous());
            ASSERT_NE(aT.ptr(),a.ptr()); // Was it deep?
            ITERATOR2(aT, ASSERT_EQ(aT(_i,_j),a(_j,_i)));

            aT = a.swapdim(0,1);
            ASSERT_EQ(aT.size(),n*(n+3));
            ASSERT_TRUE(a.iscontiguous());
            ASSERT_FALSE(aT.iscontiguous());
            ASSERT_EQ(aT.ptr(),a.ptr()); // Was it shallow?
            ITERATOR2(aT, ASSERT_EQ(aT(_i,_j),a(_j,_i)));
        }
    }


    TYPED_TEST(TensorTest, Reshaping) {
        madness::Tensor<TypeParam>  a(4,6,10);
        a.fillrandom();

        EXPECT_THROW(a.reshape(4), madness::TensorException);

        madness::Tensor<TypeParam> b = a.reshape(3,2,2,5,2,2);
        long i=0;
        ITERATOR6(b,ASSERT_EQ(b(IND6), a.ptr()[i++]));

        b = a.reshape(30,8);
        i = 0;
        ITERATOR2(b,ASSERT_EQ(b(IND2), a.ptr()[i++]));

        b = a.flat();
        i = 0;
        ITERATOR1(b,ASSERT_EQ(b(IND1), a.ptr()[i++]));

        b = a.reshape(240);
        i = 0;
        ITERATOR1(b,ASSERT_EQ(b(IND1), a.ptr()[i++]));

        b = a.swapdim(1,2);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_i,_k,_j)));

        EXPECT_TRUE(a.iscontiguous());
        EXPECT_FALSE(b.iscontiguous());

        b = a.swapdim(0,2);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_k,_j,_i)));

        EXPECT_THROW(b.reshape(240), madness::TensorException);
        EXPECT_THROW(b.flat(), madness::TensorException);

        b = a.cycledim(1,0,-1);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_j,_k,_i)));

        b = a.cycledim(-1,0,-1);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_k,_i,_j)));

        b = b.cycledim(-1,0,-1);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_j,_k,_i)));

        b = b.cycledim(-1,0,-1);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_i,_j,_k)));

        b = a.cycledim(1,0,1);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_j,_i,_k)));

        b = b.cycledim(1,0,1);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_i,_j,_k)));

        b = a.cycledim(2,0,-1);
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_k,_i,_j)));

        b = a.mapdim(madness::vector_factory<long>(1,0,2));
        ITERATOR3(b,ASSERT_EQ(b(_i,_j,_k), a(_j,_i,_k)));
    }

//     TYPED_TEST(TensorTest, Container) {
//         typedef madness::ConcurrentHashMap< int, Tensor<TypeParam> > containerT;
//         static const int N = 100;
//         containerT c;
//         Tensor<TypeParam> a[N];
//         for (int i=0; i<N; ++i) {

//         }
//     }

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the tensor test code\n";
    return 0;
}

#endif
