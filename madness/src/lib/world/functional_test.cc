#include <madness_config.h>

#ifdef MADNESS_HAS_GOOGLE_TEST

#include <world/world.h>
#include <world/functional.h>
#include <gtest/gtest.h>

using namespace madness;


int f0() { return 0; }
int f1(int) { return 1; }
int f2(int,int) { return 2; }
int f3(int,int,int) { return 3; }
int f4(int,int,int,int) { return 4; }
int f5(int,int,int,int,int) { return 5; }
int f6(int,int,int,int,int,int) { return 6; }
int f7(int,int,int,int,int,int,int) { return 7; }
int f8(int,int,int,int,int,int,int,int) { return 8; }
int f9(int,int,int,int,int,int,int,int,int) { return 9; }
int f10(int,int,int,int,int,int,int,int,int,int) { return 10; }

void vf0() { }
void vf1(int) { }
void vf2(int,int) { }
void vf3(int,int,int) { }
void vf4(int,int,int,int) { }
void vf5(int,int,int,int,int) { }
void vf6(int,int,int,int,int,int) { }
void vf7(int,int,int,int,int,int,int) { }
void vf8(int,int,int,int,int,int,int,int) { }
void vf9(int,int,int,int,int,int,int,int,int) { }
void vf10(int,int,int,int,int,int,int,int,int,int) { }

struct Tester {
    int f0() { return 0; }
    int f1(int) { return 1; }
    int f2(int,int) { return 2; }
    int f3(int,int,int) { return 3; }
    int f4(int,int,int,int) { return 4; }
    int f5(int,int,int,int,int) { return 5; }
    int f6(int,int,int,int,int,int) { return 6; }
    int f7(int,int,int,int,int,int,int) { return 7; }
    int f8(int,int,int,int,int,int,int,int) { return 8; }
    int f9(int,int,int,int,int,int,int,int,int) { return 9; }
    int f10(int,int,int,int,int,int,int,int,int,int) { return 10; }


    int f0c() const { return 0; }
    int f1c(int) const { return 1; }
    int f2c(int,int) const { return 2; }
    int f3c(int,int,int) const { return 3; }
    int f4c(int,int,int,int) const { return 4; }
    int f5c(int,int,int,int,int) const { return 5; }
    int f6c(int,int,int,int,int,int) const { return 6; }
    int f7c(int,int,int,int,int,int,int) const { return 7; }
    int f8c(int,int,int,int,int,int,int,int) const { return 8; }
    int f9c(int,int,int,int,int,int,int,int,int) const { return 9; }
    int f10c(int,int,int,int,int,int,int,int,int,int) const { return 10; }

    void vf0() { }
    void vf1(int) { }
    void vf2(int,int) { }
    void vf3(int,int,int) { }
    void vf4(int,int,int,int) { }
    void vf5(int,int,int,int,int) { }
    void vf6(int,int,int,int,int,int) { }
    void vf7(int,int,int,int,int,int,int) { }
    void vf8(int,int,int,int,int,int,int,int) { }
    void vf9(int,int,int,int,int,int,int,int,int) { }
    void vf10(int,int,int,int,int,int,int,int,int,int) { }

    void vf0c() const { }
    void vf1c(int) const { }
    void vf2c(int,int) const { }
    void vf3c(int,int,int) const { }
    void vf4c(int,int,int,int) const { }
    void vf5c(int,int,int,int,int) const { }
    void vf6c(int,int,int,int,int,int) const { }
    void vf7c(int,int,int,int,int,int,int) const { }
    void vf8c(int,int,int,int,int,int,int,int) const { }
    void vf9c(int,int,int,int,int,int,int,int,int) const { }
    void vf10c(int,int,int,int,int,int,int,int,int,int) const { }
};

namespace {
    class FunctionalTest : public ::testing::Test {
    public:
    };

    TEST_F(FunctionalTest, PtrFnConstruct) {

        // Test constructor for function pointers with a return type
        EXPECT_NO_THROW(PtrFn<int(*)()> f(&f0));
        EXPECT_NO_THROW(PtrFn<int(*)(int)> f(&f1));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int)> f(&f2));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int,int)> f(&f3));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int,int,int)> f(&f4));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int,int,int,int)> f(&f5));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int,int,int,int,int)> f(&f6));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int,int,int,int,int,int)> f(&f7));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int,int,int,int,int,int,int)> f(&f8));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int,int,int,int,int,int,int,int)> f(&f9));
        EXPECT_NO_THROW(PtrFn<int(*)(int,int,int,int,int,int,int,int,int,int)> f(&f10));

        // Test constructor for functions pointers with a void return type
        EXPECT_NO_THROW(PtrFn<void(*)()> f(&vf0));
        EXPECT_NO_THROW(PtrFn<void(*)(int)> f(&vf1));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int)> f(&vf2));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int,int)> f(&vf3));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int,int,int)> f(&vf4));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int,int,int,int)> f(&vf5));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int,int,int,int,int)> f(&vf6));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int,int,int,int,int,int)> f(&vf7));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int,int,int,int,int,int,int)> f(&vf8));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int,int,int,int,int,int,int,int)> f(&vf9));
        EXPECT_NO_THROW(PtrFn<void(*)(int,int,int,int,int,int,int,int,int,int)> f(&vf10));
    }

    TEST_F(FunctionalTest, PtrFnFactory) {

        // Test factory functions for function pointers with a return type
        EXPECT_NO_THROW(function(&f0));
        EXPECT_NO_THROW(function(&f1));
        EXPECT_NO_THROW(function(&f2));
        EXPECT_NO_THROW(function(&f3));
        EXPECT_NO_THROW(function(&f4));
        EXPECT_NO_THROW(function(&f5));
        EXPECT_NO_THROW(function(&f6));
        EXPECT_NO_THROW(function(&f7));
        EXPECT_NO_THROW(function(&f8));
        EXPECT_NO_THROW(function(&f9));
        EXPECT_NO_THROW(function(&f10));

        // Test factory functions for function pointers with a void return type
        EXPECT_NO_THROW(function(&vf0));
        EXPECT_NO_THROW(function(&vf1));
        EXPECT_NO_THROW(function(&vf2));
        EXPECT_NO_THROW(function(&vf3));
        EXPECT_NO_THROW(function(&vf4));
        EXPECT_NO_THROW(function(&vf5));
        EXPECT_NO_THROW(function(&vf6));
        EXPECT_NO_THROW(function(&vf7));
        EXPECT_NO_THROW(function(&vf8));
        EXPECT_NO_THROW(function(&vf9));
        EXPECT_NO_THROW(function(&vf10));
    }

    TEST_F(FunctionalTest, CallPtrFn) {

        // Test call for function pointers with a return type
        EXPECT_EQ(0, function(&f0)());
        EXPECT_EQ(1, function(&f1)(1));
        EXPECT_EQ(2, function(&f2)(1,1));
        EXPECT_EQ(3, function(&f3)(1,1,1));
        EXPECT_EQ(4, function(&f4)(1,1,1,1));
        EXPECT_EQ(5, function(&f5)(1,1,1,1,1));
        EXPECT_EQ(6, function(&f6)(1,1,1,1,1,1));
        EXPECT_EQ(7, function(&f7)(1,1,1,1,1,1,1));
        EXPECT_EQ(8, function(&f8)(1,1,1,1,1,1,1,1));
        EXPECT_EQ(9, function(&f9)(1,1,1,1,1,1,1,1,1));
        EXPECT_EQ(10, function(&f10)(1,1,1,1,1,1,1,1,1,1));

        // Test call for function pointers with a void return type
        EXPECT_NO_THROW(function(&vf0)());
        EXPECT_NO_THROW(function(&vf1)(1));
        EXPECT_NO_THROW(function(&vf2)(1,1));
        EXPECT_NO_THROW(function(&vf3)(1,1,1));
        EXPECT_NO_THROW(function(&vf4)(1,1,1,1));
        EXPECT_NO_THROW(function(&vf5)(1,1,1,1,1));
        EXPECT_NO_THROW(function(&vf6)(1,1,1,1,1,1));
        EXPECT_NO_THROW(function(&vf7)(1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(&vf8)(1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(&vf9)(1,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(&vf10)(1,1,1,1,1,1,1,1,1,1));
    }

    TEST_F(FunctionalTest, MemFnConstruct) {

        // Test member functions with a return type

        // Test constructor with just a member function pointer
        EXPECT_NO_THROW(MemFn<int(Tester::*)()> f(& Tester::f0));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int)> f(& Tester::f1));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int)> f(& Tester::f2));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int)> f(& Tester::f3));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int)> f(& Tester::f4));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int)> f(& Tester::f5));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int)> f(& Tester::f6));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int,int)> f(& Tester::f7));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int,int,int)> f(& Tester::f8));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int,int,int,int)> f(& Tester::f9));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int,int,int,int,int)> f(& Tester::f10));


        // Test member functions with void return type

        // Test constructor with just a member function pointer
        EXPECT_NO_THROW(MemFn<void(Tester::*)()> f(& Tester::vf0));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int)> f(& Tester::vf1));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int)> f(& Tester::vf2));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int)> f(& Tester::vf3));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int)> f(& Tester::vf4));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int)> f(& Tester::vf5));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int)> f(& Tester::vf6));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int,int)> f(& Tester::vf7));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int,int,int)> f(& Tester::vf8));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int,int,int,int)> f(& Tester::vf9));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int,int,int,int,int)> f(& Tester::vf10));


        // Test const member functions with a return type

        // Test constructor with just a member function pointer
        EXPECT_NO_THROW(MemFn<int(Tester::*)() const> f(& Tester::f0c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int) const> f(& Tester::f1c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int) const> f(& Tester::f2c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int) const> f(& Tester::f3c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int) const> f(& Tester::f4c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int) const> f(& Tester::f5c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int) const> f(& Tester::f6c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int,int) const> f(& Tester::f7c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int,int,int) const> f(& Tester::f8c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int,int,int,int) const> f(& Tester::f9c));
        EXPECT_NO_THROW(MemFn<int(Tester::*)(int,int,int,int,int,int,int,int,int,int) const> f(& Tester::f10c));


        // Test const member functions with void return type

        // Test constructor with just a member function pointer
        EXPECT_NO_THROW(MemFn<void(Tester::*)() const> f(& Tester::vf0c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int) const> f(& Tester::vf1c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int) const> f(& Tester::vf2c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int) const> f(& Tester::vf3c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int) const> f(& Tester::vf4c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int) const> f(& Tester::vf5c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int) const> f(& Tester::vf6c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int,int) const> f(& Tester::vf7c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int,int,int) const> f(& Tester::vf8c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int,int,int,int) const> f(& Tester::vf9c));
        EXPECT_NO_THROW(MemFn<void(Tester::*)(int,int,int,int,int,int,int,int,int,int) const> f(& Tester::vf10c));
    }

    TEST_F(FunctionalTest, MemFnCall) {

        Tester t;

        // Test calls for member functions with a return type

        // Test calls for objects with an object reference
        EXPECT_EQ(0, function(& Tester::f0)(t));
        EXPECT_EQ(1, function(& Tester::f1)(t,1));
        EXPECT_EQ(2, function(& Tester::f2)(t,1,1));
        EXPECT_EQ(3, function(& Tester::f3)(t,1,1,1));
        EXPECT_EQ(4, function(& Tester::f4)(t,1,1,1,1));
        EXPECT_EQ(5, function(& Tester::f5)(t,1,1,1,1,1));
        EXPECT_EQ(6, function(& Tester::f6)(t,1,1,1,1,1,1));
        EXPECT_EQ(7, function(& Tester::f7)(t,1,1,1,1,1,1,1));
        EXPECT_EQ(8, function(& Tester::f8)(t,1,1,1,1,1,1,1,1));
        EXPECT_EQ(9, function(& Tester::f9)(t,1,1,1,1,1,1,1,1,1));
        EXPECT_EQ(10, function(& Tester::f10)(t,1,1,1,1,1,1,1,1,1,1));

        // Test calls for objects with an object pointer
        EXPECT_EQ(0, function(& Tester::f0)(&t));
        EXPECT_EQ(1, function(& Tester::f1)(&t,1));
        EXPECT_EQ(2, function(& Tester::f2)(&t,1,1));
        EXPECT_EQ(3, function(& Tester::f3)(&t,1,1,1));
        EXPECT_EQ(4, function(& Tester::f4)(&t,1,1,1,1));
        EXPECT_EQ(5, function(& Tester::f5)(&t,1,1,1,1,1));
        EXPECT_EQ(6, function(& Tester::f6)(&t,1,1,1,1,1,1));
        EXPECT_EQ(7, function(& Tester::f7)(&t,1,1,1,1,1,1,1));
        EXPECT_EQ(8, function(& Tester::f8)(&t,1,1,1,1,1,1,1,1));
        EXPECT_EQ(9, function(& Tester::f9)(&t,1,1,1,1,1,1,1,1,1));
        EXPECT_EQ(10, function(& Tester::f10)(&t,1,1,1,1,1,1,1,1,1,1));
    }

    TEST_F(FunctionalTest, MemFnCallConst) {

        Tester t;
        // Test calls for const member functions with a return type

        // Test calls for objects with an object reference
        EXPECT_EQ(0, function(& Tester::f0c)(t));
        EXPECT_EQ(1, function(& Tester::f1c)(t,1));
        EXPECT_EQ(2, function(& Tester::f2c)(t,1,1));
        EXPECT_EQ(3, function(& Tester::f3c)(t,1,1,1));
        EXPECT_EQ(4, function(& Tester::f4c)(t,1,1,1,1));
        EXPECT_EQ(5, function(& Tester::f5c)(t,1,1,1,1,1));
        EXPECT_EQ(6, function(& Tester::f6c)(t,1,1,1,1,1,1));
        EXPECT_EQ(7, function(& Tester::f7c)(t,1,1,1,1,1,1,1));
        EXPECT_EQ(8, function(& Tester::f8c)(t,1,1,1,1,1,1,1,1));
        EXPECT_EQ(9, function(& Tester::f9c)(t,1,1,1,1,1,1,1,1,1));
        EXPECT_EQ(10, function(& Tester::f10c)(t,1,1,1,1,1,1,1,1,1,1));

        // Test calls for objects with an object pointer
        EXPECT_EQ(0, function(& Tester::f0c)(&t));
        EXPECT_EQ(1, function(& Tester::f1c)(&t,1));
        EXPECT_EQ(2, function(& Tester::f2c)(&t,1,1));
        EXPECT_EQ(3, function(& Tester::f3c)(&t,1,1,1));
        EXPECT_EQ(4, function(& Tester::f4c)(&t,1,1,1,1));
        EXPECT_EQ(5, function(& Tester::f5c)(&t,1,1,1,1,1));
        EXPECT_EQ(6, function(& Tester::f6c)(&t,1,1,1,1,1,1));
        EXPECT_EQ(7, function(& Tester::f7c)(&t,1,1,1,1,1,1,1));
        EXPECT_EQ(8, function(& Tester::f8c)(&t,1,1,1,1,1,1,1,1));
        EXPECT_EQ(9, function(& Tester::f9c)(&t,1,1,1,1,1,1,1,1,1));
        EXPECT_EQ(10, function(& Tester::f10c)(&t,1,1,1,1,1,1,1,1,1,1));
    }

    TEST_F(FunctionalTest, MemFnCallVoid) {

        Tester t;

        // Test calls for member functions with a void type

        // Test calls for objects with an object reference
        EXPECT_NO_THROW(function(& Tester::vf0)(t));
        EXPECT_NO_THROW(function(& Tester::vf1)(t,1));
        EXPECT_NO_THROW(function(& Tester::vf2)(t,1,1));
        EXPECT_NO_THROW(function(& Tester::vf3)(t,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf4)(t,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf5)(t,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf6)(t,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf7)(t,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf8)(t,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf9)(t,1,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf10)(t,1,1,1,1,1,1,1,1,1,1));

        // Test calls for objects with an object pointer
        EXPECT_NO_THROW(function(& Tester::vf0)(&t));
        EXPECT_NO_THROW(function(& Tester::vf1)(&t,1));
        EXPECT_NO_THROW(function(& Tester::vf2)(&t,1,1));
        EXPECT_NO_THROW(function(& Tester::vf3)(&t,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf4)(&t,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf5)(&t,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf6)(&t,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf7)(&t,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf8)(&t,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf9)(&t,1,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf10)(&t,1,1,1,1,1,1,1,1,1,1));
    }

    TEST_F(FunctionalTest, MemFnCallVoidConst) {

        Tester t;

        // Test calls for const member functions with a void type

        // Test calls for objects with an object reference
        EXPECT_NO_THROW(function(& Tester::vf0c)(t));
        EXPECT_NO_THROW(function(& Tester::vf1c)(t,1));
        EXPECT_NO_THROW(function(& Tester::vf2c)(t,1,1));
        EXPECT_NO_THROW(function(& Tester::vf3c)(t,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf4c)(t,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf5c)(t,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf6c)(t,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf7c)(t,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf8c)(t,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf9c)(t,1,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf10c)(t,1,1,1,1,1,1,1,1,1,1));

        // Test calls for objects with an object pointer
        EXPECT_NO_THROW(function(& Tester::vf0c)(&t));
        EXPECT_NO_THROW(function(& Tester::vf1c)(&t,1));
        EXPECT_NO_THROW(function(& Tester::vf2c)(&t,1,1));
        EXPECT_NO_THROW(function(& Tester::vf3c)(&t,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf4c)(&t,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf5c)(&t,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf6c)(&t,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf7c)(&t,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf8c)(&t,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf9c)(&t,1,1,1,1,1,1,1,1,1));
        EXPECT_NO_THROW(function(& Tester::vf10c)(&t,1,1,1,1,1,1,1,1,1,1));
    }
} // namespace

int main(int argc, char **argv) {
    madness::initialize(argc,argv);
    int status;
    {
        madness::World world(MPI::COMM_WORLD);

        if (world.rank()) madness::redirectio(world);
        world.args(argc,argv);
        world.gop.fence();

        ::testing::InitGoogleTest(&argc, argv);
        status = RUN_ALL_TESTS();

        world.gop.fence();
    }
    madness::finalize();
    return status;
}


#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the world test code\n";
    return 0;
}

#endif
