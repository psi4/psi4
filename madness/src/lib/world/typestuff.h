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


  $Id: typestuff.h 2387 2011-06-23 15:53:29Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_TYPESTUFF_H__INCLUDED
#define MADNESS_WORLD_TYPESTUFF_H__INCLUDED

/// \file typestuff.h
/// \brief Grossly simplified Boost-like type traits and templates

// Aims to be a compatible subset of Boost ... wish we could rely on
// Boost being there but between bjam (theirs), bloat (theirs), brain
// deficiences (mine, not theirs), bleeding edge massively parallel
// systems, and compiler deficiences, we have to be independent for the
// time being.

#include <cstddef>
#include <stdint.h>
#include <madness_config.h>
#include <world/enable_if.h>

// Select header that contains type_traits
#if defined(MADNESS_USE_TYPE_TRAITS)
#include <type_traits>
#elif defined(MADNESS_USE_TR1_TYPE_TRAITS)
#include <tr1/type_traits>
#elif defined(MADNESS_USE_BOOST_TR1_TYPE_TRAITS_HPP)
#include <boost/tr1/type_traits.hpp>
#else
#error No acceptable type_traits include directive was found.
#endif // TYPE_TRAITS

#if defined(MADNESS_HAS_STD_TR1_TYPE_TRAITS) && !defined(MADNESS_HAS_STD_TYPE_TRAITS)
#define MADNESS_HAS_STD_TYPE_TRAITS 1

// Insert the tr1 type traits into the std namespace.
namespace std {

    using ::std::tr1::integral_constant;
    using ::std::tr1::true_type;
    using ::std::tr1::false_type;
    using ::std::tr1::is_void;
    using ::std::tr1::is_integral;
    using ::std::tr1::is_floating_point;
    using ::std::tr1::is_array;
    using ::std::tr1::is_pointer;
    using ::std::tr1::is_reference;
    using ::std::tr1::is_member_object_pointer;
    using ::std::tr1::is_member_function_pointer;
    using ::std::tr1::is_enum;
    using ::std::tr1::is_union;
    using ::std::tr1::is_class;
    using ::std::tr1::is_function;
    using ::std::tr1::is_arithmetic;
    using ::std::tr1::is_fundamental;
    using ::std::tr1::is_object;
    using ::std::tr1::is_scalar;
    using ::std::tr1::is_compound;
    using ::std::tr1::is_member_pointer;
    using ::std::tr1::is_const;
    using ::std::tr1::is_volatile;
    using ::std::tr1::is_pod;
    using ::std::tr1::is_empty;
    using ::std::tr1::is_polymorphic;
    using ::std::tr1::is_abstract;
    using ::std::tr1::has_trivial_constructor;
    using ::std::tr1::has_trivial_copy;
    using ::std::tr1::has_trivial_assign;
    using ::std::tr1::has_trivial_destructor;
    using ::std::tr1::has_nothrow_constructor;
    using ::std::tr1::has_nothrow_copy;
    using ::std::tr1::has_nothrow_assign;
    using ::std::tr1::has_virtual_destructor;
    using ::std::tr1::is_signed;
    using ::std::tr1::is_unsigned;
    using ::std::tr1::alignment_of;
    using ::std::tr1::rank;
    using ::std::tr1::extent;
    using ::std::tr1::is_same;
    using ::std::tr1::is_base_of;
    using ::std::tr1::is_convertible;
    using ::std::tr1::remove_const;
    using ::std::tr1::remove_volatile;
    using ::std::tr1::remove_cv;
    using ::std::tr1::add_const;
    using ::std::tr1::add_volatile;
    using ::std::tr1::add_cv;
    using ::std::tr1::remove_reference;
    using ::std::tr1::add_reference;
    using ::std::tr1::remove_extent;
    using ::std::tr1::remove_all_extents;
    using ::std::tr1::remove_pointer;
    using ::std::tr1::add_pointer;
    using ::std::tr1::aligned_storage;

} // namespace std
#endif


namespace madness {

    template <bool Cond, typename T1, typename T2>
    struct if_c {
        typedef T1 type;
    };

    template <typename T1, typename T2>
    struct if_c<false, T1, T2> {
        typedef T2 type;
    };

    template <typename Cond, typename T1, typename T2>
    struct if_ : public if_c<Cond::value, T1, T2> {};

    /// type_or_c<bool A, bool B>::value will be true if (A || B)
    template <bool A, bool B>
    class type_or_c : public std::true_type { };

    /// type_or_c<bool A, bool B>::value will be true if (A || B)
    template<> class type_or_c<false,false> : public std::false_type { };

    /// type_or<CondA,CondB>::value  will be true if (CondA::value || CondB::value)
    template <class CondA, class CondB>
    class type_or : public type_or_c<CondA::value, CondB::value> { };


    /// type_and_c<bool A, bool B>::value will be true if (A && B)
    template <bool A, bool B>
    class type_and_c : public std::false_type { };

    /// type_and_c<bool A, bool B>::value will be true if (A && B)
    template<> class type_and_c<true, true> : public std::true_type { };

    /// type_and<CondA,CondB>::value  will be true if (CondA::value && CondB::value)
    template <class CondA, class CondB>
    class type_and: public type_and_c<CondA::value, CondB::value> { };

    /// is_eq<A,B> returns true if A and B are the same integers
    template <int A, int B>
    struct is_eq  : public std::false_type { };

    /// is_eq<A,B> returns true if A and B are the same integers
    template <int A>
    struct is_eq<A,A> : public std::true_type { };

    /// True if A is derived from B and is not B
    template <class A, class B>
    struct is_derived_from {
        typedef char yes;
        typedef int no;
        static no f(...);
        static yes f(B*);
        static const bool value = (sizeof(f((A*)0)) == sizeof(yes));
    };

    /// True if A is derived from B and is not B
    template <class A>
    struct is_derived_from<A,A> {
        static const bool value = false;
    };


    /// Function traits in the spirt of boost function traits
    template <typename functionT>
    struct function_traits {
        static const bool value = false;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT>
    struct function_traits<returnT(*)()> {
        static const bool value = true;
        static const int arity = 0;
        typedef returnT result_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T>
    struct function_traits<returnT(*)(arg1T)> {
        static const bool value = true;
        static const int arity = 1;
        typedef returnT result_type;
        typedef arg1T arg1_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T>
    struct function_traits<returnT(*)(arg1T,arg2T)> {
        static const bool value = true;
        static const int arity = 2;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T>
    struct function_traits<returnT(*)(arg1T,arg2T,arg3T)> {
        static const bool value = true;
        static const int arity = 3;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
    struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T)> {
        static const bool value = true;
        static const int arity = 4;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
    struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T)> {
        static const bool value = true;
        static const int arity = 5;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
    };


    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
    struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)> {
        static const bool value = true;
        static const int arity = 6;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
    };


    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
    struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)> {
        static const bool value = true;
        static const int arity = 7;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
    struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T)> {
        static const bool value = true;
        static const int arity = 8;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
    struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T)> {
        static const bool value = true;
        static const int arity = 9;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
        typedef arg9T arg9_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename memfuncT>
    struct memfunc_traits {
        static const bool value = false;
    };

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT>
    struct memfunc_traits<returnT(objT::*)()> {
        static const bool value = true;
        static const int arity = 0;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
    };

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T>
    struct memfunc_traits<returnT(objT::*)(arg1T)> {
        static const bool value = true;
        static const int arity = 1;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
    };

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T)> {
        static const bool value = true;
        static const int arity = 2;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
    };

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T)> {
        static const bool value = true;
        static const int arity = 3;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
    };

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T)> {
        static const bool value = true;
        static const int arity = 4;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T)> {
        static const bool value = true;
        static const int arity = 5;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)> {
        static const bool value = true;
        static const int arity = 6;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)> {
        static const bool value = true;
        static const int arity = 7;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T)> {
        static const bool value = true;
        static const int arity = 8;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T)> {
        static const bool value = true;
        static const int arity = 9;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
        typedef arg9T arg9_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT>
    struct memfunc_traits<returnT(objT::*)() const> {
        static const bool value = true;
        static const int arity = 0;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
    };

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T>
    struct memfunc_traits<returnT(objT::*)(arg1T) const> {
        static const bool value = true;
        static const int arity = 1;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
    };

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T) const> {
        static const bool value = true;
        static const int arity = 2;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
    };

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T) const> {
        static const bool value = true;
        static const int arity = 3;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T) const> {
        static const bool value = true;
        static const int arity = 4;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
    };



    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T) const> {
        static const bool value = true;
        static const int arity = 5;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T) const> {
        static const bool value = true;
        static const int arity = 6;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T) const> {
        static const bool value = true;
        static const int arity = 7;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
    };



    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T) const> {
        static const bool value = true;
        static const int arity = 8;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
    struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T) const> {
        static const bool value = true;
        static const int arity = 9;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
        typedef arg9T arg9_type;
    };


    /// is_function_ptr<T>::value is true if T is a function pointer
//    template <typename T>
//    struct is_function_ptr {
//        static const bool value = function_traits<T>::value;
//    };

    template <typename fnT, typename Enabler = void>
    struct result_of {
        typedef typename fnT::result_type type;
    };

    template <typename fnT>
    struct result_of<fnT, typename enable_if_c<function_traits<fnT>::value>::type> {
        typedef typename function_traits<fnT>::result_type type;
    };

    template <typename fnT>
    struct result_of<fnT, typename enable_if_c<memfunc_traits<fnT>::value>::type> {
        typedef typename memfunc_traits<fnT>::result_type type;
    };

    template <typename> class Future;
    template <typename> struct remove_future;

    // Remove Future, const, volatile, and reference qualifiers from the type
    template <typename T>
    struct remove_fcvr {
        typedef typename remove_future<typename std::remove_cv<
                typename std::remove_reference<T>::type>::type>::type type;
    };

    /// This defines stuff that is serialiable by default rules ... basically anything contiguous

    /// Fundamental types, member function pointers, and function pointers
    /// are handled by default.
    template <typename T>
    struct is_serializable : public std::integral_constant<bool,
            std::is_fundamental<T>::value ||
            std::is_member_function_pointer<T>::value ||
            (std::is_function<typename std::remove_pointer<T>::type>::value &&
            ::std::is_pointer<T>::value)>
    {};


    /// Simple binder for member functions with no arguments
    template <class T, typename resultT>
    class BindNullaryMemFun {
    private:
        T* t;
        resultT(T::*op)();
    public:
        BindNullaryMemFun(T* t, resultT(T::*op)()) : t(t), op(op) {};
        resultT operator()() {
            return (t->*op)();
        };
    };

    /// Specialization of BindNullaryMemFun for void return
    template <class T>
    class BindNullaryMemFun<T,void> {
    private:
        T* t;
        void (T::*op)();
    public:
        BindNullaryMemFun(T* t, void (T::*op)()) : t(t), op(op) {};
        void operator()() {
            (t->*op)();
        };
    };


    /// Simple binder for const member functions with no arguments
    template <class T, typename resultT>
    class BindNullaryConstMemFun {
    private:
        const T* t;
        resultT(T::*op)() const;
    public:
        BindNullaryConstMemFun(const T* t, resultT(T::*op)() const) : t(t), op(op) {};
        resultT operator()() const {
            return (t->*op)();
        };
    };

    /// Specialization of BindNullaryConstMemFun for void return
    template <class T>
    class BindNullaryConstMemFun<T,void> {
    private:
        const T* t;
        void (T::*op)() const;
    public:
        BindNullaryConstMemFun(const T* t, void (T::*op)() const) : t(t), op(op) {};
        void operator()() const {
            (t->*op)();
        };
    };

    /// Factory function for BindNullaryMemFun
    template <class T, typename resultT>
    inline BindNullaryMemFun<T,resultT> bind_nullary_mem_fun(T* t, resultT(T::*op)()) {
        return BindNullaryMemFun<T,resultT>(t,op);
    }

    /// Factory function for BindNullaryConstMemFun
    template <class T, typename resultT>
    inline BindNullaryConstMemFun<T,resultT> bind_nullary_mem_fun(const T* t, resultT(T::*op)() const) {
        return BindNullaryConstMemFun<T,resultT>(t,op);
    }

    /// A type you can return when you want to return void ... use "return None"
    struct Void {};

    /// None, a la Python
    static const Void None = Void();

    /// Wrapper so that can return something even if returning void
    template <typename T>
    struct ReturnWrapper {
        typedef T type;
    };

    /// Wrapper so that can return something even if returning void
    template <>
    struct ReturnWrapper<void> {
        typedef Void type;
    };


    /// Used to provide rvalue references to support move semantics
    template <typename T> class Reference {
        T* p;
    public:
        Reference(T* x) : p(x) {}
        T& operator*() const {return *p;}
        T* operator->() const {return p;}
    };


    /* Macros to make some of this stuff more readable */

    /**
       \def REMREF(TYPE)
       \brief Macro to make std::remove_reference<T> easier to use

       \def REMCONST(TYPE)
       \brief Macro to make std::remove_const<T> easier to use

       \def MEMFUN_RETURNT(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_OBJT(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARITY(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARG1T(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARG2T(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARG3T(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARG4T(TYPE)
       \brief Macro to make member function type traits easier to use

       \def FUNCTION_RETURNT(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARITY(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARG1T(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARG2T(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARG3T(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARG4T(TYPE)
       \brief Macro to make function type traits easier to use

       \def IS_SAME(TYPEA, TYPEB)
       \brief Macro to make is_same<T> template easier to use

       \def RETURN_WRAPPERT(TYPE)
       \brief Returns TYPE unless TYPE is void when returns ReturnWrapper<void>

    */


#define REMREF(TYPE)    typename std::remove_reference< TYPE >::type
#define REMCONST(TYPE)  typename std::remove_const< TYPE >::type
#define REMCONSTX(TYPE) std::remove_const< TYPE >::type
#define RETURN_WRAPPERT(TYPE) typename ReturnWrapper< TYPE >::type

#define MEMFUN_RETURNT(MEMFUN) typename memfunc_traits< MEMFUN >::result_type
#define MEMFUN_CONSTNESS(MEMFUN) memfunc_traits< MEMFUN >::constness
#define MEMFUN_OBJT(MEMFUN)    typename memfunc_traits< MEMFUN >::obj_type
#define MEMFUN_ARITY(MEMFUN)   memfunc_traits< MEMFUN >::arity
#define MEMFUN_ARG1T(MEMFUN)   typename memfunc_traits< MEMFUN >::arg1_type
#define MEMFUN_ARG2T(MEMFUN)   typename memfunc_traits< MEMFUN >::arg2_type
#define MEMFUN_ARG3T(MEMFUN)   typename memfunc_traits< MEMFUN >::arg3_type
#define MEMFUN_ARG4T(MEMFUN)   typename memfunc_traits< MEMFUN >::arg4_type

#define FUNCTION_RETURNT(FUNCTION) typename function_traits< FUNCTION >::result_type
#define FUNCTION_ARITY(FUNCTION)   function_traits< FUNCTION >::arity
#define FUNCTION_ARG1T(FUNCTION)   typename function_traits< FUNCTION >::arg1_type
#define FUNCTION_ARG2T(FUNCTION)   typename function_traits< FUNCTION >::arg2_type
#define FUNCTION_ARG3T(FUNCTION)   typename function_traits< FUNCTION >::arg3_type
#define FUNCTION_ARG4T(FUNCTION)   typename function_traits< FUNCTION >::arg4_type

#define RESULT_OF(FUNCTION) typename madness::result_of< FUNCTION >::type

#define IS_SAME(A, B) std::is_same< A, B >
#define IS_EQ(A, B) is_eq< A, B >

} // end of namespace madness
#endif // MADNESS_WORLD_TYPESTUFF_H__INCLUDED
