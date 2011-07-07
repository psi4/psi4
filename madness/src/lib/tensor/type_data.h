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


  $Id: type_data.h 1861 2010-04-13 15:03:52Z justus.c79 $
*/


#ifndef MADNESS_TENSOR_TYPE_DATA_H__INCLUDED
#define MADNESS_TENSOR_TYPE_DATA_H__INCLUDED

/// \file type_data.h
/// \brief Defines and implements TensorTypeData, a type traits class.

namespace madness {

    /// Traits class to specify support of numeric types.

    /// This traits class is used to specify which numeric types are
    /// supported by the tensor library and also their unique integer id.
    /// Unsupported types will default to an entry supported=false.
    ///
    /// To add a new type, append it to the definitions below using
    /// TYPEINFO and tensor_type_names, incrementing the id, and set
    /// TENSOR_MAX_TYPE_ID accordingly.  You might also have to specialize
    /// some of the methods in tensor.cc, and will have to add additional
    /// instantiations at the end of tensor.cc and tensoriter.cc.
    template <class T>
    class TensorTypeData {
    public:
        enum {id = -1};
        enum {supported = false};
        enum {iscomplex = false};
        enum {memcopyok = false};
        typedef T type;
        //typedef T scalar_type;
    };

    /// This provides the reverse mapping from integer id to type name
    template <int id>
    class TensorTypeFromId {
    public:
        typedef long type;
    };

    // id=unique and sequential identifier for each type
    //
    // supported=true for all supported scalar numeric types
    //
    // iscomplex=true if a complex type
    //
    // memcopyok=true if memcpy can be used to copy an array
    // of the type ... it should be true for all native types
    // and probably false for all types that require a
    // special constructor or assignment operator.
    //
    // type=the actual type
    //
    // scalar_type = is the type of abs, normf, absmin, absmax, real, imag.
    // Not all of these functions are defined for all types.
    // Unfortunately, in the current version, some misuses will only be
    // detected at run time.


#define TYPEINFO(num, T, iscmplx, mcpyok, realT,floatrealT) \
template<> class TensorTypeData<T> {\
public: \
  enum {id = num}; \
  enum {supported = true}; \
  enum {iscomplex = iscmplx}; \
  enum {memcopyok = mcpyok}; \
  typedef T type; \
  typedef realT scalar_type; \
  typedef floatrealT float_scalar_type; \
}; \
template<> class TensorTypeFromId<num> {\
public: \
  typedef T type; \
}

    TYPEINFO(0,int,false,true,int,double);
    TYPEINFO(1,long,false,true,long,double);
    TYPEINFO(2,float,false,true,float,float);
    TYPEINFO(3,double,false,true,double,double);
    TYPEINFO(4,float_complex,true,true,float,float);
    TYPEINFO(5,double_complex,true,true,double,double);
#define TENSOR_MAX_TYPE_ID 5
#undef TYPEINFO

#ifdef TENSOR_CC
    const char *tensor_type_names[TENSOR_MAX_TYPE_ID+1] = {
                "int","long","float","double","float_complex","double_complex"
            };
#else
    extern const char *tensor_type_names[];
#endif

    /// The template IsSupported is used to constrain instantiation of
    /// templates to the supported scalar types.  It is only implemented if
    /// the type is supported, in which case it evaluates to the return
    /// type.
    ///
    /// E.g., to restrict operator+ to supported types T that can be added
    /// to type A, with a return type of A.
    ///
    /// template <typename T> typename IsSupported < TensorTypeData<T>, A >::type
    /// operator+(A const &v, T const &w) {
    ///    ...
    ///    return something of type A
    /// };

    template <typename TypeData, typename, bool = TypeData::supported>
    struct IsSupported;

    template <typename TypeData, typename ReturnType>
    struct IsSupported <TypeData, ReturnType, true> {
        typedef ReturnType type;
    };

    // This macro embodies the above for a single parameter template.
    // Look in tensor.h for example of how to use in multi-parameter
    // templates.
    //
    // *** !!! *** Unfortunately, this macro was confusing Doxygen
    // so we have removed all usages.
    //#define ISSUPPORTED(T,RETURNTYPE)
    //template <typename T> typename IsSupported < TensorTypeData<T>, RETURNTYPE >::type

    /// TensorResultType<L,R>::type is the type of (L op R) where op is nominally multiplication
    template <typename leftT, typename rightT>
    struct TensorResultType {};

#define SPEC(L,R,T) \
    template <> struct TensorResultType<L,R> {typedef T type;}; \
    template <> struct TensorResultType<R,L> {typedef T type;}
#define DPEC(L,R,T) \
    template <> struct TensorResultType<L,L> {typedef T type;}

    DPEC(int,int,int);
    SPEC(int,long,long);
    SPEC(int,float,float);
    SPEC(int,double,double);
    SPEC(int,float_complex,float_complex);
    SPEC(int,double_complex,double_complex);
    DPEC(long,long,long);
    SPEC(long,float,float);
    SPEC(long,double,double);
    SPEC(long,float_complex,float_complex);
    SPEC(long,double_complex,double_complex);
    DPEC(float,float,float);
    SPEC(float,double,double);
    SPEC(float,float_complex,float_complex);
    SPEC(float,double_complex,double_complex);
    DPEC(double,double,double);
    SPEC(double,float_complex,float_complex);
    SPEC(double,double_complex,double_complex);
    DPEC(float_complex,float_complex,float_complex);
    SPEC(float_complex,double_complex,double_complex);
    DPEC(double_complex,double_complex,double_complex);

    /// This macro simplifies access to TensorResultType
#define TENSOR_RESULT_TYPE(L,R) typename TensorResultType<L,R>::type


#undef DPEC
#undef SPEC
}
#endif // MADNESS_TENSOR_TYPE_DATA_H__INCLUDED
