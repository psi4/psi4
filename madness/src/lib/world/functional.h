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


 $Id $
 */

#ifndef MADNESS_WORLD_FUNCTIONAL_H__INCLUDED
#define MADNESS_WORLD_FUNCTIONAL_H__INCLUDED

#include <world/typestuff.h>
#include <world/enable_if.h>
#include <functional>

namespace madness {

    namespace {
        // Parameter types a la Boost

        template <typename T>
        struct param_type {
        private:
            template <typename U, bool isa>
            struct impl {
                typedef const U & type;
            };

            template <typename U>
            struct impl<U, true> {
                typedef const T type;
            };
        public:
            typedef typename impl<T, std::is_arithmetic<T>::value>::type type;
        };

        template <typename T>
        struct param_type<T *> {
            typedef T * const type;
        };

        template <typename T>
        struct param_type<T &> {
            typedef T & type;
        };

        template <typename T, std::size_t N>
        struct param_type<T [N]> {
           typedef const T * const type;
        };

        template <typename T, std::size_t N>
        struct param_type<const T [N]> {
           typedef const T * const type;
        };
    } // namespace

// Macros for generating comma separated lists with a prefix element
// It can create a list of up to 10 elements
#define MADNESS_LISTP0( P , L ) P
#define MADNESS_LISTP1( P , L ) P , L( 0 )
#define MADNESS_LISTP2( P , L ) MADNESS_LISTP1( P , L ) , L( 1 )
#define MADNESS_LISTP3( P , L ) MADNESS_LISTP2( P , L ) , L( 2 )
#define MADNESS_LISTP4( P , L ) MADNESS_LISTP3( P , L ) , L( 3 )
#define MADNESS_LISTP5( P , L ) MADNESS_LISTP4( P , L ) , L( 4 )
#define MADNESS_LISTP6( P , L ) MADNESS_LISTP5( P , L ) , L( 5 )
#define MADNESS_LISTP7( P , L ) MADNESS_LISTP6( P , L ) , L( 6 )
#define MADNESS_LISTP8( P , L ) MADNESS_LISTP7( P , L ) , L( 7 )
#define MADNESS_LISTP9( P , L ) MADNESS_LISTP8( P , L ) , L( 8 )
#define MADNESS_LISTP10( P , L ) MADNESS_LISTP9( P , L ) , L( 9 )

#define MADNESS_LIST_PREFIX( P , L , N ) MADNESS_LISTP##N( P , L )


// Macros for generating comma separated lists
// It can create a list of up to 10 elements
#define MADNESS_LIST0( L )
#define MADNESS_LIST1( L ) L( 0 )
#define MADNESS_LIST2( L ) MADNESS_LIST1( L ) , L( 1 )
#define MADNESS_LIST3( L ) MADNESS_LIST2( L ) , L( 2 )
#define MADNESS_LIST4( L ) MADNESS_LIST3( L ) , L( 3 )
#define MADNESS_LIST5( L ) MADNESS_LIST4( L ) , L( 4 )
#define MADNESS_LIST6( L ) MADNESS_LIST5( L ) , L( 5 )
#define MADNESS_LIST7( L ) MADNESS_LIST6( L ) , L( 6 )
#define MADNESS_LIST8( L ) MADNESS_LIST7( L ) , L( 7 )
#define MADNESS_LIST9( L ) MADNESS_LIST8( L ) , L( 8 )
#define MADNESS_LIST10( L ) MADNESS_LIST9( L ) , L( 9 )

#define MADNESS_LIST( L , N ) MADNESS_LIST##N( L )

// Macros for generating numbered list elements
#define MADNESS_FN_TYPE( N )  a##N##T
#define MADNESS_FN_ARG( N )  a##N
#define MADNESS_FN_PARAM( N ) arg##N##_type MADNESS_FN_ARG( N )
#define MADNESS_FN_TPARAM( N ) typename MADNESS_FN_TYPE( N )
#define MADNESS_FN_TYPEDEF( N ) typedef typename param_type< MADNESS_FN_TYPE( N ) >::type arg##N##_type;

// Macros for generating a list of numbered typedefs
#define MADNESS_FN_TYPEDEF0
#define MADNESS_FN_TYPEDEF1 \
    MADNESS_FN_TYPEDEF( 0 )
#define MADNESS_FN_TYPEDEF2 \
    MADNESS_FN_TYPEDEF1 \
    MADNESS_FN_TYPEDEF( 1 )
#define MADNESS_FN_TYPEDEF3 \
    MADNESS_FN_TYPEDEF2 \
    MADNESS_FN_TYPEDEF( 2 )
#define MADNESS_FN_TYPEDEF4 \
    MADNESS_FN_TYPEDEF3 \
    MADNESS_FN_TYPEDEF( 3 )
#define MADNESS_FN_TYPEDEF5 \
    MADNESS_FN_TYPEDEF4 \
    MADNESS_FN_TYPEDEF( 4 )
#define MADNESS_FN_TYPEDEF6 \
    MADNESS_FN_TYPEDEF5 \
    MADNESS_FN_TYPEDEF( 5 )
#define MADNESS_FN_TYPEDEF7 \
    MADNESS_FN_TYPEDEF6 \
    MADNESS_FN_TYPEDEF( 6 )
#define MADNESS_FN_TYPEDEF8 \
    MADNESS_FN_TYPEDEF7 \
    MADNESS_FN_TYPEDEF( 7 )
#define MADNESS_FN_TYPEDEF9 \
    MADNESS_FN_TYPEDEF8 \
    MADNESS_FN_TYPEDEF( 8 )
#define MADNESS_FN_TYPEDEF10 \
    MADNESS_FN_TYPEDEF9 \
    MADNESS_FN_TYPEDEF( 9 )

#define MADNESS_FN_TYPEDEF_LIST( N ) MADNESS_FN_TYPEDEF##N

#define MADNESS_FN_CONST const
#define MADNESS_FN_NO_CONST


    template <typename fnT, typename Enable = void>
    struct PtrFn {
        static const bool value = false;
    };

// Macros for generating PtrFn objects
#define MADNESS_PTR_FN( ARITY ) \
    template < MADNESS_LIST_PREFIX( typename resT , MADNESS_FN_TPARAM , ARITY ) > \
    struct PtrFn< resT(*)( MADNESS_LIST( MADNESS_FN_TYPE , ARITY ) ),           \
            typename disable_if<std::is_void<resT> >::type  >                   \
    {                                                                           \
        static const int arity = ARITY;                                         \
        static const bool value = true;                                         \
        typedef resT result_type;                                               \
        MADNESS_FN_TYPEDEF_LIST( ARITY )                                        \
        typedef resT(*fn_type)( MADNESS_LIST( MADNESS_FN_TYPE , ARITY ) );      \
                                                                                \
    private:                                                                    \
        fn_type fn_;                                                            \
                                                                                \
    public:                                                                     \
        PtrFn() : fn_(0) { }                                                    \
        PtrFn(fn_type fn) : fn_(fn) { }                                         \
                                                                                \
        resT operator()( MADNESS_LIST( MADNESS_FN_PARAM , ARITY )  ) const {    \
            MADNESS_ASSERT(fn_);                                                \
            return fn_( MADNESS_LIST( MADNESS_FN_ARG , ARITY ) );               \
        }                                                                       \
                                                                                \
        template <typename Archive>                                             \
        void serialize(Archive& ar) { ar & archive::wrap_opaque(fn_); }         \
    };                                                                          \
                                                                                \
    template < MADNESS_LIST( MADNESS_FN_TPARAM , ARITY ) >                      \
    struct PtrFn< void(*)( MADNESS_LIST( MADNESS_FN_TYPE , ARITY ) ) >          \
    {                                                                           \
        static const int arity = ARITY;                                         \
        static const bool value = true;                                         \
        typedef void result_type;                                               \
        MADNESS_FN_TYPEDEF_LIST( ARITY )                                        \
        typedef void(*fn_type)( MADNESS_LIST( MADNESS_FN_TYPE , ARITY ) );      \
                                                                                \
    private:                                                                    \
        fn_type fn_;                                                            \
                                                                                \
    public:                                                                     \
        PtrFn() : fn_(0) { }                                                    \
        PtrFn(fn_type fn) : fn_(fn) { }                                         \
                                                                                \
        void operator()( MADNESS_LIST( MADNESS_FN_PARAM , ARITY )  ) const {    \
            MADNESS_ASSERT(fn_);                                                \
            fn_( MADNESS_LIST( MADNESS_FN_ARG , ARITY ) );                      \
        }                                                                       \
                                                                                \
        template <typename Archive>                                             \
        void serialize(Archive& ar) { ar & archive::wrap_opaque(fn_); }         \
    };

    // Generate PtrFn objects that take functions with 0 to 10 parameters.
    MADNESS_PTR_FN( 0 )
    MADNESS_PTR_FN( 1 )
    MADNESS_PTR_FN( 2 )
    MADNESS_PTR_FN( 3 )
    MADNESS_PTR_FN( 4 )
    MADNESS_PTR_FN( 5 )
    MADNESS_PTR_FN( 6 )
    MADNESS_PTR_FN( 7 )
    MADNESS_PTR_FN( 8 )
    MADNESS_PTR_FN( 9 )
    MADNESS_PTR_FN( 10 )


    template <typename fnT, typename Enable = void>
    struct MemFn {
        static const bool value = false;
    };

// Macros for generating MemFn object types
#define MADNESS_MEM_FN( ARITY , CONST ) \
    template <typename resT, MADNESS_LIST_PREFIX( typename objT , MADNESS_FN_TPARAM , ARITY ) > \
    struct MemFn< resT(objT::*)( MADNESS_LIST( MADNESS_FN_TYPE , ARITY ) ) CONST, \
            typename disable_if<std::is_void<resT> >::type >                    \
    {                                                                           \
        static const int arity = ARITY;                                         \
        static const bool value = true;                                         \
        typedef resT result_type;                                               \
        MADNESS_FN_TYPEDEF_LIST( ARITY )                                        \
        typedef CONST objT object_type;                                         \
        typedef resT(objT::*fn_type)( MADNESS_LIST( MADNESS_FN_TYPE , ARITY ) ) CONST; \
                                                                                \
    private:                                                                    \
        fn_type fn_;                                                            \
                                                                                \
    public:                                                                     \
        MemFn() : fn_(0) { }                                                    \
        MemFn(fn_type fn) : fn_(fn) { }                                         \
                                                                                \
        resT operator()(MADNESS_LIST_PREFIX(object_type& obj , MADNESS_FN_PARAM , ARITY ) ) const { \
            MADNESS_ASSERT(fn_);                                                \
            return (obj.*fn_)( MADNESS_LIST( MADNESS_FN_ARG , ARITY ) );        \
        }                                                                       \
                                                                                \
        resT operator()(MADNESS_LIST_PREFIX(object_type* obj , MADNESS_FN_PARAM , ARITY ) ) const { \
            MADNESS_ASSERT(obj);                                                \
            MADNESS_ASSERT(fn_);                                                \
            return (obj->*fn_)( MADNESS_LIST( MADNESS_FN_ARG , ARITY ) );       \
        }                                                                       \
                                                                                \
        template <typename Archive>                                             \
        void serialize(Archive& ar) { ar & archive::wrap_opaque(fn_); }         \
    };                                                                          \
                                                                                \
    template < MADNESS_LIST_PREFIX( typename objT , MADNESS_FN_TPARAM , ARITY ) > \
    struct MemFn< void(objT::*)( MADNESS_LIST( MADNESS_FN_TYPE , ARITY ) ) CONST > \
    {                                                                           \
        static const int arity = ARITY;                                         \
        static const bool value = true;                                         \
        typedef void result_type;                                               \
        MADNESS_FN_TYPEDEF_LIST( ARITY )                                        \
        typedef CONST objT object_type;                                         \
        typedef void(objT::*fn_type)( MADNESS_LIST( MADNESS_FN_TYPE , ARITY ) ) CONST; \
                                                                                \
    private:                                                                    \
        fn_type fn_;                                                            \
                                                                                \
    public:                                                                     \
        MemFn() : fn_(0) { }                                                    \
        MemFn(fn_type fn) : fn_(fn) { }                                         \
                                                                                \
        void operator()(MADNESS_LIST_PREFIX(object_type& obj , MADNESS_FN_PARAM , ARITY ) ) const { \
            MADNESS_ASSERT(fn_);                                                \
            (obj.*fn_)( MADNESS_LIST( MADNESS_FN_ARG , ARITY ) );               \
        }                                                                       \
                                                                                \
        void operator()(MADNESS_LIST_PREFIX(object_type* obj , MADNESS_FN_PARAM , ARITY ) ) const { \
            MADNESS_ASSERT(obj);                                                \
            MADNESS_ASSERT(fn_);                                                \
            (obj->*fn_)( MADNESS_LIST( MADNESS_FN_ARG , ARITY ) );              \
        }                                                                       \
                                                                                \
        template <typename Archive>                                             \
        void serialize(Archive& ar) { ar & archive::wrap_opaque(fn_); }         \
    };

    // Generate MemFn object types for member functions that take 0 to 10
    // parameters and return a type
    MADNESS_MEM_FN( 0 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 1 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 2 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 3 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 4 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 5 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 6 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 7 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 8 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 9 , MADNESS_FN_NO_CONST )
    MADNESS_MEM_FN( 10 , MADNESS_FN_NO_CONST )

    // Generate MemFn object types for const member functions that take 0 to 10
    // parameters and return a type
    MADNESS_MEM_FN( 0 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 1 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 2 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 3 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 4 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 5 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 6 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 7 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 8 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 9 , MADNESS_FN_CONST )
    MADNESS_MEM_FN( 10 , MADNESS_FN_CONST )

    /// Member function factory

    /// Creates a function object from member function pointer. The functor
    /// behaves similarly to the original function, except it adds a parameter
    /// to the beginning of the function call which is a pointer or reference
    /// to the object type. For example:
    /// \code
    /// objectT obj;
    /// obj.func(a1, a2); // Call member function of obj
    /// function(& objectT::func)(obj, a1, a2); // Equivalent call with the constructed functor
    /// \endcode
    /// \tparam fnT The member function pointer type
    /// \param fn The member function pointer for the object
    /// \return A functor object that turns a member function pointer into a functor
    template <typename fnT>
    typename enable_if_c<MemFn<fnT>::value, MemFn<fnT> >::type
    function(fnT fn) { return MemFn<fnT>(fn); }


    /// Function factory

    /// Create a functor object that behaves like a function pointer
    /// \tparam fnT The function pointer type
    /// \param fn The function pointer
    /// \return An functor that acts like the original function
    template <typename fnT>
    typename enable_if_c<PtrFn<fnT>::value, PtrFn<fnT> >::type
    function(fnT fn) { return PtrFn<fnT>(fn); }


// Remove the macros
#undef MADNESS_LISTP0
#undef MADNESS_LISTP1
#undef MADNESS_LISTP2
#undef MADNESS_LISTP3
#undef MADNESS_LISTP4
#undef MADNESS_LISTP5
#undef MADNESS_LISTP6
#undef MADNESS_LISTP7
#undef MADNESS_LISTP8
#undef MADNESS_LISTP9
#undef MADNESS_LISTP10
#undef MADNESS_LIST_PREFIX
#undef MADNESS_LIST0
#undef MADNESS_LIST1
#undef MADNESS_LIST2
#undef MADNESS_LIST3
#undef MADNESS_LIST4
#undef MADNESS_LIST5
#undef MADNESS_LIST6
#undef MADNESS_LIST7
#undef MADNESS_LIST8
#undef MADNESS_LIST9
#undef MADNESS_LIST10
#undef MADNESS_LIST
#undef MADNESS_FN_TYPE
#undef MADNESS_FN_ARG
#undef MADNESS_FN_PARAM
#undef MADNESS_FN_TPARAM
#undef MADNESS_FN_TYPEDEF
#undef MADNESS_FN_TYPEDEF0
#undef MADNESS_FN_TYPEDEF1
#undef MADNESS_FN_TYPEDEF2
#undef MADNESS_FN_TYPEDEF3
#undef MADNESS_FN_TYPEDEF4
#undef MADNESS_FN_TYPEDEF5
#undef MADNESS_FN_TYPEDEF6
#undef MADNESS_FN_TYPEDEF7
#undef MADNESS_FN_TYPEDEF8
#undef MADNESS_FN_TYPEDEF9
#undef MADNESS_FN_TYPEDEF10
#undef MADNESS_FN_TYPEDEF_LIST
#undef MADNESS_FN_CONST
#undef MADNESS_FN_NO_CONST
#undef MADNESS_PTR_FN
#undef MADNESS_MEM_FN

    // Provide serialization for standard library functional types
    namespace archive {

        template <class, class>
        struct ArchiveSerializeImpl;

        template <typename Archive, typename T>
        struct ArchiveSerializeImpl<Archive, std::plus<T> > {
            static inline void serialize(const Archive&, std::plus<T>&) { }
        };

        template <typename Archive, typename T>
        struct ArchiveSerializeImpl<Archive, std::minus<T> > {
            static inline void serialize(const Archive&, std::minus<T>&) { }
        };

        template <typename Archive, typename T>
        struct ArchiveSerializeImpl<Archive, std::multiplies<T> > {
            static inline void serialize(const Archive&, std::multiplies<T>&) { }
        };


        template <typename Archive, typename T>
        struct ArchiveSerializeImpl<Archive, std::divides<T> > {
            static inline void serialize(const Archive&, std::divides<T>&) { }
        };


        template <typename Archive, typename T>
        struct ArchiveSerializeImpl<Archive, std::negate<T> > {
            static inline void serialize(const Archive&, std::negate<T>&) { }
        };

        template <typename Archive, typename T>
        struct ArchiveSerializeImpl<Archive, std::modulus<T> > {
            static inline void serialize(const Archive&, std::modulus<T>&) { }
        };

    } // namespace archive

} // namespace madness

#endif // MADNESS_WORLD_FUNCTIONAL_H__INCLUDED
