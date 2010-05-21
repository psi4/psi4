//  (C) Copyright Gennadiy Rozental 2001-2010.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision: 62016 $
//
//  Description : contains definition for all test tools in test toolbox
// ***************************************************************************

#ifndef BOOST_TEST_TEST_TOOLS_HPP_012705GER
#define BOOST_TEST_TEST_TOOLS_HPP_012705GER

// Boost.Test
#include <boost/test/predicate_result.hpp>
#ifndef BOOST_TEST_PROD
#include <boost/test/unit_test_log.hpp>
#define BOOST_TEST_TOOLS_STRINGIZE( arg ) BOOST_TEST_STRINGIZE( arg )
#else
#define BOOST_TEST_PASSPOINT()
#define BOOST_TEST_TOOLS_STRINGIZE( arg ) BOOST_STRINGIZE( arg )
#endif
#include <boost/test/floating_point_comparison.hpp>

#include <boost/test/detail/config.hpp>
#include <boost/test/detail/global_typedef.hpp>
#include <boost/test/detail/workaround.hpp>

#include <boost/test/utils/wrap_stringstream.hpp>
#include <boost/test/utils/basic_cstring/io.hpp>
#include <boost/test/utils/lazy_ostream.hpp>

// Boost
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>

#include <boost/limits.hpp>

#include <boost/type_traits/is_array.hpp>
#include <boost/type_traits/is_function.hpp>
#include <boost/type_traits/is_abstract.hpp>

#include <boost/mpl/or.hpp>

// STL
#include <cstddef>          // for std::size_t
#include <iosfwd>
#include <ios>              // for std::boolalpha
#include <climits>          // for CHAR_BIT

#ifdef BOOST_MSVC
# pragma warning(disable: 4127) // conditional expression is constant
#endif

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

// ************************************************************************** //
// **************                    TOOL BOX                  ************** //
// ************************************************************************** //

// In macros below following argument abbreviations are used:
// P - predicate
// M - message
// S - statement
// E - exception
// L - left argument
// R - right argument
// TL - tool level
// CT - check type
// ARGS - arguments list (as PP sequence)

// frwd_type:
// 0 - args exists and need to be forwarded; call check_frwd
// 1 - args exists, but do not need to be forwarded; call check_impl directly
// 2 - no arguments; call check_impl directly

#ifndef BOOST_TEST_PROD 

#define BOOST_TEST_TOOL_PASS_PRED0( P, ARGS ) P
#define BOOST_TEST_TOOL_PASS_PRED1( P, ARGS ) P BOOST_PP_SEQ_TO_TUPLE(ARGS)
#define BOOST_TEST_TOOL_PASS_PRED2( P, ARGS ) P

#define BOOST_TEST_TOOL_PASS_ARG( r, _, arg ) , arg, BOOST_STRINGIZE( arg )
#define BOOST_TEST_TOOL_PASS_ARG_DSCR( r, _, arg ) , BOOST_STRINGIZE( arg )

#define BOOST_TEST_TOOL_PASS_ARGS0( ARGS ) \
    BOOST_PP_SEQ_FOR_EACH( BOOST_TEST_TOOL_PASS_ARG, _, ARGS )
#define BOOST_TEST_TOOL_PASS_ARGS1( ARGS ) \
    , BOOST_PP_SEQ_SIZE(ARGS) BOOST_PP_SEQ_FOR_EACH( BOOST_TEST_TOOL_PASS_ARG_DSCR, _, ARGS )
#define BOOST_TEST_TOOL_PASS_ARGS2( ARGS ) \
    , 0

#define BOOST_TEST_TOOL_IMPL(frwd_type, P, assertion_descr,TL,CT,ARGS)  \
do {                                                                    \
    BOOST_TEST_PASSPOINT();                                             \
    ::boost::test_tools::tt_detail::                                    \
    BOOST_PP_IF( frwd_type, check_impl, check_frwd ) (                  \
        BOOST_JOIN( BOOST_TEST_TOOL_PASS_PRED, frwd_type )( P, ARGS ),  \
        ::boost::unit_test::lazy_ostream::instance() << assertion_descr,\
        BOOST_TEST_L(__FILE__),                                         \
        static_cast<std::size_t>(__LINE__),                             \
        ::boost::test_tools::tt_detail::TL,                             \
        ::boost::test_tools::tt_detail::CT                              \
        BOOST_JOIN( BOOST_TEST_TOOL_PASS_ARGS, frwd_type )( ARGS ) );   \
} while( ::boost::test_tools::dummy_cond )                              \
/**/
#endif

//____________________________________________________________________________//

#define BOOST_WARN( P )                     BOOST_TEST_TOOL_IMPL( 2, (P), BOOST_TEST_TOOLS_STRINGIZE( P ), WARN, CHECK_PRED, _ )
#define BOOST_CHECK( P )                    BOOST_TEST_TOOL_IMPL( 2, (P), BOOST_TEST_TOOLS_STRINGIZE( P ), CHECK, CHECK_PRED, _ )
#define BOOST_REQUIRE( P )                  BOOST_TEST_TOOL_IMPL( 2, (P), BOOST_TEST_TOOLS_STRINGIZE( P ), REQUIRE, CHECK_PRED, _ )

//____________________________________________________________________________//

#define BOOST_WARN_MESSAGE( P, M )          BOOST_TEST_TOOL_IMPL( 2, (P), M, WARN, CHECK_MSG, _ )
#define BOOST_CHECK_MESSAGE( P, M )         BOOST_TEST_TOOL_IMPL( 2, (P), M, CHECK, CHECK_MSG, _ )
#define BOOST_REQUIRE_MESSAGE( P, M )       BOOST_TEST_TOOL_IMPL( 2, (P), M, REQUIRE, CHECK_MSG, _ )

//____________________________________________________________________________//

#define BOOST_ERROR( M )                    BOOST_CHECK_MESSAGE( false, M )
#define BOOST_FAIL( M )                     BOOST_REQUIRE_MESSAGE( false, M )

//____________________________________________________________________________//

#define BOOST_CHECK_THROW_IMPL( S, E, P, prefix, TL )                                                   \
    try {                                                                                               \
        BOOST_TEST_PASSPOINT();                                                                         \
        S;                                                                                              \
        BOOST_TEST_TOOL_IMPL( 2,false,"exception " BOOST_STRINGIZE(E) " is expected",TL,CHECK_MSG,_);  \
    } catch( E const& ex ) {                                                                            \
        ::boost::unit_test::ut_detail::ignore_unused_variable_warning( ex );                            \
        BOOST_TEST_TOOL_IMPL( 2, P, prefix BOOST_STRINGIZE( E ) " is caught", TL, CHECK_MSG, _ );      \
    }                                                                                                   \
/**/

//____________________________________________________________________________//

#define BOOST_WARN_THROW( S, E )            BOOST_CHECK_THROW_IMPL( S, E, true, "exception ", WARN )
#define BOOST_CHECK_THROW( S, E )           BOOST_CHECK_THROW_IMPL( S, E, true, "exception ", CHECK )
#define BOOST_REQUIRE_THROW( S, E )         BOOST_CHECK_THROW_IMPL( S, E, true, "exception ", REQUIRE )

//____________________________________________________________________________//

#define BOOST_WARN_EXCEPTION( S, E, P )     BOOST_CHECK_THROW_IMPL( S, E, P( ex ), "incorrect exception ", WARN )
#define BOOST_CHECK_EXCEPTION( S, E, P )    BOOST_CHECK_THROW_IMPL( S, E, P( ex ), "incorrect exception ", CHECK )
#define BOOST_REQUIRE_EXCEPTION( S, E, P )  BOOST_CHECK_THROW_IMPL( S, E, P( ex ), "incorrect exception ", REQUIRE )

//____________________________________________________________________________//

#define BOOST_CHECK_NO_THROW_IMPL( S, TL )                                                              \
    try {                                                                                               \
        S;                                                                                              \
        BOOST_TEST_TOOL_IMPL( 2,true,"no exceptions thrown by " BOOST_STRINGIZE( S ),TL,CHECK_MSG,_);  \
    } catch( ... ) {                                                                                    \
        BOOST_TEST_TOOL_IMPL( 2,false,"exception thrown by " BOOST_STRINGIZE( S ), TL, CHECK_MSG, _ ); \
    }                                                                                                   \
/**/

#define BOOST_WARN_NO_THROW( S )            BOOST_CHECK_NO_THROW_IMPL( S, WARN )
#define BOOST_CHECK_NO_THROW( S )           BOOST_CHECK_NO_THROW_IMPL( S, CHECK )
#define BOOST_REQUIRE_NO_THROW( S )         BOOST_CHECK_NO_THROW_IMPL( S, REQUIRE )

//____________________________________________________________________________//

#define BOOST_WARN_EQUAL( L, R )            BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::equal_impl_frwd(), "", WARN, CHECK_EQUAL, (L)(R) )
#define BOOST_CHECK_EQUAL( L, R )           BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::equal_impl_frwd(), "", CHECK, CHECK_EQUAL, (L)(R) )
#define BOOST_REQUIRE_EQUAL( L, R )         BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::equal_impl_frwd(), "", REQUIRE, CHECK_EQUAL, (L)(R) )

//____________________________________________________________________________//

#define BOOST_WARN_NE( L, R )               BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::ne_impl(), "", WARN, CHECK_NE, (L)(R) )
#define BOOST_CHECK_NE( L, R )              BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::ne_impl(), "", CHECK, CHECK_NE, (L)(R) )
#define BOOST_REQUIRE_NE( L, R )            BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::ne_impl(), "", REQUIRE, CHECK_NE, (L)(R) )

//____________________________________________________________________________//

#define BOOST_WARN_LT( L, R )               BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::lt_impl(), "", WARN, CHECK_LT, (L)(R) )
#define BOOST_CHECK_LT( L, R )              BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::lt_impl(), "", CHECK, CHECK_LT, (L)(R) )
#define BOOST_REQUIRE_LT( L, R )            BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::lt_impl(), "", REQUIRE, CHECK_LT, (L)(R) )

//____________________________________________________________________________//

#define BOOST_WARN_LE( L, R )               BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::le_impl(), "", WARN, CHECK_LE, (L)(R) )
#define BOOST_CHECK_LE( L, R )              BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::le_impl(), "", CHECK, CHECK_LE, (L)(R) )
#define BOOST_REQUIRE_LE( L, R )            BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::le_impl(), "", REQUIRE, CHECK_LE, (L)(R) )

//____________________________________________________________________________//

#define BOOST_WARN_GT( L, R )               BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::gt_impl(), "", WARN, CHECK_GT, (L)(R) )
#define BOOST_CHECK_GT( L, R )              BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::gt_impl(), "", CHECK, CHECK_GT, (L)(R) )
#define BOOST_REQUIRE_GT( L, R )            BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::gt_impl(), "", REQUIRE, CHECK_GT, (L)(R) )

//____________________________________________________________________________//

#define BOOST_WARN_GE( L, R )               BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::ge_impl(), "", WARN, CHECK_GE, (L)(R) )
#define BOOST_CHECK_GE( L, R )              BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::ge_impl(), "", CHECK, CHECK_GE, (L)(R) )
#define BOOST_REQUIRE_GE( L, R )            BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::tt_detail::ge_impl(), "", REQUIRE, CHECK_GE, (L)(R) )

//____________________________________________________________________________//

#define BOOST_WARN_CLOSE( L, R, T )         BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_close, "", WARN, CHECK_CLOSE, (L)(R)(::boost::math::fpc::percent_tolerance(T)) )
#define BOOST_CHECK_CLOSE( L, R, T )        BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_close, "", CHECK, CHECK_CLOSE, (L)(R)(::boost::math::fpc::percent_tolerance(T)) )
#define BOOST_REQUIRE_CLOSE( L, R, T )      BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_close, "", REQUIRE, CHECK_CLOSE, (L)(R)(::boost::math::fpc::percent_tolerance(T)) )

//____________________________________________________________________________//

#define BOOST_WARN_CLOSE_FRACTION(L, R, T)  BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_close, "", WARN, CHECK_CLOSE_FRACTION, (L)(R)(T) )
#define BOOST_CHECK_CLOSE_FRACTION(L, R, T) BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_close, "", CHECK, CHECK_CLOSE_FRACTION, (L)(R)(T) )
#define BOOST_REQUIRE_CLOSE_FRACTION(L,R,T) BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_close, "", REQUIRE, CHECK_CLOSE_FRACTION, (L)(R)(T) )

//____________________________________________________________________________//

#define BOOST_WARN_SMALL( FPV, T )          BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_small, "", WARN, CHECK_SMALL, (FPV)(T) )
#define BOOST_CHECK_SMALL( FPV, T )         BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_small, "", CHECK, CHECK_SMALL, (FPV)(T) )
#define BOOST_REQUIRE_SMALL( FPV, T )       BOOST_TEST_TOOL_IMPL( 0, \
    ::boost::test_tools::check_is_small, "", REQUIRE, CHECK_SMALL, (FPV)(T) )

//____________________________________________________________________________//

#define BOOST_WARN_PREDICATE( P, ARGS )     BOOST_TEST_TOOL_IMPL( 0, \
    P, BOOST_TEST_TOOLS_STRINGIZE( P ), WARN, CHECK_PRED_WITH_ARGS, ARGS )
#define BOOST_CHECK_PREDICATE( P, ARGS )    BOOST_TEST_TOOL_IMPL( 0, \
    P, BOOST_TEST_TOOLS_STRINGIZE( P ), CHECK, CHECK_PRED_WITH_ARGS, ARGS )
#define BOOST_REQUIRE_PREDICATE( P, ARGS )  BOOST_TEST_TOOL_IMPL( 0, \
    P, BOOST_TEST_TOOLS_STRINGIZE( P ), REQUIRE, CHECK_PRED_WITH_ARGS, ARGS )

//____________________________________________________________________________//

#define BOOST_WARN_EQUAL_COLLECTIONS( L_begin, L_end, R_begin, R_end )          \
    BOOST_TEST_TOOL_IMPL( 1, ::boost::test_tools::tt_detail::equal_coll_impl(), \
        "", WARN, CHECK_EQUAL_COLL, (L_begin)(L_end)(R_begin)(R_end) )          \
/**/
#define BOOST_CHECK_EQUAL_COLLECTIONS( L_begin, L_end, R_begin, R_end )         \
    BOOST_TEST_TOOL_IMPL( 1, ::boost::test_tools::tt_detail::equal_coll_impl(), \
        "", CHECK, CHECK_EQUAL_COLL, (L_begin)(L_end)(R_begin)(R_end) )         \
/**/
#define BOOST_REQUIRE_EQUAL_COLLECTIONS( L_begin, L_end, R_begin, R_end )       \
    BOOST_TEST_TOOL_IMPL( 1, ::boost::test_tools::tt_detail::equal_coll_impl(), \
        "", REQUIRE, CHECK_EQUAL_COLL, (L_begin)(L_end)(R_begin)(R_end) )       \
/**/

//____________________________________________________________________________//

#define BOOST_WARN_BITWISE_EQUAL( L, R )    BOOST_TEST_TOOL_IMPL( 1, \
    ::boost::test_tools::tt_detail::bitwise_equal_impl(), "", WARN, CHECK_BITWISE_EQUAL, (L)(R) )
#define BOOST_CHECK_BITWISE_EQUAL( L, R )   BOOST_TEST_TOOL_IMPL( 1, \
    ::boost::test_tools::tt_detail::bitwise_equal_impl(), "", CHECK, CHECK_BITWISE_EQUAL, (L)(R) )
#define BOOST_REQUIRE_BITWISE_EQUAL( L, R ) BOOST_TEST_TOOL_IMPL( 1, \
    ::boost::test_tools::tt_detail::bitwise_equal_impl(), "", REQUIRE, CHECK_BITWISE_EQUAL, (L)(R) )

//____________________________________________________________________________//

#define BOOST_IS_DEFINED( symb ) ::boost::test_tools::tt_detail::is_defined_impl( #symb, BOOST_STRINGIZE(= symb) )

//____________________________________________________________________________//

// ***************************** //
// deprecated interface

#define BOOST_BITWISE_EQUAL( L, R )         BOOST_CHECK_BITWISE_EQUAL( L, R )
#define BOOST_MESSAGE( M )                  BOOST_TEST_MESSAGE( M )
#define BOOST_CHECKPOINT( M )               BOOST_TEST_CHECKPOINT( M )

namespace boost {

namespace test_tools {

typedef unit_test::const_string      const_string;

namespace { bool dummy_cond = false; }

// ************************************************************************** //
// **************                print_log_value               ************** //
// ************************************************************************** //

template<typename T>
struct print_log_value {
    void    operator()( std::ostream& ostr, T const& t )
    {
        // avoid warning: 'boost::test_tools::<unnamed>::dummy_cond' defined but not used 
        if (::boost::test_tools::dummy_cond) {}

        typedef typename mpl::or_<is_array<T>,is_function<T>,is_abstract<T> >::type cant_use_nl;

        set_precision( ostr, cant_use_nl() );

        ostr << t; // by default print the value
    }

    void set_precision( std::ostream& ostr, mpl::false_ )
    {
        if( std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::radix == 2 )
            ostr.precision( 2 + std::numeric_limits<T>::digits * 301/1000 ); 
    }

    void set_precision( std::ostream&, mpl::true_ ) {}
};

//____________________________________________________________________________//

#define BOOST_TEST_DONT_PRINT_LOG_VALUE( the_type )         \
namespace boost { namespace test_tools {                    \
template<>                                                  \
struct print_log_value<the_type > {                         \
    void operator()( std::ostream&, the_type const& ) {}    \
};                                                          \
}}                                                          \
/**/

//____________________________________________________________________________//

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
template<typename T, std::size_t N >
struct print_log_value< T[N] > {
    void    operator()( std::ostream& ostr, T const* t )
    {
        ostr << t;
    }
};
#endif

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<bool> {
    void    operator()( std::ostream& ostr, bool t )
    {
         ostr << std::boolalpha << t;
    }
};

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<char> {
    void    operator()( std::ostream& ostr, char t );
};

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<unsigned char> {
    void    operator()( std::ostream& ostr, unsigned char t );
};

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<char const*> {
    void    operator()( std::ostream& ostr, char const* t );
};

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<wchar_t const*> {
    void    operator()( std::ostream& ostr, wchar_t const* t );
};

//____________________________________________________________________________//

namespace tt_detail {

// ************************************************************************** //
// **************              tools classification            ************** //
// ************************************************************************** //

enum check_type {
    CHECK_PRED, 
    CHECK_MSG,
    CHECK_EQUAL,
    CHECK_NE,
    CHECK_LT,
    CHECK_LE,
    CHECK_GT,
    CHECK_GE,
    CHECK_CLOSE,
    CHECK_CLOSE_FRACTION,
    CHECK_SMALL,
    CHECK_BITWISE_EQUAL,
    CHECK_PRED_WITH_ARGS,
    CHECK_EQUAL_COLL
};

enum tool_level {
    WARN, CHECK, REQUIRE, PASS
};

// ************************************************************************** //
// **************                 print_helper                 ************** //
// ************************************************************************** //
// Adds level of indirection to the output operation, allowing us to customize 
// it for types that do not support operator << directly or for any other reason

template<typename T>
struct print_helper_t {
    explicit    print_helper_t( T const& t ) : m_t( t ) {}

    T const&    m_t;
};

//____________________________________________________________________________//

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564)) 
// Borland suffers premature pointer decay passing arrays by reference
template<typename T, std::size_t N >
struct print_helper_t< T[N] > {
    explicit    print_helper_t( T const * t ) : m_t( t ) {}

    T const *   m_t;
};
#endif

//____________________________________________________________________________//

template<typename T>
inline print_helper_t<T> print_helper( T const& t )
{
    return print_helper_t<T>( t );
}

//____________________________________________________________________________//

template<typename T>
inline std::ostream& 
operator<<( std::ostream& ostr, print_helper_t<T> const& ph )
{
    print_log_value<T>()( ostr, ph.m_t );

    return ostr;
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************            TOOL BOX Implementation           ************** //
// ************************************************************************** //

BOOST_TEST_DECL 
bool check_impl( predicate_result const& pr, ::boost::unit_test::lazy_ostream const& assertion_descr,
                 const_string file_name, std::size_t line_num,
                 tool_level tl, check_type ct,
                 std::size_t num_args, ... );

//____________________________________________________________________________//
// This function adds level of indirection, but it makes sure we evaluate predicate 
// arguments only once

#ifndef BOOST_TEST_PROD
#define TEMPL_PARAMS( z, m, dummy ) , typename BOOST_JOIN( Arg, m )

#define FUNC_PARAMS( z, m, dummy )                                                  \
 , BOOST_JOIN( Arg, m ) const& BOOST_JOIN( arg, m )                                 \
 , char const* BOOST_JOIN( BOOST_JOIN( arg, m ), _descr )                           \
/**/

#define PRED_PARAMS( z, m, dummy ) BOOST_PP_COMMA_IF( m ) BOOST_JOIN( arg, m ) 

#define ARG_INFO( z, m, dummy )                                                     \
 , BOOST_JOIN( BOOST_JOIN( arg, m ), _descr )                                       \
 , &static_cast<const unit_test::lazy_ostream&>(unit_test::lazy_ostream::instance() \
        << ::boost::test_tools::tt_detail::print_helper( BOOST_JOIN( arg, m ) ))    \
/**/

#define IMPL_FRWD( z, n, dummy )                                                    \
template<typename Pred                                                              \
         BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), TEMPL_PARAMS, _ )>            \
inline bool                                                                         \
check_frwd( Pred P, unit_test::lazy_ostream const& assertion_descr,                 \
            const_string file_name, std::size_t line_num,                           \
            tool_level tl, check_type ct                                            \
            BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), FUNC_PARAMS, _ )           \
)                                                                                   \
{                                                                                   \
    return                                                                          \
    check_impl( P( BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), PRED_PARAMS, _ ) ), \
                assertion_descr, file_name, line_num, tl, ct,                       \
                BOOST_PP_ADD( n, 1 )                                                \
                BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), ARG_INFO, _ )          \
    );                                                                              \
}                                                                                   \
/**/

#ifndef BOOST_TEST_MAX_PREDICATE_ARITY
#define BOOST_TEST_MAX_PREDICATE_ARITY 5
#endif

BOOST_PP_REPEAT( BOOST_TEST_MAX_PREDICATE_ARITY, IMPL_FRWD, _ )

#undef TEMPL_PARAMS
#undef FUNC_PARAMS
#undef PRED_INFO
#undef ARG_INFO
#undef IMPL_FRWD

#endif

//____________________________________________________________________________//

template <class Left, class Right>
predicate_result equal_impl( Left const& left, Right const& right )
{
    return left == right;
}

//____________________________________________________________________________//

predicate_result        BOOST_TEST_DECL equal_impl( char const* left, char const* right );
inline predicate_result equal_impl( char* left, char const* right ) { return equal_impl( static_cast<char const*>(left), static_cast<char const*>(right) ); }
inline predicate_result equal_impl( char const* left, char* right ) { return equal_impl( static_cast<char const*>(left), static_cast<char const*>(right) ); }
inline predicate_result equal_impl( char* left, char* right )       { return equal_impl( static_cast<char const*>(left), static_cast<char const*>(right) ); }

#if !defined( BOOST_NO_CWCHAR )
predicate_result        BOOST_TEST_DECL equal_impl( wchar_t const* left, wchar_t const* right );
inline predicate_result equal_impl( wchar_t* left, wchar_t const* right ) { return equal_impl( static_cast<wchar_t const*>(left), static_cast<wchar_t const*>(right) ); }
inline predicate_result equal_impl( wchar_t const* left, wchar_t* right ) { return equal_impl( static_cast<wchar_t const*>(left), static_cast<wchar_t const*>(right) ); }
inline predicate_result equal_impl( wchar_t* left, wchar_t* right )       { return equal_impl( static_cast<wchar_t const*>(left), static_cast<wchar_t const*>(right) ); }
#endif

//____________________________________________________________________________//

struct equal_impl_frwd {
    template <typename Left, typename Right>
    inline predicate_result
    call_impl( Left const& left, Right const& right, mpl::false_ ) const
    {
        return equal_impl( left, right );
    }

    template <typename Left, typename Right>
    inline predicate_result
    call_impl( Left const& left, Right const& right, mpl::true_ ) const
    {
        return (*this)( right, &left[0] );
    }

    template <typename Left, typename Right>
    inline predicate_result
    operator()( Left const& left, Right const& right ) const
    {
        typedef typename is_array<Left>::type left_is_array;
        return call_impl( left, right, left_is_array() );
    }
};

//____________________________________________________________________________//

struct ne_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return !equal_impl_frwd()( left, right );
    }
};

//____________________________________________________________________________//

struct lt_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return left < right;
    }
};

//____________________________________________________________________________//

struct le_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return left <= right;
    }
};

//____________________________________________________________________________//

struct gt_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return left > right;
    }
};

//____________________________________________________________________________//

struct ge_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return left >= right;
    }
};

//____________________________________________________________________________//

struct equal_coll_impl {
    template <typename Left, typename Right>
    predicate_result operator()( Left left_begin, Left left_end, Right right_begin, Right right_end )
    {
        predicate_result    pr( true );
        std::size_t         pos = 0;

        for( ; left_begin != left_end && right_begin != right_end; ++left_begin, ++right_begin, ++pos ) {
            if( *left_begin != *right_begin ) {
                pr = false;
                pr.message() << "\nMismatch in a position " << pos << ": "  << *left_begin << " != " << *right_begin;
            }
        }

        if( left_begin != left_end ) {
            std::size_t r_size = pos;
            while( left_begin != left_end ) {
                ++pos;
                ++left_begin;
            }

            pr = false;
            pr.message() << "\nCollections size mismatch: " << pos << " != " << r_size;
        }

        if( right_begin != right_end ) {
            std::size_t l_size = pos;
            while( right_begin != right_end ) {
                ++pos;
                ++right_begin;
            }

            pr = false;
            pr.message() << "\nCollections size mismatch: " << l_size << " != " << pos;
        }

        return pr;
    }
};

//____________________________________________________________________________//

struct bitwise_equal_impl {
    template <class Left, class Right>
    predicate_result    operator()( Left const& left, Right const& right )
    {
        predicate_result    pr( true );

        std::size_t left_bit_size  = sizeof(Left)*CHAR_BIT;
        std::size_t right_bit_size = sizeof(Right)*CHAR_BIT;

        static Left const leftOne( 1 );
        static Right const rightOne( 1 );

        std::size_t total_bits = left_bit_size < right_bit_size ? left_bit_size : right_bit_size;

        for( std::size_t counter = 0; counter < total_bits; ++counter ) {
            if( ( left & ( leftOne << counter ) ) != ( right & ( rightOne << counter ) ) ) {
                pr = false;
                pr.message() << "\nMismatch in a position " << counter;
            }
        }

        if( left_bit_size != right_bit_size ) {
            pr = false;
            pr.message() << "\nOperands bit sizes mismatch: " << left_bit_size << " != " << right_bit_size;
        }

        return pr;
    }
};

//____________________________________________________________________________//

bool BOOST_TEST_DECL is_defined_impl( const_string symbol_name, const_string symbol_value );

//____________________________________________________________________________//

} // namespace tt_detail

} // namespace test_tools

namespace test_toolbox = test_tools;

} // namespace boost

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

#endif // BOOST_TEST_TEST_TOOLS_HPP_012705GER
