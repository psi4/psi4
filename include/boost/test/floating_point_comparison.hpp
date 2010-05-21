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
//  Description : defines algoirthms for comparing 2 floating point values
// ***************************************************************************

#ifndef BOOST_TEST_FLOATING_POINT_COMPARISON_HPP_071894GER
#define BOOST_TEST_FLOATING_POINT_COMPARISON_HPP_071894GER

// Boost.Test
#include <boost/test/detail/global_typedef.hpp>
#include <boost/test/predicate_result.hpp>

// Boost
#include <boost/limits.hpp>  // for std::numeric_limits
#include <boost/numeric/conversion/conversion_traits.hpp> // for numeric::conversion_traits
#include <boost/static_assert.hpp>
#include <boost/assert.hpp>

// STL
#include <iosfwd>

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

namespace boost {

namespace math { namespace fpc {

// ************************************************************************** //
// **************                 fpc::strength                ************** //
// ************************************************************************** //

enum strength {
    FPC_STRONG, // "Very close"   - equation 1' in docs, the default
    FPC_WEAK    // "Close enough" - equation 2' in docs.
};

// ************************************************************************** //
// **************                    details                   ************** //
// ************************************************************************** //

namespace fpc_detail {

// FPT is Floating-Point Type: float, double, long double or User-Defined.
template<typename FPT>
inline FPT
fpt_abs( FPT fpv ) 
{
    return fpv < static_cast<FPT>(0) ? -fpv : fpv;
}

//____________________________________________________________________________//

template<typename FPT>
struct fpt_limits {
    static FPT    min_value()
    {
        return std::numeric_limits<FPT>::is_specialized
                    ? (std::numeric_limits<FPT>::min)()
                    : 0;
    }
    static FPT    max_value()
    {
        return std::numeric_limits<FPT>::is_specialized
                    ? (std::numeric_limits<FPT>::max)()
                    : static_cast<FPT>(1000000); // for the our purposes it doesn't really matter what value is returned here
    }
};

//____________________________________________________________________________//

// both f1 and f2 are unsigned here
template<typename FPT>
inline FPT
safe_fpt_division( FPT f1, FPT f2 )
{
    // Avoid overflow.
    if( (f2 < static_cast<FPT>(1))  && (f1 > f2*fpt_limits<FPT>::max_value()) )
        return fpt_limits<FPT>::max_value();

    // Avoid underflow.
    if( (f1 == static_cast<FPT>(0)) ||
        ((f2 > static_cast<FPT>(1)) && (f1 < f2*fpt_limits<FPT>::min_value())) )
        return static_cast<FPT>(0);

    return f1/f2;
}

//____________________________________________________________________________//

} // namespace fpc_detail

// ************************************************************************** //
// **************         tolerance presentation types         ************** //
// ************************************************************************** //

template<typename ToleranceType>
struct tolerance_traits {
    template<typename FPT>
    static ToleranceType    actual_tolerance( FPT fraction_tolerance )
    {
        return static_cast<ToleranceType>( fraction_tolerance );
    } 
    template<typename FPT>
    static FPT              fraction_tolerance( ToleranceType tolerance )
    {
        return static_cast<FPT>(tolerance);
    } 
};

//____________________________________________________________________________//

template<typename FPT>
struct percent_tolerance_t {
    explicit    percent_tolerance_t( FPT v ) : m_value( v ) {}

    FPT m_value;
};

//____________________________________________________________________________//

template<typename FPT>
struct tolerance_traits<percent_tolerance_t<FPT> > {
    template<typename FPT2>
    static percent_tolerance_t<FPT> actual_tolerance( FPT2 fraction_tolerance )
    {
        return percent_tolerance_t<FPT>( fraction_tolerance * static_cast<FPT2>(100.) );
    }

    template<typename FPT2>
    static FPT2 fraction_tolerance( percent_tolerance_t<FPT> tolerance )
    {
        return static_cast<FPT2>(tolerance.m_value)*static_cast<FPT2>(0.01); 
    }
};

//____________________________________________________________________________//

template<typename FPT>
std::ostream& operator<<( std::ostream& out, percent_tolerance_t<FPT> t )
{
    return out << t.m_value;
}

//____________________________________________________________________________//

template<typename FPT>
inline percent_tolerance_t<FPT>
percent_tolerance( FPT v )
{
    return percent_tolerance_t<FPT>( v );
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************             close_at_tolerance               ************** //
// ************************************************************************** //

template<typename FPT>
class close_at_tolerance {
public:
    // Public typedefs
    typedef bool result_type;

    // Constructor
    template<typename ToleranceType>
    explicit    close_at_tolerance( ToleranceType tolerance, fpc::strength fpc_strength = FPC_STRONG ) 
    : m_fraction_tolerance( tolerance_traits<ToleranceType>::template fraction_tolerance<FPT>( tolerance ) )
    , m_strength( fpc_strength )
    {
        BOOST_ASSERT( m_fraction_tolerance >= 0 ); // no reason for tolerance to be negative
    }

    // Access methods
    FPT                 fraction_tolerance() const  { return m_fraction_tolerance; }
    fpc::strength       strength() const            { return m_strength; }
    FPT                 failed_fraction() const     { return m_failed_fraction; }

    // Action method
    bool                operator()( FPT left, FPT right ) const
    {
        FPT diff              = fpc_detail::fpt_abs( left - right );
        FPT fraction_of_right = fpc_detail::safe_fpt_division( diff, fpc_detail::fpt_abs( right ) );
        FPT fraction_of_left  = fpc_detail::safe_fpt_division( diff, fpc_detail::fpt_abs( left ) );
        
        bool res( m_strength == FPC_STRONG
            ? (fraction_of_right <= m_fraction_tolerance && fraction_of_left <= m_fraction_tolerance) 
            : (fraction_of_right <= m_fraction_tolerance || fraction_of_left <= m_fraction_tolerance) );

        if( !res )
            m_failed_fraction = (fraction_of_right > m_fraction_tolerance ? fraction_of_right : fraction_of_left);

        return res;
    }

private:
    // Data members
    FPT                 m_fraction_tolerance;
    fpc::strength       m_strength;
	mutable FPT         m_failed_fraction;
};

// ************************************************************************** //
// **************                 is_close_to                  ************** //
// ************************************************************************** //

template<typename FPT1, typename FPT2, typename ToleranceType>
bool
is_close_to( FPT1 left, FPT2 right, ToleranceType tolerance )
{
    // deduce "better" type from types of arguments being compared
    // if one type is floating and the second integral we use floating type and 
    // value of integral type is promoted to the floating. The same for float and double
    // But we don't want to compare two values of integral types using this tool.
    typedef typename numeric::conversion_traits<FPT1,FPT2>::supertype FPT;
    BOOST_STATIC_ASSERT( !is_integral<FPT>::value );

    return fpc::close_at_tolerance<FPT>( tolerance, FPC_STRONG )( left, right );
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************             close_at_tolerance               ************** //
// ************************************************************************** //

template<typename FPT>
class small_with_tolerance {
public:
    // Public typedefs
    typedef bool result_type;

    // Constructor
    explicit    small_with_tolerance( FPT tolerance ) 
    : m_tolerance( tolerance )
    {
        BOOST_ASSERT( m_tolerance >= 0 ); // no reason for the tolerance to be negative
    }

    // Action method
    bool        operator()( FPT fpv ) const
    {
        return fpc::fpc_detail::fpt_abs( fpv ) < m_tolerance;
    }

private:
    // Data members
    FPT         m_tolerance;
};

// ************************************************************************** //
// **************                  is_small                    ************** //
// ************************************************************************** //

template<typename FPT>
bool
is_small( FPT fpv, FPT tolerance )
{
    return small_with_tolerance<FPT>( tolerance )( fpv );
}

//____________________________________________________________________________//

} // namespace fpc
} // namespace math

namespace test_tools {

namespace fpc = math::fpc;

// ************************************************************************** //
// **************               check_is_close                 ************** //
// ************************************************************************** //

struct BOOST_TEST_DECL check_is_close_t {
    // Public typedefs
    typedef bool result_type;

    template<typename FPT1, typename FPT2, typename ToleranceType>
    predicate_result
    operator()( FPT1 left, FPT2 right, ToleranceType tolerance ) const
    {
        predicate_result pr( fpc::is_close_to( left, right, tolerance ) );

        if( !pr )
            pr.message() << tolerance;

        return pr;
    }
};

namespace {
check_is_close_t const& check_is_close = unit_test::ut_detail::static_constant<check_is_close_t>::value;
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************               check_is_small                 ************** //
// ************************************************************************** //

struct BOOST_TEST_DECL check_is_small_t {
    // Public typedefs
    typedef bool result_type;

    template<typename FPT>
    bool
    operator()( FPT fpv, FPT tolerance ) const
    {
        return fpc::is_small( fpv, tolerance );
    }
};

namespace {
check_is_small_t const& check_is_small = unit_test::ut_detail::static_constant<check_is_small_t>::value;
}

//____________________________________________________________________________//

} // namespace test_tools

} // namespace boost

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

#endif // BOOST_FLOATING_POINT_COMAPARISON_HPP_071894GER
