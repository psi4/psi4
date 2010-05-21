
//  (C) Copyright John maddock 1999. 
//  (C) David Abrahams 2002.  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// use this header as a workaround for missing <limits>

//  See http://www.boost.org/libs/compatibility/index.html for documentation.

#ifndef BOOST_LIMITS
#define BOOST_LIMITS

#include <boost/config.hpp>

#ifdef BOOST_NO_LIMITS
# include <boost/detail/limits.hpp>
#else
# include <limits>
#endif

#include <boost/detail/extended_integer.hpp>  // for BOOST_HAS_XINT, etc.

#if (defined(BOOST_HAS_LONG_LONG) && defined(BOOST_NO_LONG_LONG_NUMERIC_LIMITS)) \
      || (defined(BOOST_HAS_MS_INT64) && defined(BOOST_NO_MS_INT64_NUMERIC_LIMITS))
// Add missing specializations for numeric_limits:
#if !defined(BOOST_HAS_XINT) || !(BOOST_HAS_XINT)
#error "Shouldn't have gotten here based on preceeding preprocessor statements"
#endif
#define BOOST_LLT    ::boost::detail::xint_t
#define BOOST_ULLT  ::boost::detail::uxint_t

#include <climits>  // for CHAR_BIT

namespace std
{
  template<>
  class numeric_limits< BOOST_LLT > 
  {
   public:

      BOOST_STATIC_CONSTANT(bool, is_specialized = true);
      static BOOST_LLT min BOOST_PREVENT_MACRO_SUBSTITUTION (){ return BOOST_XINT_MIN; }
      static BOOST_LLT max BOOST_PREVENT_MACRO_SUBSTITUTION (){ return BOOST_XINT_MAX; }
      BOOST_STATIC_CONSTANT(int, digits = sizeof(BOOST_LLT) * CHAR_BIT -1);
      BOOST_STATIC_CONSTANT(int, digits10 = (CHAR_BIT * sizeof (BOOST_LLT) - 1) * 301L / 1000);
      BOOST_STATIC_CONSTANT(bool, is_signed = true);
      BOOST_STATIC_CONSTANT(bool, is_integer = true);
      BOOST_STATIC_CONSTANT(bool, is_exact = true);
      BOOST_STATIC_CONSTANT(int, radix = 2);
      static BOOST_LLT epsilon() throw() { return 0; };
      static BOOST_LLT round_error() throw() { return 0; };

      BOOST_STATIC_CONSTANT(int, min_exponent = 0);
      BOOST_STATIC_CONSTANT(int, min_exponent10 = 0);
      BOOST_STATIC_CONSTANT(int, max_exponent = 0);
      BOOST_STATIC_CONSTANT(int, max_exponent10 = 0);

      BOOST_STATIC_CONSTANT(bool, has_infinity = false);
      BOOST_STATIC_CONSTANT(bool, has_quiet_NaN = false);
      BOOST_STATIC_CONSTANT(bool, has_signaling_NaN = false);
      BOOST_STATIC_CONSTANT(bool, has_denorm = false);
      BOOST_STATIC_CONSTANT(bool, has_denorm_loss = false);
      static BOOST_LLT infinity() throw() { return 0; };
      static BOOST_LLT quiet_NaN() throw() { return 0; };
      static BOOST_LLT signaling_NaN() throw() { return 0; };
      static BOOST_LLT denorm_min() throw() { return 0; };

      BOOST_STATIC_CONSTANT(bool, is_iec559 = false);
      BOOST_STATIC_CONSTANT(bool, is_bounded = true);
      BOOST_STATIC_CONSTANT(bool, is_modulo = true);

      BOOST_STATIC_CONSTANT(bool, traps = false);
      BOOST_STATIC_CONSTANT(bool, tinyness_before = false);
      BOOST_STATIC_CONSTANT(float_round_style, round_style = round_toward_zero);
      
  };

  template<>
  class numeric_limits< BOOST_ULLT > 
  {
   public:

      BOOST_STATIC_CONSTANT(bool, is_specialized = true);
      static BOOST_ULLT min BOOST_PREVENT_MACRO_SUBSTITUTION (){ return (BOOST_ULLT) 0u; }
      static BOOST_ULLT max BOOST_PREVENT_MACRO_SUBSTITUTION (){ return BOOST_UXINT_MAX; }
      BOOST_STATIC_CONSTANT(int, digits = sizeof(BOOST_LLT) * CHAR_BIT);
      BOOST_STATIC_CONSTANT(int, digits10 = (CHAR_BIT * sizeof (BOOST_LLT)) * 301L / 1000);
      BOOST_STATIC_CONSTANT(bool, is_signed = false);
      BOOST_STATIC_CONSTANT(bool, is_integer = true);
      BOOST_STATIC_CONSTANT(bool, is_exact = true);
      BOOST_STATIC_CONSTANT(int, radix = 2);
      static BOOST_ULLT epsilon() throw() { return 0; };
      static BOOST_ULLT round_error() throw() { return 0; };

      BOOST_STATIC_CONSTANT(int, min_exponent = 0);
      BOOST_STATIC_CONSTANT(int, min_exponent10 = 0);
      BOOST_STATIC_CONSTANT(int, max_exponent = 0);
      BOOST_STATIC_CONSTANT(int, max_exponent10 = 0);

      BOOST_STATIC_CONSTANT(bool, has_infinity = false);
      BOOST_STATIC_CONSTANT(bool, has_quiet_NaN = false);
      BOOST_STATIC_CONSTANT(bool, has_signaling_NaN = false);
      BOOST_STATIC_CONSTANT(bool, has_denorm = false);
      BOOST_STATIC_CONSTANT(bool, has_denorm_loss = false);
      static BOOST_ULLT infinity() throw() { return 0; };
      static BOOST_ULLT quiet_NaN() throw() { return 0; };
      static BOOST_ULLT signaling_NaN() throw() { return 0; };
      static BOOST_ULLT denorm_min() throw() { return 0; };

      BOOST_STATIC_CONSTANT(bool, is_iec559 = false);
      BOOST_STATIC_CONSTANT(bool, is_bounded = true);
      BOOST_STATIC_CONSTANT(bool, is_modulo = true);

      BOOST_STATIC_CONSTANT(bool, traps = false);
      BOOST_STATIC_CONSTANT(bool, tinyness_before = false);
      BOOST_STATIC_CONSTANT(float_round_style, round_style = round_toward_zero);
      
  };
}
#endif 

#endif

