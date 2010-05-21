// Boost detail/extended_integer.hpp header file  ----------------------------//

// (C) Copyright Daryle Walker 2008.  Distributed under the Boost Software
// License, Version 1.0.  (See the accompanying file LICENSE_1_0.txt or a copy
// at <http://www.boost.org/LICENSE_1_0.txt>.)

// Encapsulates the double-long and __int64 type families as a single family,
// as they are mutually exclusive.

/** \file
    \brief  Common definition of extended integer types.

    Instead of other Boost headers making separate \#defines for the double-long
    and __int64 type families, since they're mutually exclusive, make a single
    set of types and macros for the family that exists (if either).
 */

#ifndef BOOST_DETAIL_EXTENDED_INTEGER_HPP
#define BOOST_DETAIL_EXTENDED_INTEGER_HPP

#include <boost/config.hpp>  // for BOOST_HAS_LONG_LONG and BOOST_HAS_MS_INT64

#include <climits>  // for CHAR_BIT, etc.


namespace boost
{
namespace detail
{


//  Extended integer type macro and alias definitions  -----------------------//

// (Unsigned) long long family
#ifdef BOOST_HAS_LONG_LONG

// Existence
#define BOOST_HAS_XINT  1

// Extents
#ifdef ULLONG_MAX
#define BOOST_XINT_MAX    LLONG_MAX
#define BOOST_XINT_MIN    LLONG_MIN
#define BOOST_UXINT_MAX  ULLONG_MAX
#elif defined(ULONG_LONG_MAX)
#define BOOST_XINT_MAX    LONG_LONG_MAX
#define BOOST_XINT_MIN    LONG_LONG_MIN
#define BOOST_UXINT_MAX  ULONG_LONG_MAX
#elif defined(ULONGLONG_MAX)
#define BOOST_XINT_MAX    LONGLONG_MAX
#define BOOST_XINT_MIN    LONGLONG_MIN
#define BOOST_UXINT_MAX  ULONGLONG_MAX
#elif defined(_LLONG_MAX) && defined(_C2)
#define BOOST_XINT_MAX    _LLONG_MAX
#define BOOST_XINT_MIN    (-_LLONG_MAX - _C2)
#define BOOST_UXINT_MAX  _ULLONG_MAX
#else  // guess
// Sometimes we get the double-long types without the corresponding constants,
// e.g. GCC in "-ansi" mode.  In this case, we'll just have to work out the
// values ourselves.  (Here we assume a two's complement representation.)
#define BOOST_XINT_MIN   (1LL << (sizeof(::boost::long_long_type) * CHAR_BIT - 1))
#define BOOST_XINT_MAX   (~ BOOST_XINT_MIN)
#define BOOST_UXINT_MAX  (~ 0uLL)
#endif

// Types
typedef ::boost:: long_long_type   xint_t;
typedef ::boost::ulong_long_type  uxint_t;

// (Unsigned) __int64 family
#elif defined(BOOST_HAS_MS_INT64)

// Existence
#define BOOST_HAS_XINT  1

// Extents
#ifdef _UI64_MAX
#define BOOST_XINT_MAX    _I64_MAX
#define BOOST_XINT_MIN    _I64_MIN
#define BOOST_UXINT_MAX  _UI64_MAX
#else  // guess
// The types are exactly 2's-compl. 64-bit, so we'll enter the values directly.
#define BOOST_XINT_MAX   0x7FFFFFFFFFFFFFFFi64
#define BOOST_XINT_MIN   0x8000000000000000i64
#define BOOST_UXINT_MAX  0xFFFFFFFFFFFFFFFFui64
#endif

// Types
typedef          __int64   xint_t;
typedef unsigned __int64  uxint_t;

// Neither
#else

// Non-existence
#define BOOST_HAS_XINT  0

// Dummy extents
#define BOOST_XINT_MAX    LONG_MAX
#define BOOST_XINT_MIN    LONG_MIN
#define BOOST_UXINT_MAX  ULONG_MAX

// Dummy types
typedef   signed long   xint_t;
typedef unsigned long  uxint_t;

#endif  // defined(BOOST_HAS_LONG_LONG)/defined(BOOST_HAS_MS_INT64)/else

/** \def  BOOST_HAS_XINT

    \brief  Flag for extended integer types.

    Indicates the presence of one of the two common extended integer type
    families, either (<code>unsigned</code>) <code>long long</code> or
    (<code>unsigned</code>) <code>__int64</code>.  \c BOOST_HAS_XINT is \c 1 if
    either type family is defined, and \c 0 if neither is.
 */

/** \def  BOOST_XINT_MAX

    \brief  Maximum value for the signed extended integer type.

    \pre  \c BOOST_HAS_XINT is \c \#defined to be \c 1.

    Macro constant representing the largest value the signed extended integer
    type supports.  Its composition may be another macro, an expression, or a
    literal.  Defaulted to \c LONG_MAX if \c BOOST_HAS_XINT is zero.
 */
/** \def  BOOST_XINT_MIN

    \brief  Minimum value for the signed extended integer type.

    \pre  \c BOOST_HAS_XINT is \c \#defined to be \c 1.

    Macro constant representing the smallest value the signed extended integer
    type supports.  Its composition may be another macro, an expression, or a
    literal.  Defaulted to \c LONG_MIN if \c BOOST_HAS_XINT is zero.
 */
/** \def  BOOST_UXINT_MAX

    \brief  Maximum value for the unsigned extended integer type.

    \pre  \c BOOST_HAS_XINT is \c \#defined to be \c 1.

    Macro constant representing the largest value the unsigned extended integer
    type supports.  Its composition may be another macro, an expression, or a
    literal.  Defaulted to \c ULONG_MAX if \c BOOST_HAS_XINT is zero.  (Use
    \c 0u for the type's minimum value.)
 */

/** \typedef  signed long  boost::detail::xint_t

    \brief  Alias for the signed extended integer type.

    \pre  \c BOOST_HAS_XINT is \c \#defined to be \c 1.

    Alias representing the signed extended integer type, no matter which type
    family it came from.  Defaulted to <code>signed long</code> if
    \c BOOST_HAS_XINT is zero.
 */
/** \typedef  unsigned long  ::boost::detail::uxint_t

    \brief  Alias for the signed extended integer type.

    \pre  \c BOOST_HAS_XINT is \c \#defined to be \c 1.

    Alias representing the unsigned extended integer type, no matter which type
    family it came from.  Defaulted to <code>unsigned long</code> if
    \c BOOST_HAS_XINT is zero.
 */


}  // namespace detail
}  // namespace boost


#endif  // BOOST_DETAIL_EXTENDED_INTEGER_HPP
