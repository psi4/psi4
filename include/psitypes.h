
#ifndef _psi_include_psitypes_h_
#define _psi_include_psitypes_h_

#include <psiconfig.h>

/*
  Define Psi data types here
*/

/* default PSI floating-point type */
typedef double PSI_FLOAT;

/* default PSI 64-bit integer */
#ifdef HAVE_STDINT_H

#include <stdint.h>
typedef int_least64_t PSI_INT_LEAST64; 

#else

#include <climits>

#if defined(ULONGLONG_MAX) && !defined(ULLONG_MAX)
#    define ULLONG_MAX ULONGLONG_MAX
#endif

#if defined(ULLONG_MAX) && !defined(ULONGLONG_MAX)
#    define ULONGLONG_MAX ULLONG_MAX
#endif

# ifdef ULLONG_MAX
#   if ULONGLONG_MAX == (0xffffffffffffffffuLL) /* uLL reqd for xlC */
     typedef long long PSI_INT_LEAST64; 
#   else
#     error defaults not correct; you must hand modify psitypes.h
#   endif
# elif ULONG_MAX != 0xffffffff

#   if ULONG_MAX == 18446744073709551615 /* 2**64 - 1 */
     typedef long PSI_INT_LEAST64;
#   else
#     error defaults not correct; you must hand modify scint.h
#   endif
# else /* assume no 64-bit integers */
#   error 64 bit integer types are required
# endif

#endif /* HAVE_STDINT_H */

#endif
