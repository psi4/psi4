/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */


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
