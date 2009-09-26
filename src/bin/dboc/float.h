/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#ifndef _psi3_DBOC_float_h_
#define _psi3_DBOC_float_h_

#include "defines.h"

#if LONG_DOUBLE
typedef long double FLOAT;
#else
typedef double FLOAT;
#endif

#if LONG_DOUBLE
#  ifdef AIX
#    define FABS fabsl
#  else
#    define FABS fabs
#  endif
#else
#  define FABS fabs
#endif

#endif
