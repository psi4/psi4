/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/

#ifndef _psi3_bin_input_float_h_
#define _psi3_bin_input_float_h_

#include "defines.h"

#if LONG_DOUBLE
#define FLOAT long double
#else
#define FLOAT double
#endif

#if LONG_DOUBLE
#  define FABS fabsl
#  define EXP expl
#  define SQRT sqrtl
#else
#  define FABS fabs
#  define EXP exp
#  define SQRT sqrt
#endif

#endif // header guard
