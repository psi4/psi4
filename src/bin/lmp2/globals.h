/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_lmp2_globals_h_
#define _psi_src_bin_lmp2_globals_h_

#include "class.h"

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

#define MAXIOFF 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

#endif /* Header guard */
