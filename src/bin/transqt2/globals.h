/*! \file
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/

#ifndef _psi_src_bin_transqt2_globals_h
#define _psi_src_bin_transqt2_globals_h

#include <psifiles.h>
#include <ccfiles.h>

#include "MOInfo.h"
#include "Params.h"

using namespace psi;

namespace psi {
  namespace transqt2 {

#define IOFF_MAX 32641

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN struct MOInfo moinfo;
EXTERN struct Params params;

#define INDEX(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))


  } //namespace transqt2
} //namespace psi

#endif // _psi_src_bin_transqt2_globals_h
