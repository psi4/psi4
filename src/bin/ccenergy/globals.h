/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/

#ifndef _psi_src_bin_ccenergy_globals_h
#define _psi_src_bin_ccenergy_globals_h

#include <ccfiles.h>
#include <psi4-dec.h>
#include <string>
#include <libciomr/libciomr.h>

namespace psi {
extern FILE* outfile;
namespace ccenergy {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

#ifdef DMALLOC
#include <dmalloc.h>
#endif

/* #define TIME_CCENERGY */

EXTERN int *ioff;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Local local;

}} // namespace psi::ccenergy

#endif // _psi_src_bin_ccenergy_globals_h
