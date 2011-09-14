/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_mp2_globals_h_
#define _psi_src_bin_mp2_globals_h_

#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <ccfiles.h>
#include <psi4-dec.h>
#include "moinfo.h"
#include "params.h"

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

#define MAXIOFF 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

namespace psi{ namespace mp2{
// BJM removing the following four lines
//extern "C" {
//  EXTERN FILE *infile, *outfile;
//  EXTERN char *psi_file_prefix;
//}
EXTERN struct moinfo mo;
EXTERN struct params params;
EXTERN int* ioff;

}} // namespaces
#endif /* Header guard */
