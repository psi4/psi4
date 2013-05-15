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

/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_mp2_globals_h_
#define _psi_src_bin_mp2_globals_h_

#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
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
