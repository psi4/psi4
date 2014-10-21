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
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** GLOBALS.H
**
** List of all the global data used by the program
**
** Note that these are given as "extern", so they must be defined
** in the main program!
**
** C. David Sherrill
** University of California, Berkeley
*/

#ifndef _psi_src_bin_detcas_globals_h
#define _psi_src_bin_detcas_globals_h

#include "calcinfo.h"
#include "MCSCF_params.h"

extern "C" {
  extern FILE *infile, *outfile;
  extern char *psi_file_prefix;
}

namespace psi { namespace detcas {

extern struct calcinfo CalcInfo;
extern struct params MCSCF_Parameters;
extern int *ioff;

}} // end namespace psi::detcas

#endif // header guard

