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
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/

#ifndef _psi_src_bin_transqt2_globals_h
#define _psi_src_bin_transqt2_globals_h

#include <psifiles.h>

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
