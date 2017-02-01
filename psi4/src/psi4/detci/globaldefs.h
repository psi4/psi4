/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here
*/
/*
** GLOBALDEFS.H
**
** Global defines for DETCAS
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
**
*/

#ifndef _psi_src_bin_detcas_globaldefs_h
#define _psi_src_bin_detcas_globaldefs_h

namespace psi { namespace detcas {

#define MAX_RAS_SPACES 4
#define IOFF_MAX       50604
#ifdef INDEX
#undef INDEX
#endif
#define INDEX(i,j) ( (i>j) ? (ioff[(i)] + (j)): (ioff[(j)] + (i)) )
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define MAX_COMMENT 10

}} // end namespace psi::detcas

#endif // header guard