/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_detci_globals_h
#define _psi_src_bin_detci_globals_h


#include <psi4-dec.h>
#include <string>

/* nice stuff to extern or not to extern properly */
#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

namespace psi { namespace detci {

//EXTERN int errcod;
EXTERN struct calcinfo CalcInfo;
EXTERN struct params Parameters;
EXTERN struct H_zero_block H0block;
//EXTERN int *ioff;
//EXTERN struct ci_blks CIblks;
//EXTERN struct olsen_graph *AlphaG;
//EXTERN struct olsen_graph *BetaG;
//EXTERN struct graph_set *AlphaGraph;
//EXTERN struct graph_set *BetaGraph;
//EXTERN int ***OV;
//EXTERN int **s1_contrib, **s2_contrib, **s3_contrib;
//EXTERN double *tmp_ras_array;
//EXTERN struct detci_timings detci_time;

}} // namespace psi::detci

#endif // header guard