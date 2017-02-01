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
**  \ingroup DETCI
**  \brief Form OV arrays of Bendazzoli and Evangelisti, JCP 98, 3141 (1993)
**
** David Sherrill
** University of Georgia
** 8 April 1996
**
*/

#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

/*
** FORM_OV()
** This will only work for Full CI's right now (where Parameters_->fci=true)
*/
void CIWavefunction::form_ov()
{

   int i, j, nirreps, norbs;
   int irrep, strnum, strsym, cnt=0;
   int fullij, idx, ovcnt;
   struct stringwr *strlist;
   int signmask,nsignmask;


   /* bitwise sign stuff */

   signmask = 1 << (sizeof(int)*8-1);
   nsignmask = ~signmask;


   /* allocate memory for OV[list][fullij][string] */

   norbs = CalcInfo_->num_ci_orbs;
   nirreps = AlphaG_->nirreps;
   OV_ = (int ***) malloc (sizeof(int **) * nirreps);
   for (i=0; i<nirreps; i++) {
      OV_[i] = (int **) malloc (sizeof(int *) * norbs * norbs);
      for (j=0; j<norbs*norbs; j++) {
         OV_[i][j] = (int *) malloc (sizeof(int) * AlphaG_->max_str_per_irrep+1);
         OV_[i][j][0] = 0;
         }
      }


   /* now fill up OV by walking through the stringwr lists */

   for (irrep=0; irrep < nirreps; irrep++) {
      strnum = AlphaG_->sg[irrep][0].num_strings;
      cnt=0;
      strlist = alplist_[irrep];
      while (cnt != strnum) {
         for (strsym=0; strsym < nirreps; strsym++) {
            for (i=0; i<strlist->cnt[strsym]; i++) {
               fullij = strlist->oij[strsym][i];
               /* idx = cnt + 1; */
               idx = cnt;
               if (strlist->sgn[strsym][i] != 1) idx = idx | signmask;
               ovcnt = OV_[irrep][fullij][0];
               ovcnt++;
               OV_[irrep][fullij][ovcnt] = idx;
               OV_[irrep][fullij][0] = ovcnt;
               }
            }
         strlist++;
         cnt++;
         }
      }


   /* print out the OV data */

   if (print_ > 3) {
      for (irrep=0; irrep < nirreps; irrep++) {
         for (fullij=0; fullij<norbs*norbs; fullij++) {
            outfile->Printf( "OV[irrep=%d][oij=%d]:  ", irrep, fullij);
            for (i=0; i<OV_[irrep][fullij][0]; i++) {
               idx = OV_[irrep][fullij][i+1];
               outfile->Printf( "%c", (idx & signmask) ? '-' : '+');
               idx = idx & nsignmask;
               outfile->Printf( "%2d ", idx);
               }
            outfile->Printf( "\n");
            }
         }
      }


}

}} // namespace psi::detci
