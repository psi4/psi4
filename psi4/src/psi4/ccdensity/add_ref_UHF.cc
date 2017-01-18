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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libiwl/iwl.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* add_ref_UHF(): This function adds the reference contributions to
** the one- and two-particle density matrices.  These contributions
** are simply the prefactors in front of the one- and two-electron
** intgegrals, respectively, in the UHF energy expression.  In the
** case of the two-pdm, however, care must be taken that only the
** permutationally unique elements be written to disk.
**
** In spin-orbitals with Mulliken-order two-electron integrals, the
** three spin contributions to the SCF energy are:
**
** E_AA(SCF) = sum_I h_II + 1/2 sum_IJ ([II|JJ] - [IJ|IJ])
** E_BB(SCF) = sum_i h_ii + 1/2 sum_ij ([ii|jj] - [ij|ij])
** E_AB(SCF) = 1/2 sum_Ij [II|jj]
**
** I use QT-standard ordering for the indices in these expressions.
*/

void add_ref_UHF(struct iwlbuf *AA, struct iwlbuf *BB, struct iwlbuf *AB)
{
  int i,j;
  int nfzc, nclsd, nopen;

  nfzc = moinfo.nfzc;
  nclsd = moinfo.nclsd;
  nopen = moinfo.nopen;

  /*** One-electron components ***/

  /* alpha-occ block */
  for(i=0; i < (nfzc + nclsd + nopen); i++)
    moinfo.opdm_a[i][i] += 1.0;

  /* beta-occ block */
  for(i=0; i < (nfzc + nclsd); i++)
    moinfo.opdm_b[i][i] += 1.0;

  /*** Two-electron components ***/

  /* AA */
  for(i=0; i < (nfzc + nclsd + nopen); i++) {
    for(j=0; j < i; j++) {
      iwl_buf_wrt_val(AA, i, i, j, j, 0.5, 0, "outfile", 0);

      iwl_buf_wrt_val(AA, i, j, i, j,-0.25, 0, "outfile", 0);
      iwl_buf_wrt_val(AA, j, i, j, i,-0.25, 0, "outfile", 0);
      iwl_buf_wrt_val(AA, i, j, j, i,-0.25, 0, "outfile", 0);
    }
  }

  /* BB */
  for(i=0; i < (nfzc + nclsd); i++) {
    for(j=0; j < i; j++) {
      iwl_buf_wrt_val(BB, i, i, j, j, 0.5, 0, "outfile", 0);

      iwl_buf_wrt_val(BB, i, j, i, j,-0.25, 0, "outfile", 0);
      iwl_buf_wrt_val(BB, j, i, j, i,-0.25, 0, "outfile", 0);
      iwl_buf_wrt_val(BB, i, j, j, i,-0.25, 0, "outfile", 0);
    }
  }

  /* AB */
  for(i=0; i < (nfzc + nclsd + nopen); i++)
    for(j=0; j < (nfzc + nclsd); j++)
      iwl_buf_wrt_val(AB, i, i, j, j, 1.0, 0, "outfile", 0);

}

}} // namespace psi::ccdensity
