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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

/* sort_pert(): Sorts the specified MO-basis one-electron property
** integrals into CC ordering for use in building the
** similarity-transformed integrals and certain components of the
** total linear response function.
**
** NB: Some integrals are antisymmetric (e.g. L or P integrals), and
** others are symmetric (e.g. Mu integrals), so we must be careful in
** this and subsequent routines.
**
** TDC, 10/05
*/

void sort_pert(const char *pert, double **pertints, int irrep)
{
  int p, q, Gp, Gq, P, Q, i;
  dpdfile2 f;
  char prefix[32], lbl[32];

  sprintf(lbl, "%s_IJ", pert);
  global_dpd_->file2_init(&f, PSIF_CC_OEI, irrep, 0, 0, lbl);
  global_dpd_->file2_mat_init(&f);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.occpi[Gq]; q++) {
        Q = moinfo.qt2pitzer[moinfo.qt_occ[q+moinfo.occ_off[Gq]]];
        f.matrix[Gp][p][q] = pertints[P][Q];
      }
    }
  }
  global_dpd_->file2_mat_wrt(&f);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);

  sprintf(lbl, "%s_AB", pert);
  global_dpd_->file2_init(&f, PSIF_CC_OEI, irrep, 1, 1, lbl);
  global_dpd_->file2_mat_init(&f);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.virtpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_vir[p+moinfo.vir_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
        Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
        f.matrix[Gp][p][q] = pertints[P][Q];
      }
    }
  }
  global_dpd_->file2_mat_wrt(&f);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);

  sprintf(lbl, "%s_IA", pert);
  global_dpd_->file2_init(&f, PSIF_CC_OEI, irrep, 0, 1, lbl);
  global_dpd_->file2_mat_init(&f);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
        Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
        f.matrix[Gp][p][q] = pertints[P][Q];
      }
    }
  }
  global_dpd_->file2_mat_wrt(&f);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);

}

}} // namespace psi::ccresponse
