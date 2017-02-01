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
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void denom1(dpdfile2 *X1, double omega);
void local_filter_T1(dpdfile2 *T1);

void X1_build(const char *pert, int irrep, double omega)
{
  dpdfile2 F, X1, X1new;
  dpdbuf4 W, X2;
  char lbl[32];
  int Gam, Gef, Gim, Gi, Ga, Gm, nrows, ncols, A, a, am;

  sprintf(lbl, "%sBAR_IA", pert);
  global_dpd_->file2_init(&X1new, PSIF_CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
  global_dpd_->file2_copy(&X1new, PSIF_CC_OEI, lbl);
  global_dpd_->file2_close(&X1new);
  global_dpd_->file2_init(&X1new, PSIF_CC_OEI, irrep, 0, 1, lbl);

  /*** S-S ***/

  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);

  global_dpd_->file2_axpy(&X1, &X1new, -omega, 0);

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "FAE");
  global_dpd_->contract222(&X1, &F, &X1new, 0, 0, 1, 1);
  global_dpd_->file2_close(&F);

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "FMI");
  global_dpd_->contract222(&F, &X1, &X1new, 1, 1, -1, 1);
  global_dpd_->file2_close(&F);

  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(jb,ME) + W(Jb,Me)");
  global_dpd_->contract422(&W, &X1, &X1new, 0, 0, 1, 1);
  global_dpd_->buf4_close(&W);

  global_dpd_->file2_close(&X1);


  /*** S-D ***/

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->dot24(&F, &X2, &X1new, 0, 0, 1, 1);
  global_dpd_->buf4_close(&X2);
  global_dpd_->file2_close(&F);

  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  /*  dpd_contract442(&X2, &W, &X1new, 0, 0, 1, 1); */
  /* ooc code below added 7/28/05, -TDC */
  global_dpd_->file2_mat_init(&X1new);
  global_dpd_->file2_mat_rd(&X1new);
  for(Gam=0; Gam < moinfo.nirreps; Gam++) {
    Gef = Gam; /* W is totally symmetric */
    Gim = Gef ^ irrep;

    global_dpd_->buf4_mat_irrep_init(&X2, Gim);
    global_dpd_->buf4_mat_irrep_rd(&X2, Gim);
    global_dpd_->buf4_mat_irrep_shift13(&X2, Gim);

    for(Gi=0; Gi < moinfo.nirreps; Gi++) {
      Ga = Gi ^ irrep;
      Gm = Ga ^ Gam;

      W.matrix[Gam] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gm], W.params->coltot[Gef]);

      nrows = moinfo.occpi[Gi];
      ncols = moinfo.occpi[Gm] * W.params->coltot[Gef];

      for(A=0; A < moinfo.virtpi[Ga]; A++) {
	a = moinfo.vir_off[Ga] + A;
	am = W.row_offset[Gam][a];

	global_dpd_->buf4_mat_irrep_rd_block(&W, Gam, am, moinfo.occpi[Gm]);

	if(nrows && ncols && moinfo.virtpi[Ga])
	  C_DGEMV('n',nrows,ncols,1,X2.shift.matrix[Gim][Gi][0],ncols,W.matrix[Gam][0],1,1,
		  &(X1new.matrix[Gi][0][A]), moinfo.virtpi[Ga]);
      }
      global_dpd_->free_dpd_block(W.matrix[Gam], moinfo.occpi[Gm], W.params->coltot[Gef]);
    }

    global_dpd_->buf4_mat_irrep_close(&X2, Gim);
  }
  global_dpd_->file2_mat_wrt(&X1new);
  global_dpd_->file2_mat_close(&X1new);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&X2);

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
  global_dpd_->contract442(&W, &X2, &X1new, 3, 3, 1, 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&X2);

  if(params.local && local.filter_singles) local_filter_T1(&X1new);
  else denom1(&X1new, omega);
  global_dpd_->file2_close(&X1new);
}

}} // namespace psi::ccresponse
