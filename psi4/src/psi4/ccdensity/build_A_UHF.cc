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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* BUILD_A_UHF(): Construct the components of the molecular orbital
** Hessian, A, for UHF orbitals. At the moment we're actually building
** all symmetry blocks of A, though for the orbital Z-vector equations
** we really only need the totally symmetric components.
**
** A(AI,JB) = delta_IJ f_AB - delta_AB f_IJ + <IJ||AB> - <JA||IB>
**
** A(ai,jb) = delta_ij f_ab - delta_ab f_ij + <ij||ab> - <ja||ib>
**
** A(AI,jb) = 2<Ij|Ab>
**
** TDC, January 2003
*/

void build_A_UHF(void)
{
  int h, nirreps, row, col;
  int a, i, b, j;
  int A, I, B, J;
  int Asym, Isym, Bsym, Jsym;
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia;
  dpdbuf4 Amat, D, C;

  nirreps = moinfo.nirreps;

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_rd(&fIJ);
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_rd(&fij);
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_rd(&fAB);
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
  global_dpd_->file2_mat_init(&fab);
  global_dpd_->file2_mat_rd(&fab);
  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_mat_init(&fIA);
  global_dpd_->file2_mat_rd(&fIA);
  global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_mat_init(&fia);
  global_dpd_->file2_mat_rd(&fia);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 1, "D <IJ|AB>");
  global_dpd_->buf4_sort(&D, PSIF_CC_MISC, rpsq, 21, 21, "A(AI,BJ)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
  global_dpd_->buf4_sort_axpy(&D, PSIF_CC_MISC, qrsp, 21, 21, "A(AI,BJ)", -1.0);
  global_dpd_->buf4_close(&D);

  global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 21, 21, 21, 21, 0, "A(AI,BJ)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Amat, h);
    global_dpd_->buf4_mat_irrep_rd(&Amat, h);

    for(row=0; row < Amat.params->rowtot[h]; row++) {
      a = Amat.params->roworb[h][row][0];
      i = Amat.params->roworb[h][row][1];
      A = fAB.params->rowidx[a];
      I = fIJ.params->rowidx[i];
      Asym = fAB.params->psym[a];
      Isym = fIJ.params->psym[i];
      for(col=0; col < Amat.params->coltot[h]; col++) {
	b = Amat.params->colorb[h][col][0];
	j = Amat.params->colorb[h][col][1];
	B = fAB.params->colidx[b];
	J = fIJ.params->colidx[j];
	Bsym = fAB.params->qsym[b];
	Jsym = fIJ.params->qsym[j];

	if((I==J) && (Asym==Bsym))
	  Amat.matrix[h][row][col] += fAB.matrix[Asym][A][B];
	if((A==B) && (Isym==Jsym))
	  Amat.matrix[h][row][col] -= fIJ.matrix[Isym][I][J];

      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&Amat, h);
    global_dpd_->buf4_mat_irrep_close(&Amat, h);
  }
  global_dpd_->buf4_close(&Amat);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 1, "D <ij|ab>");
  global_dpd_->buf4_sort(&D, PSIF_CC_MISC, rpsq, 31, 31, "A(ai,bj)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
  global_dpd_->buf4_sort_axpy(&D, PSIF_CC_MISC, qrsp, 31, 31, "A(ai,bj)", -1.0);
  global_dpd_->buf4_close(&D);

  global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 31, 31, 31, 31, 0, "A(ai,bj)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Amat, h);
    global_dpd_->buf4_mat_irrep_rd(&Amat, h);

    for(row=0; row < Amat.params->rowtot[h]; row++) {
      a = Amat.params->roworb[h][row][0];
      i = Amat.params->roworb[h][row][1];
      A = fab.params->rowidx[a];
      I = fij.params->rowidx[i];
      Asym = fab.params->psym[a];
      Isym = fij.params->psym[i];
      for(col=0; col < Amat.params->coltot[h]; col++) {
	b = Amat.params->colorb[h][col][0];
	j = Amat.params->colorb[h][col][1];
	B = fab.params->colidx[b];
	J = fij.params->colidx[j];
	Bsym = fab.params->qsym[b];
	Jsym = fij.params->qsym[j];

	if((I==J) && (Asym==Bsym))
	  Amat.matrix[h][row][col] += fab.matrix[Asym][A][B];
	if((A==B) && (Isym==Jsym))
	  Amat.matrix[h][row][col] -= fij.matrix[Isym][I][J];

      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&Amat, h);
    global_dpd_->buf4_mat_irrep_close(&Amat, h);
  }
  global_dpd_->buf4_close(&Amat);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  global_dpd_->buf4_sort(&D, PSIF_CC_MISC, rpsq, 21, 31, "A(AI,bj)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 21, 31, 21, 31, 0, "A(AI,bj)");
  global_dpd_->buf4_scm(&Amat, 2);
  global_dpd_->buf4_close(&Amat);

  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_close(&fab);
  global_dpd_->file2_mat_close(&fIA);
  global_dpd_->file2_close(&fIA);
  global_dpd_->file2_mat_close(&fia);
  global_dpd_->file2_close(&fia);
}


}} // namespace psi::ccdensity
