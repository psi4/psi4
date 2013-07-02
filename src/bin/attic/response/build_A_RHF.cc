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
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

/* build_A_RHF_singlet(): Builds the RHF A matrix for RPA calculations.
** In spin orbitals, the A matrix is:
**
** A(ai,bj) = delta_ij f_ab - delta_ab f_ij - <ja||ib>
**
** RHF references and singlet eigenstates:
**  A(AI,BJ) = delta_IJ f_AB - delta_AB f_IJ + 2 <IJ|AB> - <IA|JB>
**
** RHF references and triplet eigenstates:
**  A(AI,BJ) = delta_IJ f_AB - delta_AB f_IJ - <IA|JB>
**
** TDC, March 2003
*/

void build_A_RHF(void)
{
  int h, nirreps;
  int a, b, i, j, ai, bj, A, B, I, J, Asym, Bsym, Isym, Jsym;
  dpdbuf4 C, D, Amat, A_AA, A_BB, A_AB;
  dpdfile2 fIJ, fij, fAB, fab;

  nirreps = moinfo.nirreps;

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_sort(&D, PSIF_MO_HESS, rpsq, 11, 11, "A(AI,BJ)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&Amat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  global_dpd_->buf4_scm(&Amat, 2);
  global_dpd_->buf4_close(&Amat);

  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  global_dpd_->buf4_sort_axpy(&C, PSIF_MO_HESS, qpsr, 11, 11, "A(AI,BJ)", -1);
  global_dpd_->buf4_sort(&C, PSIF_MO_HESS, qpsr, 11, 11, "A(AI,BJ) triplet");
  global_dpd_->buf4_close(&C);

  global_dpd_->buf4_init(&Amat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");
  global_dpd_->buf4_scm(&Amat, -1);
  global_dpd_->buf4_close(&Amat);

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_rd(&fIJ);
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_rd(&fij);
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_rd(&fAB);
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
  global_dpd_->file2_mat_init(&fab);
  global_dpd_->file2_mat_rd(&fab);

  global_dpd_->buf4_init(&Amat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Amat, h);
    global_dpd_->buf4_mat_irrep_rd(&Amat, h);
    for(ai=0; ai < Amat.params->rowtot[h]; ai++) {
      a = Amat.params->roworb[h][ai][0];
      i = Amat.params->roworb[h][ai][1];
      A = fAB.params->rowidx[a];
      I = fIJ.params->rowidx[i];
      Asym = fAB.params->psym[a];
      Isym = fIJ.params->psym[i];
      for(bj=0; bj < Amat.params->coltot[h]; bj++) {
	b = Amat.params->colorb[h][bj][0];
	j = Amat.params->colorb[h][bj][1];
	B = fAB.params->colidx[b];
	J = fIJ.params->colidx[j];
	Bsym = fAB.params->qsym[b];
	Jsym = fIJ.params->qsym[j];
	if((A==B) && (Isym==Jsym)) Amat.matrix[h][ai][bj] -= fIJ.matrix[Isym][I][J];
	if((I==J) && (Asym==Bsym)) Amat.matrix[h][ai][bj] += fAB.matrix[Asym][A][B];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Amat, h);
    global_dpd_->buf4_mat_irrep_close(&Amat, h);
  }
  global_dpd_->buf4_close(&Amat);

  global_dpd_->buf4_init(&Amat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Amat, h);
    global_dpd_->buf4_mat_irrep_rd(&Amat, h);
    for(ai=0; ai < Amat.params->rowtot[h]; ai++) {
      a = Amat.params->roworb[h][ai][0];
      i = Amat.params->roworb[h][ai][1];
      A = fAB.params->rowidx[a];
      I = fIJ.params->rowidx[i];
      Asym = fAB.params->psym[a];
      Isym = fIJ.params->psym[i];
      for(bj=0; bj < Amat.params->coltot[h]; bj++) {
	b = Amat.params->colorb[h][bj][0];
	j = Amat.params->colorb[h][bj][1];
	B = fAB.params->colidx[b];
	J = fIJ.params->colidx[j];
	Bsym = fAB.params->qsym[b];
	Jsym = fIJ.params->qsym[j];
	if((A==B) && (Isym==Jsym)) Amat.matrix[h][ai][bj] -= fIJ.matrix[Isym][I][J];
	if((I==J) && (Asym==Bsym)) Amat.matrix[h][ai][bj] += fAB.matrix[Asym][A][B];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Amat, h);
    global_dpd_->buf4_mat_irrep_close(&Amat, h);
  }
  global_dpd_->buf4_close(&Amat);

  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_close(&fab);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_close(&fIJ);
}

}} // namespace psi::response
