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
    \ingroup STABLE
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace stable {

/* build_A_RHF_singlet(): Builds the RHF orbital Hessian for computing
** spatial- and spin-symmetry instabilities.  In spin orbitals the Hessian
** matrix is:
**
** A(ai,bj) = delta_ij f_ab - delta_ab f_ij + <ij||ab> - <ja||ib>
**
** RHF references and singlet eigenstates:
**  A(AI,BJ) = delta_IJ f_AB - delta_AB f_IJ + 4 <IJ|AB> - <IJ|BA> - <IA|JB>
**
** RHF references and triplet eigenstates:
**  A(AI,BJ) = delta_IJ f_AB - delta_AB f_IJ - <IJ|BA> - <IA|JB>
**
** This routine will build each spin component of the entire matrix 
** for later in-core diagonalization.
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

  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_->buf4_sort(&D, PSIF_MO_HESS, rpsq, 11, 11, "A(AI,BJ)");
  dpd_->buf4_close(&D);
  dpd_->buf4_init(&Amat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_->buf4_scm(&Amat, 4);
  dpd_->buf4_close(&Amat);

  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_->buf4_sort_axpy(&D, PSIF_MO_HESS, sprq, 11, 11, "A(AI,BJ)", -1);
  dpd_->buf4_sort(&D, PSIF_MO_HESS, sprq, 11, 11, "A(AI,BJ) triplet");
  dpd_->buf4_close(&D);

  dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_->buf4_sort_axpy(&C, PSIF_MO_HESS, qpsr, 11, 11, "A(AI,BJ)", -1);
  dpd_->buf4_sort_axpy(&C, PSIF_MO_HESS, qpsr, 11, 11, "A(AI,BJ) triplet", 1);
  dpd_->buf4_close(&C);

  dpd_->buf4_init(&Amat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");
  dpd_->buf4_scm(&Amat, -1);
  dpd_->buf4_close(&Amat);

  dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  dpd_->file2_mat_init(&fIJ);
  dpd_->file2_mat_rd(&fIJ);
  dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
  dpd_->file2_mat_init(&fij);
  dpd_->file2_mat_rd(&fij);
  dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  dpd_->file2_mat_init(&fAB);
  dpd_->file2_mat_rd(&fAB);
  dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
  dpd_->file2_mat_init(&fab);
  dpd_->file2_mat_rd(&fab);

  dpd_->buf4_init(&Amat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  for(h=0; h < nirreps; h++) {
    dpd_->buf4_mat_irrep_init(&Amat, h);
    dpd_->buf4_mat_irrep_rd(&Amat, h);
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
    dpd_->buf4_mat_irrep_wrt(&Amat, h);
    dpd_->buf4_mat_irrep_close(&Amat, h);
  }
  dpd_->buf4_close(&Amat);

  dpd_->buf4_init(&Amat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");
  for(h=0; h < nirreps; h++) {
    dpd_->buf4_mat_irrep_init(&Amat, h);
    dpd_->buf4_mat_irrep_rd(&Amat, h);
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
    dpd_->buf4_mat_irrep_wrt(&Amat, h);
    dpd_->buf4_mat_irrep_close(&Amat, h);
  }
  dpd_->buf4_close(&Amat);

  dpd_->file2_mat_close(&fab);
  dpd_->file2_close(&fab);
  dpd_->file2_mat_close(&fAB);
  dpd_->file2_close(&fAB);
  dpd_->file2_mat_close(&fij);
  dpd_->file2_close(&fij);
  dpd_->file2_mat_close(&fIJ);
  dpd_->file2_close(&fIJ);
}

}} // namespace psi::stable
