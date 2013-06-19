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
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

/* build_A(): Builds the CIS response matrix A, which is expressed in
** spin-orbitals as:
**
** A(ai,bj) = delta_ij f_ab - delta_ab f_ij - <ja||ib>
**
** This routine will build each spin component of the entire matrix 
** for later in-core diagonalization.
**
** RHF references:
**  A(AI,BJ) = delta_IJ f_AB - delta_AB f_IJ + 2 <IJ|AB> - <IA|JB>
**
** UHF and ROHF references:
**  A(AI,BJ) = delta_IJ f_AB - delta_AB f_IJ - <JA||IB>
**  A(ai,bj) = delta_ij f_ab - delta_ab f_ij - <ja||ib>
**  A(AI,bj) = <Ij|Ab>
**
** TDC, July 2002
*/

void build_A(void)
{
  int h, nirreps;
  int a, b, i, j, ai, bj, A, B, I, J, Asym, Bsym, Isym, Jsym;
  dpdbuf4 C, D, Amat, A_AA, A_BB, A_AB;
  dpdfile2 fIJ, fij, fAB, fab;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_MISC, rpsq, 11, 11, "A(AI,BJ)");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
    global_dpd_->buf4_scm(&Amat, 2);
    global_dpd_->buf4_close(&Amat);
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_sort_axpy(&C, PSIF_CC_MISC, qpsr, 11, 11, "A(AI,BJ)", -1);
    global_dpd_->buf4_sort(&C, PSIF_CC_MISC, qpsr, 11, 11, "A(AI,BJ) triplet");
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");
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

    global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
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

    global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");
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
  else if(params.ref == 2) { /** UHF **/

    /** <JA||IB> --> A(AI,BJ) **/
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->buf4_sort(&C, PSIF_CC_MISC, qrsp, 21, 21, "A(AI,BJ)");
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 21, 21, 21, 21, 0, "A(AI,BJ)");
    global_dpd_->buf4_scm(&Amat, -1);
    global_dpd_->buf4_close(&Amat);

    /** <ja||ib> --> A(ai,bj) **/
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_MISC, qrsp, 31, 31, "A(ai,bj)");
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 31, 31, 31, 31, 0, "A(ai,bj)");
    global_dpd_->buf4_scm(&Amat, -1);
    global_dpd_->buf4_close(&Amat);

    /** <Ij|Ab> --> A(AI,bj) **/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_MISC, rpsq, 21, 31, "A(AI,bj)");
    global_dpd_->buf4_close(&D);

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

    global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 21, 21, 21, 21, 0, "A(AI,BJ)");
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

    global_dpd_->buf4_init(&Amat, PSIF_CC_MISC, 0, 31, 31, 31, 31, 0, "A(ai,bj)");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&Amat, h);
      global_dpd_->buf4_mat_irrep_rd(&Amat, h);
      for(ai=0; ai < Amat.params->rowtot[h]; ai++) {
	a = Amat.params->roworb[h][ai][0];
	i = Amat.params->roworb[h][ai][1];
	A = fab.params->rowidx[a];
	I = fij.params->rowidx[i];
	Asym = fab.params->psym[a];
	Isym = fij.params->psym[i];
	for(bj=0; bj < Amat.params->coltot[h]; bj++) {
	  b = Amat.params->colorb[h][bj][0];
	  j = Amat.params->colorb[h][bj][1];
	  B = fab.params->colidx[b];
	  J = fij.params->colidx[j];
	  Bsym = fab.params->qsym[b];
	  Jsym = fij.params->qsym[j];
	  if((A==B) && (Isym==Jsym)) Amat.matrix[h][ai][bj] -= fij.matrix[Isym][I][J];
	  if((I==J) && (Asym==Bsym)) Amat.matrix[h][ai][bj] += fab.matrix[Asym][A][B];
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

}

}} // namespace psi::cis
