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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void overlap(int L_irr)
{
  int h, nirreps;
  int row, col;
  int i,j,a,b,I,J,A,B,Isym,Jsym,Asym,Bsym;
  dpdfile2 T1, L1, T1A, T1B;
  dpdbuf4 T2, L2;
  double value = 1.0;
  double ST1A, ST1B, ST2AA, ST2BB, ST2AB, ST12AA, ST12BB, ST12AB;

  nirreps = moinfo.nirreps;

  global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  ST1A = global_dpd_->file2_dot(&T1, &L1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&T1);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  }
  ST1B = global_dpd_->file2_dot(&T1, &L1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&T1);

  global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  ST2AA = global_dpd_->buf4_dot(&L2, &T2);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&L2);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
  }
  ST2BB = global_dpd_->buf4_dot(&L2, &T2);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&L2);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  }
  ST2AB = global_dpd_->buf4_dot(&L2, &T2);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&L2);

  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1A);
  global_dpd_->file2_mat_rd(&T1A);
  if(params.ref == 0 || params.ref == 1) /** RHF/ROHF **/
    global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  else if(params.ref == 2) /** UHF **/
    global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&T1B);
  global_dpd_->file2_mat_rd(&T1B);

  ST12AA = 0.0;
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&L2, h);
    global_dpd_->buf4_mat_irrep_rd(&L2, h);
    for(row=0; row < L2.params->rowtot[h]; row++) {
      i = L2.params->roworb[h][row][0];
      j = L2.params->roworb[h][row][1];
      I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      for(col=0; col < L2.params->coltot[h^L_irr]; col++) {
	a = L2.params->colorb[h^L_irr][col][0];
	b = L2.params->colorb[h^L_irr][col][1];
	A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
	B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
	if((Isym == Asym) && (Jsym == Bsym))
	  ST12AA += L2.matrix[h][row][col] *
	    T1A.matrix[Isym][I][A] * T1A.matrix[Jsym][J][B];
	if((Isym == Bsym) && (Jsym == Asym))
	  ST12AA -= L2.matrix[h][row][col] *
	    T1A.matrix[Isym][I][B] * T1A.matrix[Jsym][J][A];
      }
    }
    global_dpd_->buf4_mat_irrep_close(&L2, h);
  }
  global_dpd_->buf4_close(&L2);

  ST12BB = 0.0;

  if(params.ref == 0 || params.ref == 1)
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
  else if(params.ref == 2)
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&L2, h);
    global_dpd_->buf4_mat_irrep_rd(&L2, h);
    for(row=0; row < L2.params->rowtot[h]; row++) {
      i = L2.params->roworb[h][row][0];
      j = L2.params->roworb[h][row][1];
      I = T1B.params->rowidx[i]; Isym = T1B.params->psym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      for(col=0; col < L2.params->coltot[h^L_irr]; col++) {
	a = L2.params->colorb[h^L_irr][col][0];
	b = L2.params->colorb[h^L_irr][col][1];
	A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
	B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
	if((Isym == Asym) && (Jsym == Bsym))
	  ST12BB += L2.matrix[h][row][col] *
	    T1B.matrix[Isym][I][A] * T1B.matrix[Jsym][J][B];
	if((Isym == Bsym) && (Jsym == Asym))
	  ST12BB -= L2.matrix[h][row][col] *
	    T1B.matrix[Isym][I][B] * T1B.matrix[Jsym][J][A];
      }
    }
    global_dpd_->buf4_mat_irrep_close(&L2, h);
  }
  global_dpd_->buf4_close(&L2);

  ST12AB = 0.0;

  if(params.ref == 0 || params.ref == 1)
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  else if(params.ref == 2)
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&L2, h);
    global_dpd_->buf4_mat_irrep_rd(&L2, h);
    for(row=0; row < L2.params->rowtot[h]; row++) {
      i = L2.params->roworb[h][row][0];
      j = L2.params->roworb[h][row][1];
      I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      for(col=0; col < L2.params->coltot[h^L_irr]; col++) {
	a = L2.params->colorb[h^L_irr][col][0];
	b = L2.params->colorb[h^L_irr][col][1];
	A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
	B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
	if((Isym == Asym) && (Jsym == Bsym))
	  ST12AB += L2.matrix[h][row][col] *
	    T1A.matrix[Isym][I][A] * T1B.matrix[Jsym][J][B];
      }
    }
    global_dpd_->buf4_mat_irrep_close(&L2, h);
  }
  global_dpd_->buf4_close(&L2);


  global_dpd_->file2_mat_close(&T1A);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_mat_close(&T1B);
  global_dpd_->file2_close(&T1B);

  /*
    outfile->Printf( "\tST1A = %20.15f\n", ST1A);
    outfile->Printf( "\tST1B = %20.15f\n", ST1B);
    outfile->Printf( "\tST2AA = %20.15f\n", ST2AA);
    outfile->Printf( "\tST2BB = %20.15f\n", ST2BB);
    outfile->Printf( "\tST2AB = %20.15f\n", ST2AB);
    outfile->Printf( "\tST12AA = %20.15f\n", ST12AA);
    outfile->Printf( "\tST12BB = %20.15f\n", ST12BB);
    outfile->Printf( "\tST12AB = %20.15f\n", ST12AB);
  */

  value = 1.0 - ST1A - ST1B - ST2AA - ST2BB - ST2AB + ST12AA + ST12BB + ST12AB;

  outfile->Printf( "\tOverlap <L|e^T> = %20.11f\n", value);
}

}} // namespace psi::cclambda
