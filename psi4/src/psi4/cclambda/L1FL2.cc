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

/* L2FL2(): Computes the contributions of the Fme HBAR matrix elements
** to the Lambda doubles equations.  These contributions are given in
** spin orbitals as:
**
** L_ij^ab <-- P(ij) P(ab) L_i^a Fjb
**
** where Fjb = fjb + t_n^f <jn||bf>
**
** TDC, July 2002
*/

void L1FL2(int L_irr)
{
  int h, nirreps;
  int row,col;
  int i,j,a,b,I,J,A,B,Isym,Jsym,Asym,Bsym;
  dpdfile2 LIA, Lia, FJB, Fjb, L, F;
  dpdbuf4 newL2;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&L, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_mat_init(&L);
    global_dpd_->file2_mat_rd(&L);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);

    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&newL2, h);
      global_dpd_->buf4_mat_irrep_rd(&newL2, h);

      for(row=0; row < newL2.params->rowtot[h]; row++) {
	i = newL2.params->roworb[h][row][0];
	j = newL2.params->roworb[h][row][1];

	for(col=0; col < newL2.params->coltot[h^L_irr]; col++) {
	  a = newL2.params->colorb[h^L_irr][col][0];
	  b = newL2.params->colorb[h^L_irr][col][1];

	  I = L.params->rowidx[i]; Isym = L.params->psym[i];
	  J = F.params->rowidx[j]; Jsym = F.params->psym[j];
	  A = L.params->colidx[a]; Asym = L.params->qsym[a];
	  B = F.params->colidx[b]; Bsym = F.params->qsym[b];
	  if(((Isym^Asym) == L_irr) && (Jsym == Bsym))
	    newL2.matrix[h][row][col] += (L.matrix[Isym][I][A] * F.matrix[Jsym][J][B]);

	  if((Isym == Asym) && ((Jsym^Bsym) == L_irr))
	    newL2.matrix[h][row][col] += (L.matrix[Jsym][J][B] * F.matrix[Isym][I][A]);
	}
      }

      global_dpd_->buf4_mat_irrep_wrt(&newL2, h);
      global_dpd_->buf4_mat_irrep_close(&newL2, h);

    }

    global_dpd_->buf4_close(&newL2);

    global_dpd_->file2_mat_close(&F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_mat_close(&L);
    global_dpd_->file2_close(&L);

  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_mat_init(&LIA);
    global_dpd_->file2_mat_rd(&LIA);
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    global_dpd_->file2_mat_init(&Lia);
    global_dpd_->file2_mat_rd(&Lia);
    global_dpd_->file2_init(&FJB, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_mat_init(&FJB);
    global_dpd_->file2_mat_rd(&FJB);
    global_dpd_->file2_init(&Fjb, PSIF_CC_OEI, 0, 0, 1, "Fme");
    global_dpd_->file2_mat_init(&Fjb);
    global_dpd_->file2_mat_rd(&Fjb);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_mat_init(&LIA);
    global_dpd_->file2_mat_rd(&LIA);
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
    global_dpd_->file2_mat_init(&Lia);
    global_dpd_->file2_mat_rd(&Lia);
    global_dpd_->file2_init(&FJB, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_mat_init(&FJB);
    global_dpd_->file2_mat_rd(&FJB);
    global_dpd_->file2_init(&Fjb, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->file2_mat_init(&Fjb);
    global_dpd_->file2_mat_rd(&Fjb);

  }

  if(params.ref == 1) /** RHF/ROHF **/
    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
  else if(params.ref == 2) /** UHF **/
    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");

  if(params.ref == 1 || params.ref == 2) {
    /* loop over row irreps of LIJAB */
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&newL2, h);
      global_dpd_->buf4_mat_irrep_rd(&newL2, h);

      /* loop over rows of irrep of LIJAB */
      for(row=0; row < newL2.params->rowtot[h]; row++) {
	i = newL2.params->roworb[h][row][0];
	j = newL2.params->roworb[h][row][1];

	/* loop over cols of irrep of LIJAB */
	for(col=0; col < newL2.params->coltot[h^L_irr]; col++) {
	  a = newL2.params->colorb[h^L_irr][col][0];
	  b = newL2.params->colorb[h^L_irr][col][1];

	  I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
	  J = FJB.params->rowidx[j]; Jsym = FJB.params->psym[j];
	  A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
	  B = FJB.params->colidx[b]; Bsym = FJB.params->qsym[b];

	  if( ((Isym^Asym) == L_irr) && (Jsym == Bsym) )
	    newL2.matrix[h][row][col] += (LIA.matrix[Isym][I][A] *
					  FJB.matrix[Jsym][J][B]);

	  J = LIA.params->rowidx[j]; Jsym = LIA.params->psym[j];
	  I = FJB.params->rowidx[i]; Isym = FJB.params->psym[i];

	  if( (Isym == Asym) && ((Jsym^Bsym) == L_irr) )
	    newL2.matrix[h][row][col] += (LIA.matrix[Jsym][J][B] *
					  FJB.matrix[Isym][I][A]);

	  I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
	  J = FJB.params->rowidx[j]; Jsym = FJB.params->psym[j];
	  B = LIA.params->colidx[b]; Bsym = LIA.params->qsym[b];
	  A = FJB.params->colidx[a]; Asym = FJB.params->qsym[a];

	  if( ((Jsym^Asym) == L_irr) && (Isym == Bsym))
	    newL2.matrix[h][row][col] -= (LIA.matrix[Jsym][J][A] *
					  FJB.matrix[Isym][I][B]);

	  J = LIA.params->rowidx[j]; Jsym = LIA.params->psym[j];
	  I = FJB.params->rowidx[i]; Isym = FJB.params->psym[i];

	  if( (Jsym == Asym) && ((Isym^Bsym) == L_irr) )
	    newL2.matrix[h][row][col] -= (LIA.matrix[Isym][I][B] *
					  FJB.matrix[Jsym][J][A]);
	}
      }

      global_dpd_->buf4_mat_irrep_wrt(&newL2, h);
      global_dpd_->buf4_mat_irrep_close(&newL2, h);

    }
    global_dpd_->buf4_close(&newL2);
  }

  if(params.ref == 1) /** RHF/ROHF **/
    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
  else if(params.ref == 2) /** UHF **/
    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");

  if(params.ref == 1 || params.ref == 2) {
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&newL2, h);
      global_dpd_->buf4_mat_irrep_rd(&newL2, h);

      for(row=0; row < newL2.params->rowtot[h]; row++) {
	i = newL2.params->roworb[h][row][0];
	j = newL2.params->roworb[h][row][1];

	for(col=0; col < newL2.params->coltot[h^L_irr]; col++) {
	  a = newL2.params->colorb[h^L_irr][col][0];
	  b = newL2.params->colorb[h^L_irr][col][1];

	  I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
	  J = Fjb.params->rowidx[j]; Jsym = Fjb.params->psym[j];
	  A = Lia.params->colidx[a]; Asym = Lia.params->qsym[a];
	  B = Fjb.params->colidx[b]; Bsym = Fjb.params->qsym[b];

	  if(((Isym^Asym) == L_irr) && (Jsym == Bsym))
	    newL2.matrix[h][row][col] += (Lia.matrix[Isym][I][A] *
					  Fjb.matrix[Jsym][J][B]);

	  J = Lia.params->rowidx[j]; Jsym = Lia.params->psym[j];
	  I = Fjb.params->rowidx[i]; Isym = Fjb.params->psym[i];

	  if((Isym == Asym) && ((Jsym^Bsym) == L_irr))
	    newL2.matrix[h][row][col] += (Lia.matrix[Jsym][J][B] *
					  Fjb.matrix[Isym][I][A]);

	  I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
	  J = Fjb.params->rowidx[j]; Jsym = Fjb.params->psym[j];
	  B = Lia.params->colidx[b]; Bsym = Lia.params->qsym[b];
	  A = Fjb.params->colidx[a]; Asym = Fjb.params->qsym[a];

	  if(((Jsym^Asym) == L_irr) && (Isym == Bsym))
	    newL2.matrix[h][row][col] -= (Lia.matrix[Jsym][J][A] *
					  Fjb.matrix[Isym][I][B]);

	  J = Lia.params->rowidx[j]; Jsym = Lia.params->psym[j];
	  I = Fjb.params->rowidx[i]; Isym = Fjb.params->psym[i];

	  if((Jsym == Asym) && ((Isym^Bsym) == L_irr))
	    newL2.matrix[h][row][col] -= (Lia.matrix[Isym][I][B] *
					  Fjb.matrix[Jsym][J][A]);
	}
      }

      global_dpd_->buf4_mat_irrep_wrt(&newL2, h);
      global_dpd_->buf4_mat_irrep_close(&newL2, h);

    }
    global_dpd_->buf4_close(&newL2);
  }

  if(params.ref == 1) /** RHF/ROHF **/
    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
  else if(params.ref == 2) /** UHF **/
    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");

  if(params.ref == 1 || params.ref == 2) {
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&newL2, h);
      global_dpd_->buf4_mat_irrep_rd(&newL2, h);

      for(row=0; row < newL2.params->rowtot[h]; row++) {
	i = newL2.params->roworb[h][row][0];
	j = newL2.params->roworb[h][row][1];

	for(col=0; col < newL2.params->coltot[h^L_irr]; col++) {
	  a = newL2.params->colorb[h^L_irr][col][0];
	  b = newL2.params->colorb[h^L_irr][col][1];

	  I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
	  J = Fjb.params->rowidx[j]; Jsym = Fjb.params->psym[j];
	  A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
	  B = Fjb.params->colidx[b]; Bsym = Fjb.params->qsym[b];

	  if(((Isym^Asym) == L_irr) && (Jsym == Bsym))
	    newL2.matrix[h][row][col] += (LIA.matrix[Isym][I][A] *
					  Fjb.matrix[Jsym][J][B]);

	  J = Lia.params->rowidx[j]; Jsym = Lia.params->psym[j];
	  I = FJB.params->rowidx[i]; Isym = FJB.params->psym[i];
	  B = Lia.params->colidx[b]; Bsym = Lia.params->qsym[b];
	  A = FJB.params->colidx[a]; Asym = FJB.params->qsym[a];

	  if((Isym == Asym) && ((Jsym^Bsym) == L_irr))
	    newL2.matrix[h][row][col] += (Lia.matrix[Jsym][J][B] *
					  FJB.matrix[Isym][I][A]);
	}
      }

      global_dpd_->buf4_mat_irrep_wrt(&newL2, h);
      global_dpd_->buf4_mat_irrep_close(&newL2, h);

    }
  }

  if(params.ref == 1 || params.ref == 2) {
    global_dpd_->buf4_close(&newL2);

    global_dpd_->file2_mat_close(&FJB);
    global_dpd_->file2_close(&FJB);
    global_dpd_->file2_mat_close(&Fjb);
    global_dpd_->file2_close(&Fjb);
    global_dpd_->file2_mat_close(&LIA);
    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_mat_close(&Lia);
    global_dpd_->file2_close(&Lia);
  }
}

}} // namespace psi::cclambda
