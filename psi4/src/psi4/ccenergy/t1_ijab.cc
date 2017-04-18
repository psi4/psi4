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
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <psi4/libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::t1_ijab(void)
{
  int h, ij, ab, i, j, a, b, I, J, A, B;
  int Isym, Jsym, Asym, Bsym;
  int nirreps;
  dpdbuf4 t1_IJAB, t1_ijab, t1_IjAb, t1_IjbA; 
  dpdfile2 tIA, tia;

  nirreps = moinfo_.nirreps;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&tIA);
    global_dpd_->file2_mat_rd(&tIA);

    global_dpd_->buf4_init(&t1_IjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_IjAb");

    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&t1_IjAb, h);

      for(ij=0; ij < t1_IjAb.params->rowtot[h]; ij++) {
	i = t1_IjAb.params->roworb[h][ij][0];
	j = t1_IjAb.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tIA.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tIA.params->psym[j];
	for(ab=0; ab < t1_IjAb.params->coltot[h]; ab++) {
	  a = t1_IjAb.params->colorb[h][ab][0];
	  b = t1_IjAb.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tIA.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tIA.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    t1_IjAb.matrix[h][ij][ab] =
	      (tIA.matrix[Isym][I][A] * tIA.matrix[Jsym][J][B]);

	}
      }

      global_dpd_->buf4_mat_irrep_wrt(&t1_IjAb, h);
      global_dpd_->buf4_mat_irrep_close(&t1_IjAb, h);
    }
    //global_dpd_->buf4_sort(&t1_IjAb, PSIF_CC_TAMPS, pqsr, 0, 5, "t1_IjbA");
    global_dpd_->buf4_close(&t1_IjAb);

    global_dpd_->file2_mat_close(&tIA);
    global_dpd_->file2_close(&tIA);
  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&tIA);
    global_dpd_->file2_mat_rd(&tIA);
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_mat_init(&tia);
    global_dpd_->file2_mat_rd(&tia);

    global_dpd_->buf4_init(&t1_IJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_IJAB");

    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&t1_IJAB, h);

      for(ij=0; ij < t1_IJAB.params->rowtot[h]; ij++) {
	i = t1_IJAB.params->roworb[h][ij][0];
	j = t1_IJAB.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tIA.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tIA.params->psym[j];
	for(ab=0; ab < t1_IJAB.params->coltot[h]; ab++) {
	  a = t1_IJAB.params->colorb[h][ab][0];
	  b = t1_IJAB.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tIA.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tIA.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    t1_IJAB.matrix[h][ij][ab] =
	      (tIA.matrix[Isym][I][A] * tIA.matrix[Jsym][J][B]);
	  if((Isym==Bsym) && (Jsym==Asym))
	    t1_IJAB.matrix[h][ij][ab] -=
	      (tIA.matrix[Isym][I][B] * tIA.matrix[Jsym][J][A]);

	}
      }

      global_dpd_->buf4_mat_irrep_wrt(&t1_IJAB, h);
      global_dpd_->buf4_mat_irrep_close(&t1_IJAB, h);
    }

    global_dpd_->buf4_close(&t1_IJAB);

    global_dpd_->buf4_init(&t1_ijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_ijab");

    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&t1_ijab, h);

      for(ij=0; ij < t1_ijab.params->rowtot[h]; ij++) {
	i = t1_ijab.params->roworb[h][ij][0];
	j = t1_ijab.params->roworb[h][ij][1];
	I = tia.params->rowidx[i];
	J = tia.params->rowidx[j];
	Isym = tia.params->psym[i];
	Jsym = tia.params->psym[j];
	for(ab=0; ab < t1_ijab.params->coltot[h]; ab++) {
	  a = t1_ijab.params->colorb[h][ab][0];
	  b = t1_ijab.params->colorb[h][ab][1];
	  A = tia.params->colidx[a];
	  B = tia.params->colidx[b];
	  Asym = tia.params->qsym[a];
	  Bsym = tia.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    t1_ijab.matrix[h][ij][ab] =
	      (tia.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);
	  if((Isym==Bsym) && (Jsym==Asym))
	    t1_ijab.matrix[h][ij][ab] -=
	      (tia.matrix[Isym][I][B] * tia.matrix[Jsym][J][A]);

	}
      }

      global_dpd_->buf4_mat_irrep_wrt(&t1_ijab, h);
      global_dpd_->buf4_mat_irrep_close(&t1_ijab, h);
    }

    global_dpd_->buf4_close(&t1_ijab);

    global_dpd_->buf4_init(&t1_IjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_IjAb");

    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&t1_IjAb, h);

      for(ij=0; ij < t1_IjAb.params->rowtot[h]; ij++) {
	i = t1_IjAb.params->roworb[h][ij][0];
	j = t1_IjAb.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tia.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tia.params->psym[j];
	for(ab=0; ab < t1_IjAb.params->coltot[h]; ab++) {
	  a = t1_IjAb.params->colorb[h][ab][0];
	  b = t1_IjAb.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tia.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tia.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    t1_IjAb.matrix[h][ij][ab] =
	      (tIA.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);

	}
      }

      global_dpd_->buf4_mat_irrep_wrt(&t1_IjAb, h);
      global_dpd_->buf4_mat_irrep_close(&t1_IjAb, h);
    }

    global_dpd_->buf4_sort(&t1_IjAb, PSIF_CC_TAMPS, pqsr, 0, 5, "t1_IjbA");
    global_dpd_->buf4_close(&t1_IjAb);
    global_dpd_->buf4_init(&t1_IjbA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_IjbA");
    global_dpd_->buf4_sort(&t1_IjbA, PSIF_CC_TAMPS, qprs,  0, 5, "t1_iJaB");
    global_dpd_->buf4_close(&t1_IjbA);

    global_dpd_->file2_mat_close(&tIA);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_mat_close(&tia);
    global_dpd_->file2_close(&tia);

  }
  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&tIA);
    global_dpd_->file2_mat_rd(&tIA);
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_mat_init(&tia);
    global_dpd_->file2_mat_rd(&tia);

    global_dpd_->buf4_init(&t1_IJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_IJAB");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&t1_IJAB, h);
      for(ij=0; ij < t1_IJAB.params->rowtot[h]; ij++) {
	i = t1_IJAB.params->roworb[h][ij][0];
	j = t1_IJAB.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tIA.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tIA.params->psym[j];
	for(ab=0; ab < t1_IJAB.params->coltot[h]; ab++) {
	  a = t1_IJAB.params->colorb[h][ab][0];
	  b = t1_IJAB.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tIA.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tIA.params->qsym[b];
	  if((Isym==Asym) && (Jsym==Bsym))
	    t1_IJAB.matrix[h][ij][ab] =
	      (tIA.matrix[Isym][I][A] * tIA.matrix[Jsym][J][B]);
	  if((Isym==Bsym) && (Jsym==Asym))
	    t1_IJAB.matrix[h][ij][ab] -=
	      (tIA.matrix[Isym][I][B] * tIA.matrix[Jsym][J][A]);
	}
      }
      global_dpd_->buf4_mat_irrep_wrt(&t1_IJAB, h);
      global_dpd_->buf4_mat_irrep_close(&t1_IJAB, h);
    }
    global_dpd_->buf4_close(&t1_IJAB);

    global_dpd_->buf4_init(&t1_ijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "t1_ijab");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&t1_ijab, h);
      for(ij=0; ij < t1_ijab.params->rowtot[h]; ij++) {
	i = t1_ijab.params->roworb[h][ij][0];
	j = t1_ijab.params->roworb[h][ij][1];
	I = tia.params->rowidx[i];
	J = tia.params->rowidx[j];
	Isym = tia.params->psym[i];
	Jsym = tia.params->psym[j];
	for(ab=0; ab < t1_ijab.params->coltot[h]; ab++) {
	  a = t1_ijab.params->colorb[h][ab][0];
	  b = t1_ijab.params->colorb[h][ab][1];
	  A = tia.params->colidx[a];
	  B = tia.params->colidx[b];
	  Asym = tia.params->qsym[a];
	  Bsym = tia.params->qsym[b];
	  if((Isym==Asym) && (Jsym==Bsym))
	    t1_ijab.matrix[h][ij][ab] =
	      (tia.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);
	  if((Isym==Bsym) && (Jsym==Asym))
	    t1_ijab.matrix[h][ij][ab] -=
	      (tia.matrix[Isym][I][B] * tia.matrix[Jsym][J][A]);
	}
      }
      global_dpd_->buf4_mat_irrep_wrt(&t1_ijab, h);
      global_dpd_->buf4_mat_irrep_close(&t1_ijab, h);
    }
    global_dpd_->buf4_close(&t1_ijab);

    global_dpd_->buf4_init(&t1_IjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "t1_IjAb");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&t1_IjAb, h);
      for(ij=0; ij < t1_IjAb.params->rowtot[h]; ij++) {
	i = t1_IjAb.params->roworb[h][ij][0];
	j = t1_IjAb.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tia.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tia.params->psym[j];
	for(ab=0; ab < t1_IjAb.params->coltot[h]; ab++) {
	  a = t1_IjAb.params->colorb[h][ab][0];
	  b = t1_IjAb.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tia.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tia.params->qsym[b];
	  if((Isym==Asym) && (Jsym==Bsym))
	    t1_IjAb.matrix[h][ij][ab] =
	      (tIA.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);
	}
      }
      global_dpd_->buf4_mat_irrep_wrt(&t1_IjAb, h);
      global_dpd_->buf4_mat_irrep_close(&t1_IjAb, h);
    }
    global_dpd_->buf4_close(&t1_IjAb);

    global_dpd_->file2_mat_close(&tIA);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_mat_close(&tia);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_init(&t1_IjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "t1_IjAb");
    global_dpd_->buf4_sort(&t1_IjAb, PSIF_CC_TAMPS, pqsr, 22, 29, "t1_IjbA");
    global_dpd_->buf4_close(&t1_IjAb);
    global_dpd_->buf4_init(&t1_IjbA, PSIF_CC_TAMPS, 0, 22, 29, 22, 29, 0, "t1_IjbA");
    global_dpd_->buf4_sort(&t1_IjbA, PSIF_CC_TAMPS, qprs,  23, 29, "t1_iJaB");
    global_dpd_->buf4_close(&t1_IjbA);



  } /*** UHF ***/
}
}} // namespace psi::ccenergy
