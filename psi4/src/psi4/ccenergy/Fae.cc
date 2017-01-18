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
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "psi4/libciomr/libciomr.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::Fae_build(void)
{
  int h,a,e,nirreps;
  int ma,fe,ef,m,f,M,A,Gm,Ga,Ge,Gf,Gma,nrows,ncols;
  double *X;
  dpdfile2 tIA, tia;
  dpdfile2 FME, Fme;
  dpdfile2 fAB, fab, fIA, fia;
  dpdfile2 FAE, Fae;
  dpdfile2 FAEt, Faet;
  dpdbuf4 F_anti, F, D_anti, D;
  dpdbuf4 tautIJAB, tautijab, tautIjAb, taut;

  nirreps = moinfo_.nirreps;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    global_dpd_->file2_close(&fAB);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    global_dpd_->file2_close(&fAB);

    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
    global_dpd_->file2_copy(&fab, PSIF_CC_OEI, "Fae");
    global_dpd_->file2_close(&fab);
  }
  else if(params_.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    global_dpd_->file2_close(&fAB);

    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
    global_dpd_->file2_copy(&fab, PSIF_CC_OEI, "Fae");
    global_dpd_->file2_close(&fab);
  }

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");

    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);

    /*
      for(h=0; h < moinfo.nirreps; h++) {
      for(a=0; a < FAE.params->rowtot[h]; a++)
      FAE.matrix[h][a][a] = 0;
      }
    */

    global_dpd_->file2_mat_wrt(&FAE);
    global_dpd_->file2_mat_close(&FAE);
    global_dpd_->file2_close(&FAE);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 1, 1, "Fae");

    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);
    global_dpd_->file2_mat_init(&Fae);
    global_dpd_->file2_mat_rd(&Fae);

    for(h=0; h < moinfo_.nirreps; h++) {

      for(a=0; a < FAE.params->rowtot[h]; a++)
	for(e=0; e < FAE.params->coltot[h]; e++)
	  FAE.matrix[h][a][e] *= (1 - (a==e));

      for(a=0; a < Fae.params->rowtot[h]; a++)
	for(e=0; e < Fae.params->coltot[h]; e++)
	  Fae.matrix[h][a][e] *= (1 - (a==e));

    }

    global_dpd_->file2_mat_wrt(&FAE);
    global_dpd_->file2_mat_close(&FAE);
    global_dpd_->file2_mat_wrt(&Fae);
    global_dpd_->file2_mat_close(&Fae);

    global_dpd_->file2_close(&FAE);
    global_dpd_->file2_close(&Fae);
  }
  else if(params_.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");

    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);
    global_dpd_->file2_mat_init(&Fae);
    global_dpd_->file2_mat_rd(&Fae);

    for(h=0; h < moinfo_.nirreps; h++) {

      for(a=0; a < FAE.params->rowtot[h]; a++)
	for(e=0; e < FAE.params->coltot[h]; e++)
	  FAE.matrix[h][a][e] *= (1 - (a==e));

      for(a=0; a < Fae.params->rowtot[h]; a++)
	for(e=0; e < Fae.params->coltot[h]; e++)
	  Fae.matrix[h][a][e] *= (1 - (a==e));

    }

    global_dpd_->file2_mat_wrt(&FAE);
    global_dpd_->file2_mat_close(&FAE);
    global_dpd_->file2_mat_wrt(&Fae);
    global_dpd_->file2_mat_close(&Fae);

    global_dpd_->file2_close(&FAE);
    global_dpd_->file2_close(&Fae);
  }

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);
    global_dpd_->file2_close(&FAE);

    /*
      dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
      dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
      dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
      dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
      dpd_dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
      dpd_dot13(&tIA, &F, &FAE, 0, 0, 1.0, 1.0);
      dpd_file2_close(&tIA);
      dpd_buf4_close(&F_anti);
      dpd_buf4_close(&F);
      dpd_file2_close(&FAE);
    */

    /* Out-of-core algorithm for F->FAE added 3/20/05 - TDC */
    /* Fae <-- t(m,f) [2 <ma|fe> - <ma|ef>] */
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&tIA);
    global_dpd_->file2_mat_rd(&tIA);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    for(Gma=0; Gma < nirreps; Gma++) {
      global_dpd_->buf4_mat_irrep_row_init(&F, Gma);
      X = init_array(F.params->coltot[Gma]);

      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {
	global_dpd_->buf4_mat_irrep_row_rd(&F, Gma, ma);
	m = F.params->roworb[Gma][ma][0];
	a = F.params->roworb[Gma][ma][1];
	Gm = F.params->psym[m];
	Ga = Ge = Gm ^ Gma;  /* Fae is totally symmetric */
	Gf = Gm; /* T1 is totally symmetric */
	M = m - F.params->poff[Gm];
	A = a - F.params->qoff[Ga];

	zero_arr(X, F.params->coltot[Gma]);

	/* build spin-adapted F-integrals for current ma */
	for(fe=0; fe < F.params->coltot[Gma]; fe++) {
	  f = F.params->colorb[Gma][fe][0];
	  e = F.params->colorb[Gma][fe][1];
	  ef = F.params->colidx[e][f];
	  X[fe] = 2.0 * F.matrix[Gma][0][fe] - F.matrix[Gma][0][ef];
	}

    nrows = moinfo_.virtpi[Gf];
    ncols = moinfo_.virtpi[Ge];
	if(nrows && ncols)
	  C_DGEMV('t',nrows,ncols,1.0,&X[F.col_offset[Gma][Gf]],ncols,
		  tIA.matrix[Gm][M],1,1.0,
		  FAE.matrix[Ga][A],1);
      }

      free(X);
      global_dpd_->buf4_mat_irrep_row_close(&F, Gma);
    }
    global_dpd_->buf4_close(&F);
    global_dpd_->file2_mat_close(&tIA);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_mat_wrt(&FAE);
    global_dpd_->file2_mat_close(&FAE);
    global_dpd_->file2_close(&FAE);

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    global_dpd_->contract442(&tautIjAb, &D, &FAE, 3, 3, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tautIjAb);

    /* Build the tilde intermediates */
    global_dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    global_dpd_->file2_close(&FAE);

    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_close(&FAEt);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 1, 1, "Fae");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&fia);

    global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tia, &F, &FAE, 0, 0, 1.0, 1.0);

    global_dpd_->dot13(&tia, &F_anti, &Fae, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tIA, &F, &Fae, 0, 0, 1.0, 1.0);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&F_anti);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");

    global_dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    global_dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautijab");

    global_dpd_->contract442(&tautIJAB, &D_anti, &FAE, 2, 2, -1, 1);
    global_dpd_->contract442(&tautijab, &D_anti, &Fae, 2, 2, -1, 1);

    global_dpd_->buf4_close(&D_anti);
    global_dpd_->buf4_close(&tautIJAB);
    global_dpd_->buf4_close(&tautijab);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    global_dpd_->contract442(&tautIjAb, &D, &Fae, 3, 3, -1, 1);
    global_dpd_->contract442(&tautIjAb, &D, &FAE, 2, 2, -1, 1);

    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tautIjAb);


    /* Build the tilde intermediates */
    global_dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    global_dpd_->file2_copy(&Fae, PSIF_CC_OEI, "Faet");

    global_dpd_->file2_close(&FAE);
    global_dpd_->file2_close(&Fae);

    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 1, 1, "Faet");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
    global_dpd_->contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&Fme);

    global_dpd_->file2_close(&FAEt);
    global_dpd_->file2_close(&Faet);
  }
  else if(params_.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->dot13(&tIA, &F, &FAE, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->dot13(&tia, &F, &FAE, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->dot13(&tia, &F, &Fae, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->dot13(&tIA, &F, &Fae, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    global_dpd_->contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    global_dpd_->buf4_close(&taut);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    global_dpd_->contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    global_dpd_->buf4_close(&taut);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tautijab");
    global_dpd_->contract442(&taut, &D, &Fae, 2, 2, -1, 1);
    global_dpd_->buf4_close(&taut);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    global_dpd_->contract442(&taut, &D, &Fae, 3, 3, -1, 1);
    global_dpd_->buf4_close(&taut);
    global_dpd_->buf4_close(&D);

    /* Build the tilde intermediates */
    global_dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    global_dpd_->file2_copy(&Fae, PSIF_CC_OEI, "Faet");

    global_dpd_->file2_close(&FAE);
    global_dpd_->file2_close(&Fae);

    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 3, 3, "Faet");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&Fme);


    global_dpd_->file2_close(&FAEt);
    global_dpd_->file2_close(&Faet);
  }
}
}} // namespace psi::ccenergy
