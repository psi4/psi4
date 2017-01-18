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
#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::t1_build(void)
{
  dpdfile2 newtIA, newtia, tIA, tia, fIA, fia;
  dpdfile2 FAE, Fae, FMI, Fmi, FME, Fme;
  dpdfile2 dIA, dia;
  dpdbuf4 tIJAB, tijab, tIjAb, tiJaB, T2;
  dpdbuf4 C, C_anti, D, F_anti, F, E_anti, E, Z;
  int Gma, Gmi, Gm, Gi, Ga, ma, m, a, A, nrows, ncols;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "New tIA");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    global_dpd_->file2_close(&FAE);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    global_dpd_->file2_close(&FMI);

    global_dpd_->file2_close(&tIA);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");

    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    global_dpd_->contract422(&tIjAb, &FME, &newtIA, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tIjAb);

    global_dpd_->file2_close(&FME);

    global_dpd_->buf4_init(&C_anti, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
    global_dpd_->dot13(&tIA, &D, &newtIA, 0, 0, 1, 1);

    global_dpd_->file2_close(&tIA);

    global_dpd_->buf4_close(&C_anti);
    global_dpd_->buf4_close(&D);

    /*
      dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(ma,mi)");
      dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F 2<ia|bc> - <ia|cb>");
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_contract444(&F, &tIjAb, &Z, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_close(&F);
      dpd_trace42_13(&Z, &newtIA, 1, 1.0, 1.0);
      dpd_buf4_close(&Z);
    */

    /* t(i,a) <-- (2 t(mi,ef) - t(mi,fe)) <ma|ef> */
    /* out-of-core version replacing the *stupid* code above 3/22/05, TDC */
    global_dpd_->file2_mat_init(&newtIA);
    global_dpd_->file2_mat_rd(&newtIA);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    for(Gma=0; Gma < moinfo_.nirreps; Gma++) {
      Gmi = Gma; /* T1 is totally symmetric */

      global_dpd_->buf4_mat_irrep_row_init(&F, Gma);
      global_dpd_->buf4_mat_irrep_init(&T2, Gmi);
      global_dpd_->buf4_mat_irrep_rd(&T2, Gmi);

      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {

	global_dpd_->buf4_mat_irrep_row_rd(&F, Gma, ma);

	m = F.params->roworb[Gma][ma][0];
	a = F.params->roworb[Gma][ma][1];
	Gm = F.params->psym[m];
	Ga = F.params->qsym[a];
	Gi = Ga; /* T1 is totally symmetric */
	A = a - F.params->qoff[Ga];

    nrows = moinfo_.occpi[Gi];
	ncols = F.params->coltot[Gma];

    if(nrows && ncols && moinfo_.virtpi[Ga])
	  C_DGEMV('n',nrows,ncols,1.0,T2.matrix[Gmi][T2.row_offset[Gmi][m]],ncols,
          F.matrix[Gma][0],1,1.0,&newtIA.matrix[Gi][0][A],moinfo_.virtpi[Ga]);

      }

      global_dpd_->buf4_mat_irrep_close(&T2, Gmi);
      global_dpd_->buf4_mat_irrep_row_close(&F, Gma);
    }
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_mat_wrt(&newtIA);
    global_dpd_->file2_mat_close(&newtIA);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&E, &tIjAb, &newtIA, 1, 3, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&tIjAb);

    if (params_.just_residuals) {
      global_dpd_->file2_close(&newtIA);
      return;
    }

//    dpd_file2_copy(&newtIA, CC_OEI, "New tIA Increment");
    global_dpd_->file2_close(&newtIA);

/*
    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA Increment");
    if(params.local && local.filter_singles) {
      local_filter_T1(&newtIA);
    }
    else {
      dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &newtIA);
      dpd_file2_close(&dIA);
    }
    dpd_file2_close(&newtIA);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_copy(&tIA, CC_OEI, "New tIA");
    dpd_file2_close(&tIA);
    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "New tIA Increment");
    dpd_file2_axpy(&tIA, &newtIA, 1, 0);
    dpd_file2_close(&tIA);
    dpd_file2_close(&newtIA);
*/
  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "New tIA");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "New tia");
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&newtia, PSIF_CC_OEI, 0, 0, 1, "New tia");


    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 1, 1, "Fae");

    global_dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    global_dpd_->contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);

    global_dpd_->file2_close(&FAE);  global_dpd_->file2_close(&Fae);

    global_dpd_->file2_close(&tIA); global_dpd_->file2_close(&tia);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 0, 0, "Fmi");

    global_dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    global_dpd_->contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);

    global_dpd_->file2_close(&FMI);  global_dpd_->file2_close(&Fmi);
    global_dpd_->file2_close(&tIA);  global_dpd_->file2_close(&tia);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");

    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");

    global_dpd_->dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
    global_dpd_->dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);

    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);

    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    global_dpd_->dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
    global_dpd_->dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);

    global_dpd_->buf4_close(&tIjAb);

    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&Fme);

    global_dpd_->buf4_init(&C_anti, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0,
		  "C <ia||jb>");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
    global_dpd_->dot13(&tia, &D, &newtIA, 0, 0, 1, 1);

    global_dpd_->dot14(&tia, &C_anti, &newtia, 0, 1, -1, 1);
    global_dpd_->dot13(&tIA, &D, &newtia, 0, 0, 1, 1);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_close(&C_anti);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");

    global_dpd_->contract442(&tIJAB, &F_anti, &newtIA, 1, 1, 1, 1);
    global_dpd_->contract442(&tijab, &F_anti, &newtia, 1, 1, 1, 1);

    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&F_anti);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

    global_dpd_->contract442(&tiJaB, &F, &newtIA, 1, 1, 1, 1);
    global_dpd_->contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);

    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&tiJaB);

    global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");

    global_dpd_->contract442(&E_anti, &tIJAB, &newtIA, 1, 3, -1, 1);
    global_dpd_->contract442(&E_anti, &tijab, &newtia, 1, 3, -1, 1);

    global_dpd_->buf4_close(&E_anti);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

    global_dpd_->contract442(&E, &tiJaB, &newtIA, 1, 3, -1, 1);
    global_dpd_->contract442(&E, &tIjAb, &newtia, 1, 3, -1, 1);

    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&tiJaB);

/*
    dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newtIA);
    dpd_file2_close(&dIA);

    dpd_file2_init(&dia, CC_OEI, 0, 0, 1, "dia");
    dpd_file2_dirprd(&dia, &newtia);
    dpd_file2_close(&dia);
*/

    global_dpd_->file2_close(&newtIA);  global_dpd_->file2_close(&newtia);
  }
  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "New tIA");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "New tia");
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&newtia, PSIF_CC_OEI, 0, 2, 3, "New tia");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    global_dpd_->file2_close(&FAE);

    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");
    global_dpd_->contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);
    global_dpd_->file2_close(&Fae);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    global_dpd_->file2_close(&FMI);

    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");
    global_dpd_->contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);
    global_dpd_->file2_close(&Fmi);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");


    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tIJAB);

    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tijab);

    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
    global_dpd_->dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tIjAb);

    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&Fme);



    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->dot14(&tIA, &C, &newtIA, 0, 1, -1, 1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->dot14(&tia, &C, &newtia, 0, 1, -1, 1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->dot13(&tia, &D, &newtIA, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->dot13(&tIA, &D, &newtia, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);


    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&tIJAB, &F, &newtIA, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->contract442(&tijab, &F, &newtia, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&tIjAb, &F, &newtIA, 0, 2, 1, 1);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&F);



    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&E, &tIJAB, &newtIA, 1, 3, -1, 1);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->contract442(&E, &tijab, &newtia, 1, 3, -1, 1);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&E, &tIjAb, &newtIA, 2, 2, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&tIjAb);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&E, &tIjAb, &newtia, 2, 2, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&tIjAb);

/*
    dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newtIA);
    dpd_file2_close(&dIA);

    dpd_file2_init(&dia, CC_OEI, 0, 2, 3, "dia");
    dpd_file2_dirprd(&dia, &newtia);
    dpd_file2_close(&dia);
*/

    global_dpd_->file2_close(&newtIA);
    global_dpd_->file2_close(&newtia);

  }
}
}} // namespace psi::ccenergy
