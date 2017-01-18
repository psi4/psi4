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
#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::Fmi_build(void)
{
  int h,m,i;
  dpdfile2 FMI, Fmi, FMIt, Fmit, fIJ, fij, fIA, fia;
  dpdfile2 tIA, tia, FME, Fme;
  dpdbuf4 E_anti, E, D_anti, D;
  dpdbuf4 tautIJAB, tautijab, tautIjAb;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    global_dpd_->file2_close(&fIJ);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    global_dpd_->file2_close(&fIJ);

    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
    global_dpd_->file2_copy(&fij, PSIF_CC_OEI, "Fmi");
    global_dpd_->file2_close(&fij);
  }
  else if(params_.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    global_dpd_->file2_close(&fIJ);

    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
    global_dpd_->file2_copy(&fij, PSIF_CC_OEI, "Fmi");
    global_dpd_->file2_close(&fij);
  }

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_mat_init(&FMI);
    global_dpd_->file2_mat_rd(&FMI);

    /*
    for(h=0; h < moinfo.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++)
	FMI.matrix[h][m][m] = 0;
    }
    */

    global_dpd_->file2_mat_wrt(&FMI);
    global_dpd_->file2_mat_close(&FMI);
    global_dpd_->file2_close(&FMI);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 0, 0, "Fmi");

    global_dpd_->file2_mat_init(&FMI);
    global_dpd_->file2_mat_rd(&FMI);
    global_dpd_->file2_mat_init(&Fmi);
    global_dpd_->file2_mat_rd(&Fmi);

    for(h=0; h < moinfo_.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++)
	FMI.matrix[h][m][m] = 0;

      for(m=0; m < Fmi.params->rowtot[h]; m++)
	Fmi.matrix[h][m][m] = 0;
    }

    global_dpd_->file2_mat_wrt(&FMI);
    global_dpd_->file2_mat_close(&FMI);
    global_dpd_->file2_mat_wrt(&Fmi);
    global_dpd_->file2_mat_close(&Fmi);

    global_dpd_->file2_close(&FMI);
    global_dpd_->file2_close(&Fmi);
  }
  else if(params_.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");

    global_dpd_->file2_mat_init(&FMI);
    global_dpd_->file2_mat_rd(&FMI);
    global_dpd_->file2_mat_init(&Fmi);
    global_dpd_->file2_mat_rd(&Fmi);

    for(h=0; h < moinfo_.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++)
	FMI.matrix[h][m][m] = 0;

      for(m=0; m < Fmi.params->rowtot[h]; m++)
	Fmi.matrix[h][m][m] = 0;
    }

    global_dpd_->file2_mat_wrt(&FMI);
    global_dpd_->file2_mat_close(&FMI);
    global_dpd_->file2_mat_wrt(&Fmi);
    global_dpd_->file2_mat_close(&Fmi);

    global_dpd_->file2_close(&FMI);
    global_dpd_->file2_close(&Fmi);
  }

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
    global_dpd_->dot13(&tIA, &E, &FMI, 1, 1, 1.0, 1.0);

    global_dpd_->file2_close(&tIA);

    global_dpd_->buf4_close(&E_anti);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    global_dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautIjAb);
    global_dpd_->buf4_close(&D);

    /* Build the tilde intermediate */
    global_dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    global_dpd_->file2_close(&FMI);

    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&tIA);

    global_dpd_->file2_close(&FMIt);
  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 0, 0, "Fmi");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&fia);

    global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
    global_dpd_->dot13(&tia, &E, &FMI, 1, 1, 1.0, 1.0);

    global_dpd_->dot13(&tia, &E_anti, &Fmi, 1, 1, 1.0, 1.0);
    global_dpd_->dot13(&tIA, &E, &Fmi, 1, 1, 1.0, 1.0);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&E_anti);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    global_dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautijab");

    global_dpd_->contract442(&D_anti, &tautIJAB, &FMI, 0, 0, 1, 1);
    global_dpd_->contract442(&D_anti, &tautijab, &Fmi, 0, 0, 1, 1);

    global_dpd_->buf4_close(&tautIJAB);
    global_dpd_->buf4_close(&tautijab);
    global_dpd_->buf4_close(&D_anti);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    global_dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    global_dpd_->contract442(&D, &tautIjAb, &Fmi, 1, 1, 1, 1);

    global_dpd_->buf4_close(&tautIjAb);
    global_dpd_->buf4_close(&D);

    /* Build the tilde intermediate */
    global_dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    global_dpd_->file2_copy(&Fmi, PSIF_CC_OEI, "Fmit");

    global_dpd_->file2_close(&FMI);
    global_dpd_->file2_close(&Fmi);

    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 0, 0, "Fmit");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&tIA);

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
    global_dpd_->contract222(&Fme, &tia, &Fmit, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&Fme);
    global_dpd_->file2_close(&tia);

    global_dpd_->file2_close(&FMIt);
    global_dpd_->file2_close(&Fmit);
  }
  else if(params_.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");

    global_dpd_->dot13(&tIA, &E_anti, &FMI, 1, 1, 1, 1);
    global_dpd_->dot24(&tia, &E, &FMI, 0, 0, 1, 1);

    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&E_anti);

    global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");

    global_dpd_->dot13(&tia, &E_anti, &Fmi, 1, 1, 1, 1);
    global_dpd_->dot13(&tIA, &E, &Fmi, 1, 1, 1, 1);

    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&E_anti);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    global_dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    global_dpd_->contract442(&D, &tautIJAB, &FMI, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautIJAB);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tautijab");
    global_dpd_->contract442(&D, &tautijab, &Fmi, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautijab);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    global_dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautIjAb);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tautiJaB");
    global_dpd_->contract442(&D, &tautIjAb, &Fmi, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautIjAb);
    global_dpd_->buf4_close(&D);


    /* Build the tilde intermediate */
    global_dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    global_dpd_->file2_copy(&Fmi, PSIF_CC_OEI, "Fmit");

    global_dpd_->file2_close(&FMI);
    global_dpd_->file2_close(&Fmi);

    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 2, 2, "Fmit");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&tIA);

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->contract222(&Fme, &tia, &Fmit, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&Fme);
    global_dpd_->file2_close(&tia);

    global_dpd_->file2_close(&FMIt);
    global_dpd_->file2_close(&Fmit);
  }
}
}} // namespace psi::ccenergy
