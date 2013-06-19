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
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void Fmi_build(void)
{
  int h,m,i;
  dpdfile2 FMI, Fmi, FMIt, Fmit, fIJ, fij, fIA, fia;
  dpdfile2 tIA, tia, FME, Fme;
  dpdbuf4 E_anti, E, D_anti, D;
  dpdbuf4 tautIJAB, tautijab, tautIjAb;

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    dpd_->file2_close(&fIJ);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    dpd_->file2_close(&fIJ);
  
    dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
    dpd_->file2_copy(&fij, PSIF_CC_OEI, "Fmi");
    dpd_->file2_close(&fij);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    dpd_->file2_close(&fIJ);
  
    dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
    dpd_->file2_copy(&fij, PSIF_CC_OEI, "Fmi");
    dpd_->file2_close(&fij);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    dpd_->file2_mat_init(&FMI);
    dpd_->file2_mat_rd(&FMI);

    /*
    for(h=0; h < moinfo.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++) 
	FMI.matrix[h][m][m] = 0;
    }
    */

    dpd_->file2_mat_wrt(&FMI);
    dpd_->file2_mat_close(&FMI);
    dpd_->file2_close(&FMI);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 0, 0, "Fmi");

    dpd_->file2_mat_init(&FMI);
    dpd_->file2_mat_rd(&FMI);
    dpd_->file2_mat_init(&Fmi);
    dpd_->file2_mat_rd(&Fmi);

    for(h=0; h < moinfo.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++) 
	FMI.matrix[h][m][m] = 0;

      for(m=0; m < Fmi.params->rowtot[h]; m++) 
	Fmi.matrix[h][m][m] = 0;
    }

    dpd_->file2_mat_wrt(&FMI);
    dpd_->file2_mat_close(&FMI);
    dpd_->file2_mat_wrt(&Fmi);
    dpd_->file2_mat_close(&Fmi);

    dpd_->file2_close(&FMI);
    dpd_->file2_close(&Fmi);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");

    dpd_->file2_mat_init(&FMI);
    dpd_->file2_mat_rd(&FMI);
    dpd_->file2_mat_init(&Fmi);
    dpd_->file2_mat_rd(&Fmi);

    for(h=0; h < moinfo.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++) 
	FMI.matrix[h][m][m] = 0;

      for(m=0; m < Fmi.params->rowtot[h]; m++) 
	Fmi.matrix[h][m][m] = 0;
    }

    dpd_->file2_mat_wrt(&FMI);
    dpd_->file2_mat_close(&FMI);
    dpd_->file2_mat_wrt(&Fmi);
    dpd_->file2_mat_close(&Fmi);

    dpd_->file2_close(&FMI);
    dpd_->file2_close(&Fmi);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&fIA);
  
    dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    dpd_->dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
    dpd_->dot13(&tIA, &E, &FMI, 1, 1, 1.0, 1.0);

    dpd_->file2_close(&tIA);

    dpd_->buf4_close(&E_anti);
    dpd_->buf4_close(&E);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_->buf4_close(&tautIjAb);
    dpd_->buf4_close(&D);

    /* Build the tilde intermediate */
    dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    dpd_->file2_close(&FMI);

    dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    dpd_->file2_close(&FME);
    dpd_->file2_close(&tIA);

    dpd_->file2_close(&FMIt);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 0, 0, "Fmi");

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&fIA);
  
    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    dpd_->file2_close(&tia);
    dpd_->file2_close(&fia);
  
    dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    dpd_->dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
    dpd_->dot13(&tia, &E, &FMI, 1, 1, 1.0, 1.0);

    dpd_->dot13(&tia, &E_anti, &Fmi, 1, 1, 1.0, 1.0);
    dpd_->dot13(&tIA, &E, &Fmi, 1, 1, 1.0, 1.0);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);
    dpd_->buf4_close(&E_anti);
    dpd_->buf4_close(&E);

    dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautijab");

    dpd_->contract442(&D_anti, &tautIJAB, &FMI, 0, 0, 1, 1);
    dpd_->contract442(&D_anti, &tautijab, &Fmi, 0, 0, 1, 1);

    dpd_->buf4_close(&tautIJAB);
    dpd_->buf4_close(&tautijab);
    dpd_->buf4_close(&D_anti);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_->contract442(&D, &tautIjAb, &Fmi, 1, 1, 1, 1);

    dpd_->buf4_close(&tautIjAb);
    dpd_->buf4_close(&D);

    /* Build the tilde intermediate */
    dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    dpd_->file2_copy(&Fmi, PSIF_CC_OEI, "Fmit");

    dpd_->file2_close(&FMI);
    dpd_->file2_close(&Fmi);

    dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 0, 0, "Fmit");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    dpd_->file2_close(&FME);
    dpd_->file2_close(&tIA);

    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
    dpd_->contract222(&Fme, &tia, &Fmit, 0, 0, 0.5, 1);
    dpd_->file2_close(&Fme);
    dpd_->file2_close(&tia);

    dpd_->file2_close(&FMIt);
    dpd_->file2_close(&Fmit);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&fIA);
  
    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    dpd_->file2_close(&tia);
    dpd_->file2_close(&fia);
  
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");

    dpd_->dot13(&tIA, &E_anti, &FMI, 1, 1, 1, 1);
    dpd_->dot24(&tia, &E, &FMI, 0, 0, 1, 1);

    dpd_->buf4_close(&E);
    dpd_->buf4_close(&E_anti);

    dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");

    dpd_->dot13(&tia, &E_anti, &Fmi, 1, 1, 1, 1);
    dpd_->dot13(&tIA, &E, &Fmi, 1, 1, 1, 1);

    dpd_->buf4_close(&E);
    dpd_->buf4_close(&E_anti);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    dpd_->contract442(&D, &tautIJAB, &FMI, 0, 0, 1, 1);
    dpd_->buf4_close(&tautIJAB);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tautijab");
    dpd_->contract442(&D, &tautijab, &Fmi, 0, 0, 1, 1);
    dpd_->buf4_close(&tautijab);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_->buf4_close(&tautIjAb);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tautiJaB");
    dpd_->contract442(&D, &tautIjAb, &Fmi, 0, 0, 1, 1);
    dpd_->buf4_close(&tautIjAb);
    dpd_->buf4_close(&D);


    /* Build the tilde intermediate */
    dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    dpd_->file2_copy(&Fmi, PSIF_CC_OEI, "Fmit");

    dpd_->file2_close(&FMI);
    dpd_->file2_close(&Fmi);

    dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 2, 2, "Fmit");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    dpd_->file2_close(&FME);
    dpd_->file2_close(&tIA);

    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    dpd_->contract222(&Fme, &tia, &Fmit, 0, 0, 0.5, 1);
    dpd_->file2_close(&Fme);
    dpd_->file2_close(&tia);

    dpd_->file2_close(&FMIt);
    dpd_->file2_close(&Fmit);
  }
}
}} // namespace psi::ccenergy
