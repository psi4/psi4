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
    dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_copy(&fIJ, CC_OEI, "FMI");
    dpd_file2_close(&fIJ);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_copy(&fIJ, CC_OEI, "FMI");
    dpd_file2_close(&fIJ);
  
    dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
    dpd_file2_copy(&fij, CC_OEI, "Fmi");
    dpd_file2_close(&fij);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_copy(&fIJ, CC_OEI, "FMI");
    dpd_file2_close(&fIJ);
  
    dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
    dpd_file2_copy(&fij, CC_OEI, "Fmi");
    dpd_file2_close(&fij);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);

    /*
    for(h=0; h < moinfo.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++) 
	FMI.matrix[h][m][m] = 0;
    }
    */

    dpd_file2_mat_wrt(&FMI);
    dpd_file2_mat_close(&FMI);
    dpd_file2_close(&FMI);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&Fmi, CC_OEI, 0, 0, 0, "Fmi");

    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);
    dpd_file2_mat_init(&Fmi);
    dpd_file2_mat_rd(&Fmi);

    for(h=0; h < moinfo.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++) 
	FMI.matrix[h][m][m] = 0;

      for(m=0; m < Fmi.params->rowtot[h]; m++) 
	Fmi.matrix[h][m][m] = 0;
    }

    dpd_file2_mat_wrt(&FMI);
    dpd_file2_mat_close(&FMI);
    dpd_file2_mat_wrt(&Fmi);
    dpd_file2_mat_close(&Fmi);

    dpd_file2_close(&FMI);
    dpd_file2_close(&Fmi);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&Fmi, CC_OEI, 0, 2, 2, "Fmi");

    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);
    dpd_file2_mat_init(&Fmi);
    dpd_file2_mat_rd(&Fmi);

    for(h=0; h < moinfo.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++) 
	FMI.matrix[h][m][m] = 0;

      for(m=0; m < Fmi.params->rowtot[h]; m++) 
	Fmi.matrix[h][m][m] = 0;
    }

    dpd_file2_mat_wrt(&FMI);
    dpd_file2_mat_close(&FMI);
    dpd_file2_mat_wrt(&Fmi);
    dpd_file2_mat_close(&Fmi);

    dpd_file2_close(&FMI);
    dpd_file2_close(&Fmi);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);
  
    dpd_buf4_init(&E_anti, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
    dpd_dot13(&tIA, &E, &FMI, 1, 1, 1.0, 1.0);

    dpd_file2_close(&tIA);

    dpd_buf4_close(&E_anti);
    dpd_buf4_close(&E);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    dpd_contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_buf4_close(&tautIjAb);
    dpd_buf4_close(&D);

    /* Build the tilde intermediate */
    dpd_file2_copy(&FMI, CC_OEI, "FMIt");
    dpd_file2_close(&FMI);

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    dpd_file2_close(&FME);
    dpd_file2_close(&tIA);

    dpd_file2_close(&FMIt);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&Fmi, CC_OEI, 0, 0, 0, "Fmi");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);
  
    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);
  
    dpd_buf4_init(&E_anti, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
    dpd_dot13(&tia, &E, &FMI, 1, 1, 1.0, 1.0);

    dpd_dot13(&tia, &E_anti, &Fmi, 1, 1, 1.0, 1.0);
    dpd_dot13(&tIA, &E, &Fmi, 1, 1, 1.0, 1.0);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
    dpd_buf4_close(&E_anti);
    dpd_buf4_close(&E);

    dpd_buf4_init(&D_anti, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_buf4_init(&tautIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    dpd_buf4_init(&tautijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautijab");

    dpd_contract442(&D_anti, &tautIJAB, &FMI, 0, 0, 1, 1);
    dpd_contract442(&D_anti, &tautijab, &Fmi, 0, 0, 1, 1);

    dpd_buf4_close(&tautIJAB);
    dpd_buf4_close(&tautijab);
    dpd_buf4_close(&D_anti);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    dpd_contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_contract442(&D, &tautIjAb, &Fmi, 1, 1, 1, 1);

    dpd_buf4_close(&tautIjAb);
    dpd_buf4_close(&D);

    /* Build the tilde intermediate */
    dpd_file2_copy(&FMI, CC_OEI, "FMIt");
    dpd_file2_copy(&Fmi, CC_OEI, "Fmit");

    dpd_file2_close(&FMI);
    dpd_file2_close(&Fmi);

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_init(&Fmit, CC_OEI, 0, 0, 0, "Fmit");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    dpd_file2_close(&FME);
    dpd_file2_close(&tIA);

    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
    dpd_contract222(&Fme, &tia, &Fmit, 0, 0, 0.5, 1);
    dpd_file2_close(&Fme);
    dpd_file2_close(&tia);

    dpd_file2_close(&FMIt);
    dpd_file2_close(&Fmit);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&Fmi, CC_OEI, 0, 2, 2, "Fmi");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);
  
    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);
  
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&E_anti, CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");

    dpd_dot13(&tIA, &E_anti, &FMI, 1, 1, 1, 1);
    dpd_dot24(&tia, &E, &FMI, 0, 0, 1, 1);

    dpd_buf4_close(&E);
    dpd_buf4_close(&E_anti);

    dpd_buf4_init(&E_anti, CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
    dpd_buf4_init(&E, CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");

    dpd_dot13(&tia, &E_anti, &Fmi, 1, 1, 1, 1);
    dpd_dot13(&tIA, &E, &Fmi, 1, 1, 1, 1);

    dpd_buf4_close(&E);
    dpd_buf4_close(&E_anti);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_buf4_init(&tautIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    dpd_contract442(&D, &tautIJAB, &FMI, 0, 0, 1, 1);
    dpd_buf4_close(&tautIJAB);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_buf4_init(&tautijab, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tautijab");
    dpd_contract442(&D, &tautijab, &Fmi, 0, 0, 1, 1);
    dpd_buf4_close(&tautijab);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_buf4_close(&tautIjAb);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tautiJaB");
    dpd_contract442(&D, &tautIjAb, &Fmi, 0, 0, 1, 1);
    dpd_buf4_close(&tautIjAb);
    dpd_buf4_close(&D);


    /* Build the tilde intermediate */
    dpd_file2_copy(&FMI, CC_OEI, "FMIt");
    dpd_file2_copy(&Fmi, CC_OEI, "Fmit");

    dpd_file2_close(&FMI);
    dpd_file2_close(&Fmi);

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_init(&Fmit, CC_OEI, 0, 2, 2, "Fmit");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
    dpd_file2_close(&FME);
    dpd_file2_close(&tIA);

    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_contract222(&Fme, &tia, &Fmit, 0, 0, 0.5, 1);
    dpd_file2_close(&Fme);
    dpd_file2_close(&tia);

    dpd_file2_close(&FMIt);
    dpd_file2_close(&Fmit);
  }
}
}} // namespace psi::ccenergy
