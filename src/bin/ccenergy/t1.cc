/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void local_filter_T1(dpdfile2 *T1);

void t1_build(void)
{
  dpdfile2 newtIA, newtia, tIA, tia, fIA, fia;
  dpdfile2 FAE, Fae, FMI, Fmi, FME, Fme;
  dpdfile2 dIA, dia;
  dpdbuf4 tIJAB, tijab, tIjAb, tiJaB, T2;
  dpdbuf4 C, C_anti, D, F_anti, F, E_anti, E, Z;
  int Gma, Gmi, Gm, Gi, Ga, ma, m, a, A, nrows, ncols;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "New tIA");
    dpd_file2_close(&fIA);

    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    dpd_file2_close(&FAE);

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    dpd_file2_close(&FMI);

    dpd_file2_close(&tIA); 

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");

    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    dpd_contract422(&tIjAb, &FME, &newtIA, 0, 0, 1, 1);
    dpd_buf4_close(&tIjAb);

    dpd_file2_close(&FME);

    dpd_buf4_init(&C_anti, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
    dpd_dot13(&tIA, &D, &newtIA, 0, 0, 1, 1);

    dpd_file2_close(&tIA);

    dpd_buf4_close(&C_anti);
    dpd_buf4_close(&D);

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
    dpd_file2_mat_init(&newtIA);
    dpd_file2_mat_rd(&newtIA);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    for(Gma=0; Gma < moinfo.nirreps; Gma++) {
      Gmi = Gma; /* T1 is totally symmetric */

      dpd_buf4_mat_irrep_row_init(&F, Gma);
      dpd_buf4_mat_irrep_init(&T2, Gmi);
      dpd_buf4_mat_irrep_rd(&T2, Gmi);

      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {

	dpd_buf4_mat_irrep_row_rd(&F, Gma, ma);

	m = F.params->roworb[Gma][ma][0];
	a = F.params->roworb[Gma][ma][1];
	Gm = F.params->psym[m];
	Ga = F.params->qsym[a];
	Gi = Ga; /* T1 is totally symmetric */
	A = a - F.params->qoff[Ga];

	nrows = moinfo.occpi[Gi];
	ncols = F.params->coltot[Gma];

	if(nrows && ncols && moinfo.virtpi[Ga])
	  C_DGEMV('n',nrows,ncols,1.0,T2.matrix[Gmi][T2.row_offset[Gmi][m]],ncols,
		  F.matrix[Gma][0],1,1.0,&newtIA.matrix[Gi][0][A],moinfo.virtpi[Ga]);

      }

      dpd_buf4_mat_irrep_close(&T2, Gmi);
      dpd_buf4_mat_irrep_row_close(&F, Gma);
    }
    dpd_buf4_close(&F);
    dpd_buf4_close(&T2);
    dpd_file2_mat_wrt(&newtIA);
    dpd_file2_mat_close(&newtIA);

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&E, &tIjAb, &newtIA, 1, 3, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&tIjAb);

    if (params.just_residuals) {
      dpd_file2_close(&newtIA);
      return;
    }

//    dpd_file2_copy(&newtIA, CC_OEI, "New tIA Increment");
    dpd_file2_close(&newtIA);

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
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "New tIA");
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_copy(&fia, CC_OEI, "New tia");
    dpd_file2_close(&fia);

    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&newtia, CC_OEI, 0, 0, 1, "New tia");


    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 1, 1, "Fae");

    dpd_contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    dpd_contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);

    dpd_file2_close(&FAE);  dpd_file2_close(&Fae);

    dpd_file2_close(&tIA); dpd_file2_close(&tia);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&Fmi, CC_OEI, 0, 0, 0, "Fmi");

    dpd_contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    dpd_contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);

    dpd_file2_close(&FMI);  dpd_file2_close(&Fmi);
    dpd_file2_close(&tIA);  dpd_file2_close(&tia);

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");

    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");

    dpd_dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
    dpd_dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);
  
    dpd_buf4_close(&tIJAB);  
    dpd_buf4_close(&tijab);

    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    dpd_dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
    dpd_dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);
  
    dpd_buf4_close(&tIjAb);
  
    dpd_file2_close(&FME);  
    dpd_file2_close(&Fme);

    dpd_buf4_init(&C_anti, CC_CINTS, 0, 10, 10, 10, 10, 0,
		  "C <ia||jb>");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
    dpd_dot13(&tia, &D, &newtIA, 0, 0, 1, 1);

    dpd_dot14(&tia, &C_anti, &newtia, 0, 1, -1, 1);
    dpd_dot13(&tIA, &D, &newtia, 0, 0, 1, 1);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);

    dpd_buf4_close(&C_anti);  
    dpd_buf4_close(&D);

    dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");

    dpd_contract442(&tIJAB, &F_anti, &newtIA, 1, 1, 1, 1);
    dpd_contract442(&tijab, &F_anti, &newtia, 1, 1, 1, 1);

    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&F_anti);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&tiJaB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

    dpd_contract442(&tiJaB, &F, &newtIA, 1, 1, 1, 1);
    dpd_contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);
  
    dpd_buf4_close(&F);  
    dpd_buf4_close(&tIjAb);  
    dpd_buf4_close(&tiJaB);

    dpd_buf4_init(&E_anti, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");

    dpd_contract442(&E_anti, &tIJAB, &newtIA, 1, 3, -1, 1);
    dpd_contract442(&E_anti, &tijab, &newtia, 1, 3, -1, 1);

    dpd_buf4_close(&E_anti);  
    dpd_buf4_close(&tIJAB);  
    dpd_buf4_close(&tijab);
  
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&tiJaB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

    dpd_contract442(&E, &tiJaB, &newtIA, 1, 3, -1, 1);
    dpd_contract442(&E, &tIjAb, &newtia, 1, 3, -1, 1);

    dpd_buf4_close(&E);  
    dpd_buf4_close(&tIjAb);  
    dpd_buf4_close(&tiJaB);

/*
    dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newtIA);
    dpd_file2_close(&dIA);

    dpd_file2_init(&dia, CC_OEI, 0, 0, 1, "dia");
    dpd_file2_dirprd(&dia, &newtia);
    dpd_file2_close(&dia);
*/

    dpd_file2_close(&newtIA);  dpd_file2_close(&newtia);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "New tIA");
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_copy(&fia, CC_OEI, "New tia");
    dpd_file2_close(&fia);

    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&newtia, CC_OEI, 0, 2, 3, "New tia");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    dpd_file2_close(&FAE);  

    dpd_file2_init(&Fae, CC_OEI, 0, 3, 3, "Fae");
    dpd_contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);
    dpd_file2_close(&Fae);

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    dpd_file2_close(&FMI);  

    dpd_file2_init(&Fmi, CC_OEI, 0, 2, 2, "Fmi");
    dpd_contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);
    dpd_file2_close(&Fmi);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");


    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
    dpd_buf4_close(&tIJAB);  

    dpd_buf4_init(&tijab, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);
    dpd_buf4_close(&tijab);

    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
    dpd_dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);
    dpd_buf4_close(&tIjAb);
  
    dpd_file2_close(&FME);  
    dpd_file2_close(&Fme);



    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_dot14(&tIA, &C, &newtIA, 0, 1, -1, 1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_dot14(&tia, &C, &newtia, 0, 1, -1, 1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_dot13(&tia, &D, &newtIA, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_dot13(&tIA, &D, &newtia, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);


    dpd_buf4_init(&F, CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&tIJAB, &F, &newtIA, 1, 1, 1, 1);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_contract442(&tijab, &F, &newtia, 1, 1, 1, 1);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&tIjAb, &F, &newtIA, 0, 2, 1, 1);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&F);  
  
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);
    dpd_buf4_close(&tIjAb);  
    dpd_buf4_close(&F);  



    dpd_buf4_init(&E, CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&E, &tIJAB, &newtIA, 1, 3, -1, 1);
    dpd_buf4_close(&E);  

    dpd_buf4_init(&E, CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_contract442(&E, &tijab, &newtia, 1, 3, -1, 1);
    dpd_buf4_close(&E);  

    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&E, &tIjAb, &newtIA, 2, 2, -1, 1);
    dpd_buf4_close(&E);  
    dpd_buf4_close(&tIjAb);  

    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&E, &tIjAb, &newtia, 2, 2, -1, 1);
    dpd_buf4_close(&E);  
    dpd_buf4_close(&tIjAb);  

/*
    dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newtIA);
    dpd_file2_close(&dIA);

    dpd_file2_init(&dia, CC_OEI, 0, 2, 3, "dia");
    dpd_file2_dirprd(&dia, &newtia);
    dpd_file2_close(&dia);
*/

    dpd_file2_close(&newtIA);  
    dpd_file2_close(&newtia);

  }
}
}} // namespace psi::ccenergy
