/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* This function computes <phi_i^a | Hbar | 0>.  This will be zero as long
as the T amplitudes were obtained from CCSD computation with the same Hbar. */

void Fme_for_Fai();
void Fae_for_Fai();
void Fmi_for_Fai();

void Fai_build(void)
{
  dpdfile2 newtIA, newtia, tIA, tia, fIA, fia;
  dpdfile2 FAE, Fae, FMI, Fmi, FME, Fme;
  dpdfile2 dIA, dia;
  dpdbuf4 tIJAB, tijab, tIjAb, tiJaB, T2;
  dpdbuf4 C, C_anti, D, F_anti, F, E_anti, E, Z;
  int Gma, Gmi, Gm, Gi, Ga, ma, m, a, A, nrows, ncols, h, e, nirreps;
  int *occpi, *virtpi, *openpi;
  double dotval;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;

  Fme_for_Fai();
  Fae_for_Fai();
  Fmi_for_Fai();

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FAI residual");
    dpd_file2_close(&fIA);

    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "FAI residual");

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

    /* t(i,a) <-- (2 t(mi,ef) - t(mi,fe)) <ma|ef> */
    /* out-of-core version replacing the *stupid* code above 3/22/05, TDC */
      dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(ma,mi)");
      dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F 2<ia|bc> - <ia|cb>");
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_contract444(&F, &tIjAb, &Z, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_close(&F);
      dpd_trace42_13(&Z, &newtIA, 1, 1.0, 1.0);
      dpd_buf4_close(&Z);

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&E, &tIjAb, &newtIA, 1, 3, -1, 1);
    dpd_buf4_close(&E);  
    dpd_buf4_close(&tIjAb);

dotval = dpd_file2_dot_self(&newtIA); 
fprintf(outfile,"\t Norm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",dotval);
    dpd_file2_close(&newtIA);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FAI residual");
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_copy(&fia, CC_OEI, "Fai residual");
    dpd_file2_close(&fia);

    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "FAI residual");
    dpd_file2_init(&newtia, CC_OEI, 0, 0, 1, "Fai residual");

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

  /* Purge FAI matrix elements */
  dpd_file2_mat_init(&newtIA);
  dpd_file2_mat_rd(&newtIA);
  for(h=0; h < nirreps; h++) {
    for(m=0; m<occpi[h]; m++)
      for(e=(virtpi[h]-openpi[h]); e<virtpi[h]; e++)
        newtIA.matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(&newtIA);
  dpd_file2_mat_close(&newtIA);

  /* Purge Fai matrix elements */
  dpd_file2_mat_init(&newtia);
  dpd_file2_mat_rd(&newtia);
  for(h=0; h < nirreps; h++) {
    for(e=0; e<virtpi[h]; e++)
      for(m=(occpi[h]-openpi[h]); m<occpi[h]; m++)
        newtia.matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(&newtia);
  dpd_file2_mat_close(&newtia);

dotval = dpd_file2_dot_self(&newtIA); 
fprintf(outfile,"\t Norm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",dotval);
dotval = dpd_file2_dot_self(&newtia); 
fprintf(outfile,"\t Norm squared of <Phi_i^a|Hbar|0> = %20.15lf\n",dotval);

    dpd_file2_close(&newtIA);  dpd_file2_close(&newtia);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FAI residual");
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_copy(&fia, CC_OEI, "Fai residual");
    dpd_file2_close(&fia);

    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "FAI residual");
    dpd_file2_init(&newtia, CC_OEI, 0, 2, 3, "Fai residual");

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

dotval = dpd_file2_dot_self(&newtIA); 
fprintf(outfile,"\t Norm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",dotval);
dotval = dpd_file2_dot_self(&newtia); 
fprintf(outfile,"\t Norm squared of <Phi_i^a|Hbar|0> = %20.15lf\n",dotval);

    dpd_file2_close(&newtIA);  dpd_file2_close(&newtia);
  }
}


void Fae_for_Fai(void)
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

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_copy(&fAB, CC_OEI, "FAE");
    dpd_file2_close(&fAB);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_copy(&fAB, CC_OEI, "FAE");
    dpd_file2_close(&fAB);

    dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
    dpd_file2_copy(&fab, CC_OEI, "Fae");
    dpd_file2_close(&fab);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_copy(&fAB, CC_OEI, "FAE");
    dpd_file2_close(&fAB);

    dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
    dpd_file2_copy(&fab, CC_OEI, "Fae");
    dpd_file2_close(&fab);
  }

/* don't remove diagonal elements here */

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);
    dpd_file2_close(&FAE);

    /* Out-of-core algorithm for F->FAE added 3/20/05 - TDC */
    /* Fae <-- t(m,f) [2 <ma|fe> - <ma|ef>] */
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&tIA);
    dpd_file2_mat_rd(&tIA);
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    for(Gma=0; Gma < nirreps; Gma++) {
      dpd_buf4_mat_irrep_row_init(&F, Gma);
      X = init_array(F.params->coltot[Gma]);

      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {
	dpd_buf4_mat_irrep_row_rd(&F, Gma, ma);
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
	
	nrows = moinfo.virtpi[Gf];
	ncols = moinfo.virtpi[Ge];
	if(nrows && ncols)
	  C_DGEMV('t',nrows,ncols,1.0,&X[F.col_offset[Gma][Gf]],ncols,
		  tIA.matrix[Gm][M],1,1.0,
		  FAE.matrix[Ga][A],1);
      }

      free(X);
      dpd_buf4_mat_irrep_row_close(&F, Gma);
    }
    dpd_buf4_close(&F);
    dpd_file2_mat_close(&tIA);
    dpd_file2_close(&tIA);
    dpd_file2_mat_wrt(&FAE);
    dpd_file2_mat_close(&FAE);
    dpd_file2_close(&FAE);

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    dpd_contract442(&tautIjAb, &D, &FAE, 3, 3, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tautIjAb);

    /* Build the tilde intermediates */
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_close(&FAE);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_close(&FAEt);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 1, 1, "Fae");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);

    dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
    dpd_dot13(&tia, &F, &FAE, 0, 0, 1.0, 1.0);

    dpd_dot13(&tia, &F_anti, &Fae, 0, 0, 1.0, 1.0);
    dpd_dot13(&tIA, &F, &Fae, 0, 0, 1.0, 1.0);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
    dpd_buf4_close(&F_anti);
    dpd_buf4_close(&F);

    dpd_buf4_init(&D_anti, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");

    dpd_buf4_init(&tautIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    dpd_buf4_init(&tautijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautijab");

    dpd_contract442(&tautIJAB, &D_anti, &FAE, 2, 2, -1, 1);
    dpd_contract442(&tautijab, &D_anti, &Fae, 2, 2, -1, 1);

    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&tautIJAB);
    dpd_buf4_close(&tautijab);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    dpd_contract442(&tautIjAb, &D, &Fae, 3, 3, -1, 1);
    dpd_contract442(&tautIjAb, &D, &FAE, 2, 2, -1, 1);

    dpd_buf4_close(&D);
    dpd_buf4_close(&tautIjAb);


    /* Build the tilde intermediates */
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_copy(&Fae, CC_OEI, "Faet");

    dpd_file2_close(&FAE);
    dpd_file2_close(&Fae);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&Faet, CC_OEI, 0, 1, 1, "Faet");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
    dpd_contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&Fme);

    dpd_file2_close(&FAEt);
    dpd_file2_close(&Faet); 
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 3, 3, "Fae");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_dot13(&tIA, &F, &FAE, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_dot13(&tia, &F, &FAE, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_dot13(&tia, &F, &Fae, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_dot13(&tIA, &F, &Fae, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    dpd_contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tautijab");
    dpd_contract442(&taut, &D, &Fae, 2, 2, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_contract442(&taut, &D, &Fae, 3, 3, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    /* Build the tilde intermediates */
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_copy(&Fae, CC_OEI, "Faet");

    dpd_file2_close(&FAE);
    dpd_file2_close(&Fae);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&Faet, CC_OEI, 0, 3, 3, "Faet");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&Fme);


    dpd_file2_close(&FAEt);
    dpd_file2_close(&Faet); 
  }
}

void Fme_for_Fai(void)
{
  dpdfile2 FME, Fme, fIA, fia, tIA, tia;
  dpdbuf4 D_anti, D;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FME");
    dpd_file2_close(&fIA);

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  
    dpd_buf4_init(&D_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    dpd_dot13(&tIA, &D, &FME, 0, 0, 1.0, 1.0);

    dpd_file2_close(&tIA);
    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&D);

    dpd_file2_close(&FME);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FME");
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_copy(&fia, CC_OEI, "Fme");
    dpd_file2_close(&fia);
  
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
  
    dpd_buf4_init(&D_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    dpd_dot13(&tia, &D, &FME, 0, 0, 1.0, 1.0);

    dpd_dot13(&tia, &D_anti, &Fme, 0, 0, 1.0, 1.0);
    dpd_dot13(&tIA, &D, &Fme, 0, 0, 1.0, 1.0);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&D);

    dpd_file2_close(&FME);
    dpd_file2_close(&Fme);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FME");
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_copy(&fia, CC_OEI, "Fme");
    dpd_file2_close(&fia);

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_contract422(&D, &tIA, &FME, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_contract422(&D, &tia, &FME, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_contract422(&D, &tia, &Fme, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_contract422(&D, &tIA, &Fme, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_file2_close(&FME);
    dpd_file2_close(&Fme);

  }
}

void Fmi_for_Fai(void)
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

/* don't remove diagonals */

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

}} // namespace psi::cchbar
