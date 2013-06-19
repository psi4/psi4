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
    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FAI residual");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "FAI residual");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    dpd_->file2_close(&FAE);

    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    dpd_->file2_close(&FMI);

    dpd_->file2_close(&tIA); 

    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");

    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    dpd_->contract422(&tIjAb, &FME, &newtIA, 0, 0, 1, 1);
    dpd_->buf4_close(&tIjAb);

    dpd_->file2_close(&FME);  

    dpd_->buf4_init(&C_anti, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    dpd_->dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
    dpd_->dot13(&tIA, &D, &newtIA, 0, 0, 1, 1);

    dpd_->file2_close(&tIA);  

    dpd_->buf4_close(&C_anti);  
    dpd_->buf4_close(&D);

    /* t(i,a) <-- (2 t(mi,ef) - t(mi,fe)) <ma|ef> */
    /* out-of-core version replacing the *stupid* code above 3/22/05, TDC */
      dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(ma,mi)");
      dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F 2<ia|bc> - <ia|cb>");
      dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_->contract444(&F, &tIjAb, &Z, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&tIjAb);
      dpd_->buf4_close(&F);
      dpd_->trace42_13(&Z, &newtIA, 1, 1.0, 1.0);
      dpd_->buf4_close(&Z);

    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->contract442(&E, &tIjAb, &newtIA, 1, 3, -1, 1);
    dpd_->buf4_close(&E);  
    dpd_->buf4_close(&tIjAb);

dotval = dpd_->file2_dot_self(&newtIA); 
fprintf(outfile,"\t Norm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",dotval);
    dpd_->file2_close(&newtIA);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FAI residual");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fai residual");
    dpd_->file2_close(&fia);

    dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "FAI residual");
    dpd_->file2_init(&newtia, PSIF_CC_OEI, 0, 0, 1, "Fai residual");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 1, 1, "Fae");

    dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    dpd_->contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);

    dpd_->file2_close(&FAE);  dpd_->file2_close(&Fae);

    dpd_->file2_close(&tIA); dpd_->file2_close(&tia);

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 0, 0, "Fmi");

    dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    dpd_->contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);

    dpd_->file2_close(&FMI);  dpd_->file2_close(&Fmi);
    dpd_->file2_close(&tIA);  dpd_->file2_close(&tia);

    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");

    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");

    dpd_->dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
    dpd_->dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);
  
    dpd_->buf4_close(&tIJAB);  
    dpd_->buf4_close(&tijab);

    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    dpd_->dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
    dpd_->dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);
  
    dpd_->buf4_close(&tIjAb);
  
    dpd_->file2_close(&FME);  
    dpd_->file2_close(&Fme);

    dpd_->buf4_init(&C_anti, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0,
		  "C <ia||jb>");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    dpd_->dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
    dpd_->dot13(&tia, &D, &newtIA, 0, 0, 1, 1);

    dpd_->dot14(&tia, &C_anti, &newtia, 0, 1, -1, 1);
    dpd_->dot13(&tIA, &D, &newtia, 0, 0, 1, 1);

    dpd_->file2_close(&tIA);  
    dpd_->file2_close(&tia);

    dpd_->buf4_close(&C_anti);  
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");

    dpd_->contract442(&tIJAB, &F_anti, &newtIA, 1, 1, 1, 1);
    dpd_->contract442(&tijab, &F_anti, &newtia, 1, 1, 1, 1);

    dpd_->buf4_close(&tIJAB);
    dpd_->buf4_close(&tijab);
    dpd_->buf4_close(&F_anti);

    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

    dpd_->contract442(&tiJaB, &F, &newtIA, 1, 1, 1, 1);
    dpd_->contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);
  
    dpd_->buf4_close(&F);  
    dpd_->buf4_close(&tIjAb);  
    dpd_->buf4_close(&tiJaB);

    dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");

    dpd_->contract442(&E_anti, &tIJAB, &newtIA, 1, 3, -1, 1);
    dpd_->contract442(&E_anti, &tijab, &newtia, 1, 3, -1, 1);

    dpd_->buf4_close(&E_anti);  
    dpd_->buf4_close(&tIJAB);  
    dpd_->buf4_close(&tijab);
  
    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

    dpd_->contract442(&E, &tiJaB, &newtIA, 1, 3, -1, 1);
    dpd_->contract442(&E, &tIjAb, &newtia, 1, 3, -1, 1);

    dpd_->buf4_close(&E);  
    dpd_->buf4_close(&tIjAb);  
    dpd_->buf4_close(&tiJaB);

  /* Purge FAI matrix elements */
  dpd_->file2_mat_init(&newtIA);
  dpd_->file2_mat_rd(&newtIA);
  for(h=0; h < nirreps; h++) {
    for(m=0; m<occpi[h]; m++)
      for(e=(virtpi[h]-openpi[h]); e<virtpi[h]; e++)
        newtIA.matrix[h][m][e] = 0.0;
  }
  dpd_->file2_mat_wrt(&newtIA);
  dpd_->file2_mat_close(&newtIA);

  /* Purge Fai matrix elements */
  dpd_->file2_mat_init(&newtia);
  dpd_->file2_mat_rd(&newtia);
  for(h=0; h < nirreps; h++) {
    for(e=0; e<virtpi[h]; e++)
      for(m=(occpi[h]-openpi[h]); m<occpi[h]; m++)
        newtia.matrix[h][m][e] = 0.0;
  }
  dpd_->file2_mat_wrt(&newtia);
  dpd_->file2_mat_close(&newtia);

dotval = dpd_->file2_dot_self(&newtIA); 
fprintf(outfile,"\t Norm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",dotval);
dotval = dpd_->file2_dot_self(&newtia); 
fprintf(outfile,"\t Norm squared of <Phi_i^a|Hbar|0> = %20.15lf\n",dotval);

    dpd_->file2_close(&newtIA);  dpd_->file2_close(&newtia);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FAI residual");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fai residual");
    dpd_->file2_close(&fia);

    dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "FAI residual");
    dpd_->file2_init(&newtia, PSIF_CC_OEI, 0, 2, 3, "Fai residual");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    dpd_->file2_close(&FAE);  

    dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");
    dpd_->contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);
    dpd_->file2_close(&Fae);

    dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
    dpd_->file2_close(&FMI);  

    dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");
    dpd_->contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);
    dpd_->file2_close(&Fmi);

    dpd_->file2_close(&tIA);  
    dpd_->file2_close(&tia);

    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");


    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_->dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
    dpd_->buf4_close(&tIJAB);  

    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_->dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);
    dpd_->buf4_close(&tijab);

    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
    dpd_->dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);
    dpd_->buf4_close(&tIjAb);
  
    dpd_->file2_close(&FME);  
    dpd_->file2_close(&Fme);



    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_->dot14(&tIA, &C, &newtIA, 0, 1, -1, 1);
    dpd_->buf4_close(&C);

    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_->dot14(&tia, &C, &newtia, 0, 1, -1, 1);
    dpd_->buf4_close(&C);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_->dot13(&tia, &D, &newtIA, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->dot13(&tIA, &D, &newtia, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&tIA);  
    dpd_->file2_close(&tia);


    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->contract442(&tIJAB, &F, &newtIA, 1, 1, 1, 1);
    dpd_->buf4_close(&tIJAB);
    dpd_->buf4_close(&F);

    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_->contract442(&tijab, &F, &newtia, 1, 1, 1, 1);
    dpd_->buf4_close(&tijab);
    dpd_->buf4_close(&F);

    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->contract442(&tIjAb, &F, &newtIA, 0, 2, 1, 1);
    dpd_->buf4_close(&tIjAb);
    dpd_->buf4_close(&F);  
  
    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);
    dpd_->buf4_close(&tIjAb);  
    dpd_->buf4_close(&F);  



    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_->contract442(&E, &tIJAB, &newtIA, 1, 3, -1, 1);
    dpd_->buf4_close(&E);  

    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_->contract442(&E, &tijab, &newtia, 1, 3, -1, 1);
    dpd_->buf4_close(&E);  

    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->contract442(&E, &tIjAb, &newtIA, 2, 2, -1, 1);
    dpd_->buf4_close(&E);  
    dpd_->buf4_close(&tIjAb);  

    dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_->contract442(&E, &tIjAb, &newtia, 2, 2, -1, 1);
    dpd_->buf4_close(&E);  
    dpd_->buf4_close(&tIjAb);  

dotval = dpd_->file2_dot_self(&newtIA); 
fprintf(outfile,"\t Norm squared of <Phi_I^A|Hbar|0> = %20.15lf\n",dotval);
dotval = dpd_->file2_dot_self(&newtia); 
fprintf(outfile,"\t Norm squared of <Phi_i^a|Hbar|0> = %20.15lf\n",dotval);

    dpd_->file2_close(&newtIA);  dpd_->file2_close(&newtia);
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
    dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    dpd_->file2_close(&fAB);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    dpd_->file2_close(&fAB);

    dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
    dpd_->file2_copy(&fab, PSIF_CC_OEI, "Fae");
    dpd_->file2_close(&fab);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    dpd_->file2_close(&fAB);

    dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
    dpd_->file2_copy(&fab, PSIF_CC_OEI, "Fae");
    dpd_->file2_close(&fab);
  }

/* don't remove diagonal elements here */

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&fIA);
    dpd_->file2_close(&FAE);

    /* Out-of-core algorithm for F->FAE added 3/20/05 - TDC */
    /* Fae <-- t(m,f) [2 <ma|fe> - <ma|ef>] */
    dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    dpd_->file2_mat_init(&FAE);
    dpd_->file2_mat_rd(&FAE);
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_mat_init(&tIA);
    dpd_->file2_mat_rd(&tIA);
    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    for(Gma=0; Gma < nirreps; Gma++) {
      dpd_->buf4_mat_irrep_row_init(&F, Gma);
      X = init_array(F.params->coltot[Gma]);

      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {
	dpd_->buf4_mat_irrep_row_rd(&F, Gma, ma);
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
      dpd_->buf4_mat_irrep_row_close(&F, Gma);
    }
    dpd_->buf4_close(&F);
    dpd_->file2_mat_close(&tIA);
    dpd_->file2_close(&tIA);
    dpd_->file2_mat_wrt(&FAE);
    dpd_->file2_mat_close(&FAE);
    dpd_->file2_close(&FAE);

    dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    dpd_->contract442(&tautIjAb, &D, &FAE, 3, 3, -1, 1);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&tautIjAb);

    /* Build the tilde intermediates */
    dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    dpd_->file2_close(&FAE);

    dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&FME);
  
    dpd_->file2_close(&FAEt);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 1, 1, "Fae");

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    dpd_->file2_close(&tia);
    dpd_->file2_close(&fia);

    dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    dpd_->dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
    dpd_->dot13(&tia, &F, &FAE, 0, 0, 1.0, 1.0);

    dpd_->dot13(&tia, &F_anti, &Fae, 0, 0, 1.0, 1.0);
    dpd_->dot13(&tIA, &F, &Fae, 0, 0, 1.0, 1.0);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);
    dpd_->buf4_close(&F_anti);
    dpd_->buf4_close(&F);

    dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");

    dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautijab");

    dpd_->contract442(&tautIJAB, &D_anti, &FAE, 2, 2, -1, 1);
    dpd_->contract442(&tautijab, &D_anti, &Fae, 2, 2, -1, 1);

    dpd_->buf4_close(&D_anti);
    dpd_->buf4_close(&tautIJAB);
    dpd_->buf4_close(&tautijab);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    dpd_->contract442(&tautIjAb, &D, &Fae, 3, 3, -1, 1);
    dpd_->contract442(&tautIjAb, &D, &FAE, 2, 2, -1, 1);

    dpd_->buf4_close(&D);
    dpd_->buf4_close(&tautIjAb);


    /* Build the tilde intermediates */
    dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    dpd_->file2_copy(&Fae, PSIF_CC_OEI, "Faet");

    dpd_->file2_close(&FAE);
    dpd_->file2_close(&Fae);

    dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 1, 1, "Faet");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&FME);
  
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
    dpd_->contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
    dpd_->file2_close(&tia);
    dpd_->file2_close(&Fme);

    dpd_->file2_close(&FAEt);
    dpd_->file2_close(&Faet); 
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    dpd_->file2_close(&tia);
    dpd_->file2_close(&fia);

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_->dot13(&tIA, &F, &FAE, 0, 0, 1, 1);
    dpd_->buf4_close(&F);

    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_->dot13(&tia, &F, &FAE, 0, 0, 1, 1);
    dpd_->buf4_close(&F);

    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_->dot13(&tia, &F, &Fae, 0, 0, 1, 1);
    dpd_->buf4_close(&F);

    dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_->dot13(&tIA, &F, &Fae, 0, 0, 1, 1);
    dpd_->buf4_close(&F);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    dpd_->contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    dpd_->buf4_close(&taut);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_->contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    dpd_->buf4_close(&taut);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tautijab");
    dpd_->contract442(&taut, &D, &Fae, 2, 2, -1, 1);
    dpd_->buf4_close(&taut);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_->contract442(&taut, &D, &Fae, 3, 3, -1, 1);
    dpd_->buf4_close(&taut);
    dpd_->buf4_close(&D);

    /* Build the tilde intermediates */
    dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    dpd_->file2_copy(&Fae, PSIF_CC_OEI, "Faet");

    dpd_->file2_close(&FAE);
    dpd_->file2_close(&Fae);

    dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 3, 3, "Faet");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_->file2_close(&tIA);
    dpd_->file2_close(&FME);
  
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    dpd_->contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
    dpd_->file2_close(&tia);
    dpd_->file2_close(&Fme);


    dpd_->file2_close(&FAEt);
    dpd_->file2_close(&Faet); 
  }
}

void Fme_for_Fai(void)
{
  dpdfile2 FME, Fme, fIA, fia, tIA, tia;
  dpdbuf4 D_anti, D;

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
  
    dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    dpd_->dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    dpd_->dot13(&tIA, &D, &FME, 0, 0, 1.0, 1.0);

    dpd_->file2_close(&tIA);
    dpd_->buf4_close(&D_anti);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&FME);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fme");
    dpd_->file2_close(&fia);
  
    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
  
    dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    dpd_->dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    dpd_->dot13(&tia, &D, &FME, 0, 0, 1.0, 1.0);

    dpd_->dot13(&tia, &D_anti, &Fme, 0, 0, 1.0, 1.0);
    dpd_->dot13(&tIA, &D, &Fme, 0, 0, 1.0, 1.0);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);
    dpd_->buf4_close(&D_anti);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&FME);
    dpd_->file2_close(&Fme);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fme");
    dpd_->file2_close(&fia);

    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_->contract422(&D, &tIA, &FME, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_->contract422(&D, &tia, &FME, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_->contract422(&D, &tia, &Fme, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_->contract422(&D, &tIA, &Fme, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);

    dpd_->file2_close(&FME);
    dpd_->file2_close(&Fme);

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

/* don't remove diagonals */

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

}} // namespace psi::cchbar
