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
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void local_filter_T1(dpdfile2 *T1);

void cc2_L1_build(struct L_Params L_params) {

  int GW, GL1, GL2, Gab, Gij, Gei, Gi, Ga, Gm;
  int a, A, i, I, ab, nlinks, nrows, ncols;
  dpdfile2 newLIA, newLia, LIA, Lia;
  dpdfile2 dIA, dia, Fme, FME;
  dpdfile2 Fae, FAE, Fmi, FMI;
  dpdfile2 GMI, Gmi, Gae, XIA, Xia;
  dpdfile2 GAE, G, GAI, Gai;
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ;
  dpdbuf4 WMBIJ, Wmbij, WMbIj, WmBiJ;
  dpdbuf4 LIJAB, Lijab, LIjAb, LiJaB, L2;
  dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE;
  dpdbuf4 WAMEF, Wamef, WAmEf, WaMeF, W;
  dpdbuf4 Z, D, E;
  dpdfile2 XLD;
  int L_irr;
  L_irr = L_params.irrep;

  /* ground state inhomogeneous term is Fme */
  if (L_params.ground) {
    if(params.ref == 0) {
      global_dpd_->file2_init(&FME,PSIF_CC_OEI, 0, 0, 1, "FME");
      global_dpd_->file2_copy(&FME, PSIF_CC_LAMBDA, "New LIA");
      global_dpd_->file2_close(&FME);
    }
    else if(params.ref == 1) {
      global_dpd_->file2_init(&Fme,PSIF_CC_OEI, 0, 0, 1, "Fme");
      global_dpd_->file2_init(&FME,PSIF_CC_OEI, 0, 0, 1, "FME");
      global_dpd_->file2_copy(&Fme, PSIF_CC_LAMBDA, "New Lia");
      global_dpd_->file2_copy(&FME, PSIF_CC_LAMBDA, "New LIA");
      global_dpd_->file2_close(&Fme);
      global_dpd_->file2_close(&FME);
    }
    else if(params.ref == 2) {
      global_dpd_->file2_init(&Fme,PSIF_CC_OEI, 0, 2, 3, "Fme");
      global_dpd_->file2_init(&FME,PSIF_CC_OEI, 0, 0, 1, "FME");
      global_dpd_->file2_copy(&Fme, PSIF_CC_LAMBDA, "New Lia");
      global_dpd_->file2_copy(&FME, PSIF_CC_LAMBDA, "New LIA");
      global_dpd_->file2_close(&Fme);
      global_dpd_->file2_close(&FME);
    }
  }
  /* excited state - no inhomogenous term, first term is -energy*L*/
  else if (!params.zeta) {
    if (params.ref == 0 || params.ref == 1) {
      global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
      global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
      global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
      global_dpd_->file2_axpy(&LIA, &newLIA, -1.0 * L_params.cceom_energy, 0);
      global_dpd_->file2_axpy(&Lia, &newLia, -1.0 * L_params.cceom_energy, 0);
      global_dpd_->file2_close(&LIA);
      global_dpd_->file2_close(&newLIA);
      global_dpd_->file2_close(&Lia);
      global_dpd_->file2_close(&newLia);
    }
    else if (params.ref == 2) {
      /* do nothing - TDC did not change to increments for the UHF case */
    }
  }

  if(params.ref == 0) {
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");

    /* L1 RHS += Lie*Fea */
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&LIA,&FAE,&newLIA, 0, 1, 1, 1);
    global_dpd_->file2_close(&FAE);

    /* L1 RHS += -Lma*Fim */
    global_dpd_->file2_init(&FMI,PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract222(&FMI,&LIA,&newLIA, 0, 1, -1, 1);
    global_dpd_->file2_close(&FMI);

    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");

    /* L1 RHS += Lie*Fea */
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");
    global_dpd_->contract222(&Lia,&Fae,&newLia, 0, 1, 1, 1);
    global_dpd_->contract222(&LIA,&FAE,&newLIA, 0, 1, 1, 1);
    global_dpd_->file2_close(&Fae);
    global_dpd_->file2_close(&FAE);

    /* L1 RHS += -Lma*Fim */
    global_dpd_->file2_init(&FMI,PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi,PSIF_CC_OEI, 0, 2, 2, "Fmi");
    global_dpd_->contract222(&Fmi,&Lia,&newLia, 0, 1, -1, 1);
    global_dpd_->contract222(&FMI,&LIA,&newLIA, 0, 1, -1, 1);
    global_dpd_->file2_close(&Fmi);
    global_dpd_->file2_close(&FMI);

    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&Lia);

    global_dpd_->file2_close(&newLIA);
    global_dpd_->file2_close(&newLia);
  }

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");

    /* L1 RHS += Lme*Wieam */
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->contract422(&W, &LIA, &newLIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");

    global_dpd_->buf4_init(&WMBEJ, PSIF_CC2_HET1, 0, 20, 21, 20, 21, 0, "CC2 WMBEJ (ME,BJ)");
    global_dpd_->contract422(&WMBEJ, &LIA, &newLIA, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);

    global_dpd_->buf4_init(&Wmbej, PSIF_CC2_HET1, 0, 30, 31, 30, 31, 0, "CC2 Wmbej (me,bj)");
    global_dpd_->contract422(&Wmbej, &Lia, &newLia, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Wmbej);

    global_dpd_->buf4_init(&WMbEj, PSIF_CC2_HET1, 0, 20, 31, 20, 31, 0, "CC2 WMbEj (ME,bj)");
    global_dpd_->contract422(&WMbEj, &Lia, &newLIA, 1, 0, -1, 1);
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WmBeJ, PSIF_CC2_HET1, 0, 30, 21, 30, 21, 0, "CC2 WmBeJ (me,BJ)");
    global_dpd_->contract422(&WmBeJ, &LIA, &newLia, 1, 0, -1, 1);
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&Lia);

    global_dpd_->file2_close(&newLIA);
    global_dpd_->file2_close(&newLia);

  }

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_TMP0, rspq, 5, 0, "Z (2 AbIj - AbjI)");
    global_dpd_->buf4_close(&L2);

    /* L1 RHS += 1/2 Limef*Wefam */
    /* Out-of-core contract442 */

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
    global_dpd_->buf4_init(&L2, PSIF_CC_TMP0, L_irr, 5, 0, 5, 0, 0, "Z (2 AbIj - AbjI)");

    /* dpd_contract442(&L2, &W, &newLIA, 0, 0, 1, 1); */

    GW = W.file.my_irrep;
    GL2 = L2.file.my_irrep;
    GL1 = newLIA.my_irrep;

    global_dpd_->file2_mat_init(&newLIA);
    global_dpd_->file2_mat_rd(&newLIA);

    for(Gab=0; Gab < moinfo.nirreps; Gab++) {

      global_dpd_->buf4_mat_irrep_row_init(&L2, Gab);
      global_dpd_->buf4_mat_irrep_row_init(&W, Gab);

      for(ab=0; ab < L2.params->rowtot[Gab]; ab++) {

	global_dpd_->buf4_mat_irrep_row_zero(&L2, Gab, ab);
	global_dpd_->buf4_mat_irrep_row_rd(&L2, Gab, ab);

	global_dpd_->buf4_mat_irrep_row_zero(&W, Gab, ab);
	global_dpd_->buf4_mat_irrep_row_rd(&W, Gab, ab);

	for(Gi=0; Gi < moinfo.nirreps; Gi++) {
	  Ga = Gi^GL1;
	  Gm = GL2^Gab^Gi;

	  nrows = L2.params->rpi[Gi];
	  ncols = W.params->rpi[Ga];
	  nlinks = L2.params->spi[Gm];

	  if(nrows && ncols && nlinks) {
	    C_DGEMM('n','t',nrows,ncols,nlinks,1.0,
		    &(L2.matrix[Gab][0][L2.col_offset[Gab][Gi]]),nlinks,
		    &(W.matrix[Gab][0][W.col_offset[Gab][Ga]]),nlinks,1.0,
		    &(newLIA.matrix[Gi][0][0]),ncols);
	  }
	}
      }
      global_dpd_->buf4_mat_irrep_row_close(&L2, Gab);
      global_dpd_->buf4_mat_irrep_row_close(&W, Gab);
    }
    global_dpd_->file2_mat_wrt(&newLIA);
    global_dpd_->file2_mat_close(&newLIA);

    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 21, 7, 21, 7, 0, "CC2 WABEI (EI,A>B)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 26, 28, 26, 28, 0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 31, 17, 31, 17, 0, "CC2 Wabei (ei,a>b)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&L2, &W, &newLia, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 25, 29, 25, 29, 0, "CC2 WaBeI (eI,aB)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&L2, &W, &newLia, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_close(&newLIA);
    global_dpd_->file2_close(&newLia);
  }

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");

    /* L1 RHS += -1/2 Lmnae*Wiemn */
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&W, &LIjAb, &newLIA, 0, 2, -1, 1);
    global_dpd_->buf4_close(&LIjAb);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");

    /* L1 RHS += -1/2 Lmnae*Wiemn */
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC2_HET1, 0, 20, 2, 20, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1, 1);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_close(&WMBIJ);

    global_dpd_->buf4_init(&WMbIj, PSIF_CC2_HET1, 0, 24, 22, 24, 22, 0, "CC2 WMbIj (Mb,Ij)");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1, 1);
    global_dpd_->buf4_close(&LIjAb);
    global_dpd_->buf4_close(&WMbIj);

    global_dpd_->buf4_init(&Wmbij, PSIF_CC2_HET1, 0, 30, 12, 30, 12, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1, 1);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_close(&Wmbij);

    global_dpd_->buf4_init(&WmBiJ, PSIF_CC2_HET1, 0, 27, 23, 27, 23, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_init(&LiJaB, PSIF_CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1, 1);
    global_dpd_->buf4_close(&LiJaB);
    global_dpd_->buf4_close(&WmBiJ);

    global_dpd_->file2_close(&newLIA);
    global_dpd_->file2_close(&newLia);
  }

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");

    /* L1 RHS += Gbj*<ij|ab> */
    global_dpd_->file2_init(&G, PSIF_CC_TMP0, L_irr, 1, 0, "CC2 GAI");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract422(&D, &G, &newLIA, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&G);

    global_dpd_->file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");

    global_dpd_->file2_init(&GAI, PSIF_CC_LAMBDA, L_irr, 1, 0, "CC2 GAI");
    global_dpd_->file2_init(&Gai, PSIF_CC_LAMBDA, L_irr, 3, 2, "CC2 Gai");

    /* L1 RHS += Gbj*<ij|ab> */
    /** AA **/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->contract422(&D, &GAI, &newLIA, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->contract422(&D, &Gai, &newLIA, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    /** BB**/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->contract422(&D, &Gai, &newLia, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->contract422(&D, &GAI, &newLia, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&Gai);
    global_dpd_->file2_close(&GAI);

    global_dpd_->file2_close(&newLIA);
    global_dpd_->file2_close(&newLia);
  }

  if(params.ref == 0) { /** RHF **/
    /* newLia * Dia */
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_copy(&newLIA, PSIF_CC_LAMBDA, "New LIA Increment");
    global_dpd_->file2_close(&newLIA);

    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA Increment");
    if(params.local && local.filter_singles) local_filter_T1(&newLIA);
    else {
      global_dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
      global_dpd_->file2_dirprd(&dIA, &newLIA);
      global_dpd_->file2_close(&dIA);
    }
    global_dpd_->file2_close(&newLIA);

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_copy(&LIA, PSIF_CC_LAMBDA, "New LIA");
    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA Increment");
    global_dpd_->file2_axpy(&LIA, &newLIA, 1, 0);
    /*dpd_file2_print(&newLIA,outfile);*/
    global_dpd_->file2_close(&LIA);

    global_dpd_->file2_copy(&newLIA, PSIF_CC_LAMBDA, "New Lia");  /* spin-adaptation for RHF */
    global_dpd_->file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/

    /* newLia * Dia */
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_copy(&newLIA, PSIF_CC_LAMBDA, "New LIA Increment");
    global_dpd_->file2_close(&newLIA);

    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA Increment");
    global_dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
    global_dpd_->file2_dirprd(&dIA, &newLIA);
    global_dpd_->file2_close(&dIA);
    global_dpd_->file2_close(&newLIA);

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_copy(&LIA, PSIF_CC_LAMBDA, "New LIA");
    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA Increment");
    global_dpd_->file2_axpy(&LIA, &newLIA, 1, 0);
    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&newLIA);

    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    global_dpd_->file2_copy(&newLia, PSIF_CC_LAMBDA, "New Lia Increment");
    global_dpd_->file2_close(&newLia);

    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia Increment");
    global_dpd_->file2_init(&dia, PSIF_CC_DENOM, L_irr, 0, 1, "dia");
    global_dpd_->file2_dirprd(&dia, &newLia);
    global_dpd_->file2_close(&dia);
    global_dpd_->file2_close(&newLia);

    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    global_dpd_->file2_copy(&Lia, PSIF_CC_LAMBDA, "New Lia");
    global_dpd_->file2_close(&Lia);
    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia Increment");
    global_dpd_->file2_axpy(&Lia, &newLia, 1, 0);
    global_dpd_->file2_close(&Lia);
    global_dpd_->file2_close(&newLia);
  }
  else if(params.ref == 2) { /** UHF **/

    /* newLia * Dia */
    global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
    global_dpd_->file2_dirprd(&dIA, &newLIA);
    global_dpd_->file2_close(&dIA);
    global_dpd_->file2_close(&newLIA);

    global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");
    global_dpd_->file2_init(&dia, PSIF_CC_DENOM, L_irr, 2, 3, "dia");
    global_dpd_->file2_dirprd(&dia, &newLia);
    global_dpd_->file2_close(&dia);
    global_dpd_->file2_close(&newLia);
  }

#ifdef EOM_DEBUG
  check_sum("after L1 build",L_irr);
#endif

  return;
}


}} // namespace psi::cclambda
