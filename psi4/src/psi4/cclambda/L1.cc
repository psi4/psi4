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
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

    void local_filter_T1(dpdfile2 *);

    void L1_build(struct L_Params L_params) {
      dpdfile2 newLIA, newLia, LIA, Lia;
      dpdfile2 dIA, dia, Fme, FME;
      dpdfile2 LFaet2, LFAEt2, LFmit2, LFMIt2;
      dpdfile2 GMI, Gmi, Gae, XIA, Xia;
      dpdfile2 GAE;
      dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ;
      dpdbuf4 WMBIJ, Wmbij, WMbIj, WmBiJ;
      dpdbuf4 LIJAB, Lijab, LIjAb, LiJaB, L2;
      dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE;
      dpdbuf4 WAMEF, Wamef, WAmEf, WaMeF, W;
      dpdbuf4 Z, D;
      dpdfile2 XLD;
      int Gim, Gi, Gm, Ga, Gam, nrows, ncols, A, a, am;
      int Gei, ei, e, i, Gef, Ge, Gf, E, I, af, fa, f;
      double *X;
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
	// Add T3 contribution to CCSD(T) lambda equations  
	if(params.wfn == "CCSD_T") {
	    global_dpd_->file2_init(&FME,PSIF_CC_OEI, 0, 0, 1, "SIA");
	    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
	    global_dpd_->file2_axpy(&FME, &LIA, 1, 0);
	    global_dpd_->file2_close(&LIA);
	    global_dpd_->file2_close(&FME);

	    global_dpd_->file2_init(&FME,PSIF_CC_OEI, 0, 2, 3, "Sia");
	    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");
	    global_dpd_->file2_axpy(&FME, &LIA, 1, 0);
	    global_dpd_->file2_close(&LIA);
	    global_dpd_->file2_close(&FME);
	  }
	}
      }
      /* excited state - no inhomogenous term, first term is -energy*L*/
      else if (!params.zeta) {
	if (params.ref == 0 || params.ref == 1) {
	  global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
	  global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
	  global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
	  global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
	  global_dpd_->file2_axpy(&LIA, &newLIA, -1.0 * L_params.cceom_energy,0.0);
	  global_dpd_->file2_axpy(&Lia, &newLia, -1.0 * L_params.cceom_energy,0.0);
	  global_dpd_->file2_close(&LIA);
	  global_dpd_->file2_close(&newLIA);
	  global_dpd_->file2_close(&Lia);
	  global_dpd_->file2_close(&newLia);
	}
	else if (params.ref == 2) {
	  /* do nothing - TDC did not change to increments for the UHF case */
	}
      }
      /* solving zeta equations; inhomogeneous term is Xi */
      else {
	if (params.ref == 0) {
	  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, 0, 0, 1, "XIA");
	  global_dpd_->file2_copy(&XIA, PSIF_CC_LAMBDA, "New LIA");
	  global_dpd_->file2_close(&XIA);
	}
	else if (params.ref == 1) {
	  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, 0, 0, 1, "XIA");
	  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, 0, 0, 1, "Xia");
	  global_dpd_->file2_copy(&XIA, PSIF_CC_LAMBDA, "New LIA");
	  global_dpd_->file2_copy(&Xia, PSIF_CC_LAMBDA, "New Lia");
	  global_dpd_->file2_close(&XIA);
	  global_dpd_->file2_close(&Xia);
	}
	else if(params.ref == 2) {
	  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, 0, 0, 1, "XIA");
	  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, 0, 2, 3, "Xia");
	  global_dpd_->file2_copy(&XIA, PSIF_CC_LAMBDA, "New LIA");
	  global_dpd_->file2_copy(&Xia, PSIF_CC_LAMBDA, "New Lia");
	  global_dpd_->file2_close(&XIA);
	  global_dpd_->file2_close(&Xia);
	}
      }

      if(params.ref == 0 || params.ref == 1) {
	global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
	global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
      }
      else if(params.ref == 2) {
	global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
	global_dpd_->file2_init(&newLia, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");
      }

      if(params.ref == 0) { /** RHF **/
	global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");

	/* L1 RHS += Lie*Fea */
	global_dpd_->file2_init(&LFAEt2, PSIF_CC_OEI, 0, 1, 1, "FAE");
	global_dpd_->contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&LFAEt2);

	/* L1 RHS += -Lma*Fim */
	global_dpd_->file2_init(&LFMIt2,PSIF_CC_OEI, 0, 0, 0, "FMI");
	global_dpd_->contract222(&LFMIt2,&LIA,&newLIA, 0, 1, -1.0, 1.0);
	global_dpd_->file2_close(&LFMIt2);

	/* L1 RHS += Lme*Wieam */
	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
	global_dpd_->contract422(&W, &LIA, &newLIA, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&W);

	global_dpd_->file2_close(&LIA);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* L1 RHS += Lie*Fea */
	global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
	global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");

	global_dpd_->file2_init(&LFAEt2, PSIF_CC_OEI, 0, 1, 1, "FAE");
	global_dpd_->file2_init(&LFaet2, PSIF_CC_OEI, 0, 1, 1, "Fae");
	global_dpd_->contract222(&Lia,&LFaet2,&newLia, 0, 1, 1.0, 1.0);
	global_dpd_->contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&LFaet2);
	global_dpd_->file2_close(&LFAEt2);

	/* L1 RHS += -Lma*Fim */
	global_dpd_->file2_init(&LFMIt2,PSIF_CC_OEI, 0, 0, 0, "FMI");
	global_dpd_->file2_init(&LFmit2,PSIF_CC_OEI, 0, 0, 0, "Fmi");
	global_dpd_->contract222(&LFmit2,&Lia,&newLia, 0, 1, -1.0, 1.0);
	global_dpd_->contract222(&LFMIt2,&LIA,&newLIA, 0, 1, -1.0, 1.0);
	global_dpd_->file2_close(&LFmit2);
	global_dpd_->file2_close(&LFMIt2);

	/* L1 RHS += Lme*Wieam */
	global_dpd_->buf4_init(&WMBEJ, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
	global_dpd_->contract422(&WMBEJ, &LIA, &newLIA, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&WMBEJ);

	global_dpd_->buf4_init(&WMbEj, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
	global_dpd_->contract422(&WMbEj, &Lia, &newLIA, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&WMbEj);

	global_dpd_->buf4_init(&Wmbej, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
	global_dpd_->contract422(&Wmbej, &Lia, &newLia, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Wmbej);

	global_dpd_->buf4_init(&WmBeJ, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
	global_dpd_->contract422(&WmBeJ, &LIA, &newLia, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&WmBeJ);

	global_dpd_->file2_close(&LIA);
	global_dpd_->file2_close(&Lia);
      }
      else if(params.ref == 2) { /** UHF **/

	/* L1 RHS += Lie*Fea */
	global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
	global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");

	global_dpd_->file2_init(&LFAEt2, PSIF_CC_OEI, 0, 1, 1, "FAEt");
	global_dpd_->file2_init(&LFaet2, PSIF_CC_OEI, 0, 3, 3, "Faet");
	global_dpd_->contract222(&Lia,&LFaet2,&newLia, 0, 1, 1, 1);
	global_dpd_->contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1, 1);
	global_dpd_->file2_close(&LFaet2);
	global_dpd_->file2_close(&LFAEt2);

	/* L1 RHS += -Lma*Fim */
	global_dpd_->file2_init(&LFMIt2,PSIF_CC_OEI, 0, 0, 0, "FMIt");
	global_dpd_->file2_init(&LFmit2,PSIF_CC_OEI, 0, 2, 2, "Fmit");
	global_dpd_->contract222(&LFmit2,&Lia,&newLia, 0, 1, -1, 1);
	global_dpd_->contract222(&LFMIt2,&LIA,&newLIA, 0, 1, -1, 1);
	global_dpd_->file2_close(&LFmit2);
	global_dpd_->file2_close(&LFMIt2);

	/* L1 RHS += Lme*Wieam */
	global_dpd_->buf4_init(&WMBEJ, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
	global_dpd_->contract422(&WMBEJ, &LIA, &newLIA, 0, 0, 1, 1);
	global_dpd_->buf4_close(&WMBEJ);

	global_dpd_->buf4_init(&WMbEj, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
	global_dpd_->contract422(&WMbEj, &Lia, &newLIA, 0, 0, 1, 1);
	global_dpd_->buf4_close(&WMbEj);

	global_dpd_->buf4_init(&Wmbej, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
	global_dpd_->contract422(&Wmbej, &Lia, &newLia, 0, 0, 1, 1);
	global_dpd_->buf4_close(&Wmbej);

	global_dpd_->buf4_init(&WmBeJ, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
	global_dpd_->contract422(&WmBeJ, &LIA, &newLia, 0, 0, 1, 1);
	global_dpd_->buf4_close(&WmBeJ);

	global_dpd_->file2_close(&LIA);
	global_dpd_->file2_close(&Lia);
      }

      /* L1 RHS += 1/2 Limef*Wefam */
      /* L(i,a) += [ 2 L(im,ef) - L(im,fe) ] * W(am,ef) */
      /* Note: W(am,ef) is really Wabei (ei,ab) */
      if(params.ref == 0) { /** RHF **/

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
	/*    dpd_contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0); */
	global_dpd_->file2_mat_init(&newLIA);
	global_dpd_->file2_mat_rd(&newLIA);
	for(Gam=0; Gam < moinfo.nirreps; Gam++) {
	  Gef = Gam; /* W is totally symmetric */
	  Gim = Gef ^ L_irr;
	  global_dpd_->buf4_mat_irrep_init(&L2, Gim);
	  global_dpd_->buf4_mat_irrep_rd(&L2, Gim);
	  global_dpd_->buf4_mat_irrep_shift13(&L2, Gim);

	  for(Gi=0; Gi < moinfo.nirreps; Gi++) {
	    Ga = Gi ^ L_irr;
	    Gm = Ga ^ Gam;
	    W.matrix[Gam] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gm],W.params->coltot[Gam]);

	    nrows = moinfo.occpi[Gi];
	    ncols = moinfo.occpi[Gm] * W.params->coltot[Gam];

	    for(A=0; A < moinfo.virtpi[Ga]; A++) {
	      a = moinfo.vir_off[Ga] + A;
	      am = W.row_offset[Gam][a];

	      global_dpd_->buf4_mat_irrep_rd_block(&W, Gam, am, moinfo.occpi[Gm]);

	      if(nrows && ncols && moinfo.virtpi[Ga])
		C_DGEMV('n',nrows,ncols,1,L2.shift.matrix[Gim][Gi][0],ncols,W.matrix[Gam][0],1,
			1, &(newLIA.matrix[Gi][0][A]), moinfo.virtpi[Ga]);

	    }

	    global_dpd_->free_dpd_block(W.matrix[Gam], moinfo.occpi[Gm], W.params->coltot[Gam]);
	  }
	  global_dpd_->buf4_mat_irrep_close(&L2, Gim);
	}
	global_dpd_->file2_mat_wrt(&newLIA);
	global_dpd_->file2_mat_close(&newLIA);
	global_dpd_->buf4_close(&L2);
	global_dpd_->buf4_close(&W);

      }
      else if(params.ref == 1) { /** ROHF **/

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
	global_dpd_->contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&W);
	global_dpd_->buf4_close(&L2);
	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	global_dpd_->contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&W);
	global_dpd_->buf4_close(&L2);

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "Lijab");
	global_dpd_->contract442(&L2, &W, &newLia, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&W);
	global_dpd_->buf4_close(&L2);
	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
	global_dpd_->contract442(&L2, &W, &newLia, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&W);
	global_dpd_->buf4_close(&L2);
      }
      else if(params.ref == 2) {

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
	global_dpd_->contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
	global_dpd_->buf4_close(&W);
	global_dpd_->buf4_close(&L2);
	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
	global_dpd_->contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
	global_dpd_->buf4_close(&W);
	global_dpd_->buf4_close(&L2);

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "Lijab");
	global_dpd_->contract442(&L2, &W, &newLia, 0, 0, 1, 1);
	global_dpd_->buf4_close(&W);
	global_dpd_->buf4_close(&L2);
	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
	global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
	global_dpd_->contract442(&L2, &W, &newLia, 0, 0, 1, 1);
	global_dpd_->buf4_close(&W);
	global_dpd_->buf4_close(&L2);

      }

      /* L1 RHS += -1/2 Lmnae*Wiemn */
      if(params.ref == 0) {
	global_dpd_->buf4_init(&WMbIj, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
	global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
	global_dpd_->contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1.0, 1.0);
	global_dpd_->buf4_close(&LIjAb);
	global_dpd_->buf4_close(&WMbIj);
      }
      else if(params.ref == 1) {

	global_dpd_->buf4_init(&WMBIJ, PSIF_CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
	global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
	global_dpd_->contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1.0, 1.0);
	global_dpd_->buf4_close(&LIJAB);
	global_dpd_->buf4_close(&WMBIJ);

	global_dpd_->buf4_init(&WMbIj, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
	global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	global_dpd_->contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1.0, 1.0);
	global_dpd_->buf4_close(&LIjAb);
	global_dpd_->buf4_close(&WMbIj);

	global_dpd_->buf4_init(&Wmbij, PSIF_CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
	global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "Lijab");
	global_dpd_->contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1.0, 1.0);
	global_dpd_->buf4_close(&Lijab);
	global_dpd_->buf4_close(&Wmbij);

	global_dpd_->buf4_init(&WmBiJ, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
	global_dpd_->buf4_init(&LiJaB, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
	global_dpd_->contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1.0, 1.0);
	global_dpd_->buf4_close(&LiJaB);
	global_dpd_->buf4_close(&WmBiJ);
      }
      else if(params.ref == 2) {

	global_dpd_->buf4_init(&WMBIJ, PSIF_CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
	global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
	global_dpd_->contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1, 1);
	global_dpd_->buf4_close(&LIJAB);
	global_dpd_->buf4_close(&WMBIJ);

	global_dpd_->buf4_init(&WMbIj, PSIF_CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
	global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
	global_dpd_->contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1, 1);
	global_dpd_->buf4_close(&LIjAb);
	global_dpd_->buf4_close(&WMbIj);

	global_dpd_->buf4_init(&Wmbij, PSIF_CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
	global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "Lijab");
	global_dpd_->contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1, 1);
	global_dpd_->buf4_close(&Lijab);
	global_dpd_->buf4_close(&Wmbij);

	global_dpd_->buf4_init(&WmBiJ, PSIF_CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
	global_dpd_->buf4_init(&LiJaB, PSIF_CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
	global_dpd_->contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1, 1);
	global_dpd_->buf4_close(&LiJaB);
	global_dpd_->buf4_close(&WmBiJ);
      }


      /* L1 RHS += -Gef*Weifa */
      if(params.ref == 0) {

	/*     dpd_file2_init(&GAE, CC_LAMBDA, L_irr, 1, 1, "GAE"); */
	/*     dpd_buf4_init(&WaMeF, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); */
	/*     dpd_dot13(&GAE,&WaMeF,&newLIA, 0, 0, -1.0, 1.0); */
	/*     dpd_buf4_close(&WaMeF); */
	/*     dpd_file2_close(&GAE); */

	/* Above code replaced to remove disk-space and memory bottlenecks 7/26/05, -TDC */
	global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");
	global_dpd_->file2_mat_init(&GAE);
	global_dpd_->file2_mat_rd(&GAE);
	global_dpd_->file2_mat_init(&newLIA);
	global_dpd_->file2_mat_rd(&newLIA);
	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
	for(Gei=0; Gei < moinfo.nirreps; Gei++) {
	  global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
	  X = init_array(W.params->coltot[Gei]);
	  for(ei=0; ei < W.params->rowtot[Gei]; ei++) {
	    global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
	    e = W.params->roworb[Gei][ei][0];
	    i = W.params->roworb[Gei][ei][1];
	    Ge = W.params->psym[e];
	    Gf = Ge ^ L_irr;
	    Gi = Ge ^ Gei;
	    Ga = Gi ^ L_irr;
	    E = e - moinfo.vir_off[Ge];
	    I = i - moinfo.occ_off[Gi];

	    zero_arr(X,W.params->coltot[Gei]);

	    for(fa=0; fa < W.params->coltot[Gei]; fa++) {
	      f = W.params->colorb[Gei][fa][0];
	      a = W.params->colorb[Gei][fa][1];
	      af = W.params->colidx[a][f];
	      X[fa] = 2.0 * W.matrix[Gei][0][fa] - W.matrix[Gei][0][af];
	    }

	    nrows = moinfo.virtpi[Gf];
	    ncols = moinfo.virtpi[Ga];
	    if(nrows && ncols)
	      C_DGEMV('t',nrows,ncols,-1,&X[W.col_offset[Gei][Gf]],ncols,
		      GAE.matrix[Ge][E],1,1,newLIA.matrix[Gi][I],1);

	  }
	  global_dpd_->buf4_mat_irrep_row_close(&W, Gei);
	  free(X);
	}
	global_dpd_->buf4_close(&W);
	global_dpd_->file2_mat_wrt(&newLIA);
	global_dpd_->file2_mat_close(&newLIA);
	global_dpd_->file2_mat_close(&GAE);
	global_dpd_->file2_close(&GAE);
      }
      else if(params.ref == 1) {

	global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");
	global_dpd_->file2_init(&Gae, PSIF_CC_LAMBDA, L_irr, 1, 1, "Gae");

	global_dpd_->buf4_init(&WAMEF, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
	global_dpd_->dot13(&GAE,&WAMEF,&newLIA, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&WAMEF);

	global_dpd_->buf4_init(&WaMeF, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
	global_dpd_->dot13(&Gae,&WaMeF,&newLIA, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&WaMeF);

	global_dpd_->buf4_init(&Wamef, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
	global_dpd_->dot13(&Gae,&Wamef,&newLia, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&Wamef);

	global_dpd_->buf4_init(&WAmEf, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
	global_dpd_->dot13(&GAE,&WAmEf,&newLia, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&WAmEf);

	/*
	  dpd_buf4_init(&WAMEF, CC_HBAR, 0, 10, 5, 10, 7, 0, "WAMEF");
	  dpd_dot23(&GAE,&WAMEF,&newLIA, 0, 0, -1.0, 1.0);
	  dpd_buf4_close(&WAMEF);
	  dpd_buf4_init(&WaMeF, CC_HBAR, 0, 10, 5, 10, 5, 0, "WaMeF");
	  dpd_dot23(&Gae,&WaMeF,&newLIA, 0, 0, -1.0, 1.0);
	  dpd_buf4_close(&WaMeF);
	  dpd_buf4_init(&Wamef, CC_HBAR, 0, 10, 5, 10, 7, 0, "Wamef");
	  dpd_dot23(&Gae,&Wamef,&newLia, 0, 0, -1.0, 1.0);
	  dpd_buf4_close(&Wamef);
	  dpd_buf4_init(&WAmEf, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
	  dpd_dot23(&GAE,&WAmEf,&newLia, 0, 0, -1.0, 1.0);
	  dpd_buf4_close(&WAmEf);
	*/

	global_dpd_->file2_close(&Gae);
	global_dpd_->file2_close(&GAE);
      }
      else if(params.ref == 2) {

	global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");
	global_dpd_->file2_init(&Gae, PSIF_CC_LAMBDA, L_irr, 3, 3, "Gae");

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 5, 21, 7, 0, "WAMEF");
	global_dpd_->dot13(&GAE,&W,&newLIA, 0, 0, -1, 1);
	global_dpd_->buf4_close(&W);

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
	global_dpd_->dot13(&Gae,&W,&newLIA, 0, 0, -1, 1);
	global_dpd_->buf4_close(&W);

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 31, 15, 31, 17, 0, "Wamef");
	global_dpd_->dot13(&Gae,&W,&newLia, 0, 0, -1, 1);
	global_dpd_->buf4_close(&W);

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
	global_dpd_->dot13(&GAE,&W,&newLia, 0, 0, -1, 1);
	global_dpd_->buf4_close(&W);

	global_dpd_->file2_close(&Gae);
	global_dpd_->file2_close(&GAE);

      }

      /* L1 RHS += -Gmn*Wmina */
      if(params.ref == 0) {
	global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");

	global_dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
	global_dpd_->dot14(&GMI, &WmNiE, &newLIA, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&WmNiE);

	global_dpd_->file2_close(&GMI);
      }
      else if(params.ref == 1) {

	global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");
	global_dpd_->file2_init(&Gmi, PSIF_CC_LAMBDA, L_irr, 0, 0, "Gmi");

	global_dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
	global_dpd_->dot14(&GMI, &WMNIE, &newLIA, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&WMNIE);

	global_dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
	global_dpd_->dot14(&Gmi, &WmNiE, &newLIA, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&WmNiE);

	global_dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
	global_dpd_->dot14(&Gmi, &Wmnie, &newLia, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&Wmnie);

	global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
	global_dpd_->dot14(&GMI, &WMnIe, &newLia, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&WMnIe);

	global_dpd_->file2_close(&Gmi);
	global_dpd_->file2_close(&GMI);

      }
      else if(params.ref == 2) {

	global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");
	global_dpd_->file2_init(&Gmi, PSIF_CC_LAMBDA, L_irr, 2, 2, "Gmi");

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 21, 2, 21, 0, "WMNIE (M>N,EI)");
	global_dpd_->dot14(&GMI, &W, &newLIA, 0, 0, -1, 1);
	global_dpd_->buf4_close(&W);

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
	global_dpd_->dot14(&Gmi, &W, &newLIA, 0, 0, -1, 1);
	global_dpd_->buf4_close(&W);

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 31, 12, 31, 0, "Wmnie (m>n,ei)");
	global_dpd_->dot14(&Gmi, &W, &newLia, 0, 0, -1, 1);
	global_dpd_->buf4_close(&W);

	global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
	global_dpd_->dot14(&GMI, &W, &newLia, 0, 0, -1, 1);
	global_dpd_->buf4_close(&W);

	global_dpd_->file2_close(&Gmi);
	global_dpd_->file2_close(&GMI);
      }

      /* CC3 T3->L1 */
      if(params.wfn == "CC3") {
	if(params.ref == 0) {

	  global_dpd_->file2_init(&XLD, PSIF_CC3_MISC, 0, 0, 1, "CC3 XLD");
	  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
	  global_dpd_->dot24(&XLD, &D, &newLIA, 0, 0, 1, 1);
	  global_dpd_->buf4_close(&D);
	  global_dpd_->file2_close(&XLD);

	  global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
	  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 10, 0, 10, 0, 0, "CC3 ZIFLN");
	  global_dpd_->contract442(&Z, &L2, &newLIA, 0, 2, -0.5, 1);
	  global_dpd_->buf4_close(&Z);
	  global_dpd_->buf4_close(&L2);

	  global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 10, 0, 10, 0, 0, "CC3 ZIfLn");
	  global_dpd_->contract442(&Z, &L2, &newLIA, 0, 2, -1, 1);
	  global_dpd_->buf4_close(&Z);
	  global_dpd_->buf4_close(&L2);

	  global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
	  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 11, 5, 11, 5, 0, "CC3 ZDFAN (AN,DF)");
	  global_dpd_->contract442(&L2, &Z, &newLIA, 0, 0, 0.5, 1);
	  global_dpd_->buf4_close(&Z);
	  global_dpd_->buf4_close(&L2);

	  global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 11, 5, 11, 5, 0, "CC3 ZDfAn (An,Df)");
	  global_dpd_->contract442(&L2, &Z, &newLIA, 0, 0, 1.0, 1);
	  global_dpd_->buf4_close(&Z);
	  global_dpd_->buf4_close(&L2);
	}
      }

      global_dpd_->file2_close(&newLIA);
      global_dpd_->file2_close(&newLia);

      /* newLia * Dia */
      if(params.ref == 0) { /** RHF **/

        if(params.wfn == "CCSD_T") {
            global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "SIA(T)");
            global_dpd_->file2_init(&newLIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
            global_dpd_->file2_axpy(&FME, &newLIA, 1, 0);
            global_dpd_->file2_close(&newLIA);
            global_dpd_->file2_close(&FME);
          }

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
	global_dpd_->file2_close(&LIA);

	global_dpd_->file2_copy(&newLIA, PSIF_CC_LAMBDA, "New Lia");  /* spin-adaptation for RHF */
	global_dpd_->file2_close(&newLIA);
      }
      else if(params.ref == 1) { /** ROHF **/

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
      else if(params.ref == 2) {

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
