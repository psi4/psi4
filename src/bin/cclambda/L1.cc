/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
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
	  dpd_file2_init(&FME,CC_OEI, 0, 0, 1, "FME");
	  dpd_file2_copy(&FME, CC_LAMBDA, "New LIA");
	  dpd_file2_close(&FME);
	}
	else if(params.ref == 1) {
	  dpd_file2_init(&Fme,CC_OEI, 0, 0, 1, "Fme");
	  dpd_file2_init(&FME,CC_OEI, 0, 0, 1, "FME");
	  dpd_file2_copy(&Fme, CC_LAMBDA, "New Lia");
	  dpd_file2_copy(&FME, CC_LAMBDA, "New LIA");
	  dpd_file2_close(&Fme);
	  dpd_file2_close(&FME);
	}
	else if(params.ref == 2) {
	  dpd_file2_init(&Fme,CC_OEI, 0, 2, 3, "Fme");
	  dpd_file2_init(&FME,CC_OEI, 0, 0, 1, "FME");
	  dpd_file2_copy(&Fme, CC_LAMBDA, "New Lia");
	  dpd_file2_copy(&FME, CC_LAMBDA, "New LIA");
	  dpd_file2_close(&Fme);
	  dpd_file2_close(&FME);
	  /* If CCSD(T) gradient, add T3 contributions */
	  if(!strcmp(params.wfn, "CCSD_T")) {
	    dpd_file2_init(&FME,CC_OEI, 0, 0, 1, "SIA");
	    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	    dpd_file2_axpy(&FME, &LIA, 1, 0);
	    dpd_file2_close(&LIA);
	    dpd_file2_close(&FME);

	    dpd_file2_init(&FME,CC_OEI, 0, 2, 3, "Sia");
	    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 2, 3, "New Lia");
	    dpd_file2_axpy(&FME, &LIA, 1, 0);
	    dpd_file2_close(&LIA);
	    dpd_file2_close(&FME);
	  }
	}
      }
      /* excited state - no inhomogenous term, first term is -energy*L*/
      else if (!params.zeta) {
	if (params.ref == 0 || params.ref == 1) {
	  dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
	  dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	  dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");
	  dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 0, 1, "New Lia");
	  dpd_file2_axpy(&LIA, &newLIA, -1.0 * L_params.cceom_energy,0.0);
	  dpd_file2_axpy(&Lia, &newLia, -1.0 * L_params.cceom_energy,0.0);
	  dpd_file2_close(&LIA);
	  dpd_file2_close(&newLIA);
	  dpd_file2_close(&Lia);
	  dpd_file2_close(&newLia);
	}
	else if (params.ref == 2) {
	  /* do nothing - TDC did not change to increments for the UHF case */
	}
      }
      /* solving zeta equations; inhomogeneous term is Xi */
      else {
	if (params.ref == 0) {
	  dpd_file2_init(&XIA, EOM_XI, 0, 0, 1, "XIA");
	  dpd_file2_copy(&XIA, CC_LAMBDA, "New LIA");
	  dpd_file2_close(&XIA);
	}
	else if (params.ref == 1) {
	  dpd_file2_init(&XIA, EOM_XI, 0, 0, 1, "XIA");
	  dpd_file2_init(&Xia, EOM_XI, 0, 0, 1, "Xia");
	  dpd_file2_copy(&XIA, CC_LAMBDA, "New LIA");
	  dpd_file2_copy(&Xia, CC_LAMBDA, "New Lia");
	  dpd_file2_close(&XIA);
	  dpd_file2_close(&Xia);
	}
	else if(params.ref == 2) {
	  dpd_file2_init(&XIA, EOM_XI, 0, 0, 1, "XIA");
	  dpd_file2_init(&Xia, EOM_XI, 0, 2, 3, "Xia");
	  dpd_file2_copy(&XIA, CC_LAMBDA, "New LIA");
	  dpd_file2_copy(&Xia, CC_LAMBDA, "New Lia");
	  dpd_file2_close(&XIA);
	  dpd_file2_close(&Xia);
	}
      }

      if(params.ref == 0 || params.ref == 1) {
	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 0, 1, "New Lia");
      }
      else if(params.ref == 2) {
	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 2, 3, "New Lia");
      }

      if(params.ref == 0) { /** RHF **/
	dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");

	/* L1 RHS += Lie*Fea */
	dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAE");
	dpd_contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1.0, 1.0);
	dpd_file2_close(&LFAEt2);

	/* L1 RHS += -Lma*Fim */
	dpd_file2_init(&LFMIt2,CC_OEI, 0, 0, 0, "FMI");
	dpd_contract222(&LFMIt2,&LIA,&newLIA, 0, 1, -1.0, 1.0);
	dpd_file2_close(&LFMIt2);

	/* L1 RHS += Lme*Wieam */
	dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
	dpd_contract422(&W, &LIA, &newLIA, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&W);

	dpd_file2_close(&LIA);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* L1 RHS += Lie*Fea */
	dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
	dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");

	dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAE");
	dpd_file2_init(&LFaet2, CC_OEI, 0, 1, 1, "Fae");
	dpd_contract222(&Lia,&LFaet2,&newLia, 0, 1, 1.0, 1.0);
	dpd_contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1.0, 1.0);
	dpd_file2_close(&LFaet2);
	dpd_file2_close(&LFAEt2);

	/* L1 RHS += -Lma*Fim */
	dpd_file2_init(&LFMIt2,CC_OEI, 0, 0, 0, "FMI");
	dpd_file2_init(&LFmit2,CC_OEI, 0, 0, 0, "Fmi");
	dpd_contract222(&LFmit2,&Lia,&newLia, 0, 1, -1.0, 1.0);
	dpd_contract222(&LFMIt2,&LIA,&newLIA, 0, 1, -1.0, 1.0);
	dpd_file2_close(&LFmit2);
	dpd_file2_close(&LFMIt2);

	/* L1 RHS += Lme*Wieam */
	dpd_buf4_init(&WMBEJ, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
	dpd_contract422(&WMBEJ, &LIA, &newLIA, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&WMBEJ);

	dpd_buf4_init(&WMbEj, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
	dpd_contract422(&WMbEj, &Lia, &newLIA, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&WMbEj);

	dpd_buf4_init(&Wmbej, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
	dpd_contract422(&Wmbej, &Lia, &newLia, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&Wmbej);

	dpd_buf4_init(&WmBeJ, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
	dpd_contract422(&WmBeJ, &LIA, &newLia, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&WmBeJ);

	dpd_file2_close(&LIA);
	dpd_file2_close(&Lia);
      }
      else if(params.ref == 2) { /** UHF **/

	/* L1 RHS += Lie*Fea */
	dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
	dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");

	dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAEt");
	dpd_file2_init(&LFaet2, CC_OEI, 0, 3, 3, "Faet");
	dpd_contract222(&Lia,&LFaet2,&newLia, 0, 1, 1, 1);
	dpd_contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1, 1);
	dpd_file2_close(&LFaet2);
	dpd_file2_close(&LFAEt2);

	/* L1 RHS += -Lma*Fim */
	dpd_file2_init(&LFMIt2,CC_OEI, 0, 0, 0, "FMIt");
	dpd_file2_init(&LFmit2,CC_OEI, 0, 2, 2, "Fmit");
	dpd_contract222(&LFmit2,&Lia,&newLia, 0, 1, -1, 1);
	dpd_contract222(&LFMIt2,&LIA,&newLIA, 0, 1, -1, 1);
	dpd_file2_close(&LFmit2);
	dpd_file2_close(&LFMIt2);

	/* L1 RHS += Lme*Wieam */
	dpd_buf4_init(&WMBEJ, CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
	dpd_contract422(&WMBEJ, &LIA, &newLIA, 0, 0, 1, 1);
	dpd_buf4_close(&WMBEJ);

	dpd_buf4_init(&WMbEj, CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
	dpd_contract422(&WMbEj, &Lia, &newLIA, 0, 0, 1, 1);
	dpd_buf4_close(&WMbEj);

	dpd_buf4_init(&Wmbej, CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
	dpd_contract422(&Wmbej, &Lia, &newLia, 0, 0, 1, 1);
	dpd_buf4_close(&Wmbej);

	dpd_buf4_init(&WmBeJ, CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
	dpd_contract422(&WmBeJ, &LIA, &newLia, 0, 0, 1, 1);
	dpd_buf4_close(&WmBeJ);

	dpd_file2_close(&LIA);
	dpd_file2_close(&Lia);
      }

      /* L1 RHS += 1/2 Limef*Wefam */
      /* L(i,a) += [ 2 L(im,ef) - L(im,fe) ] * W(am,ef) */
      /* Note: W(am,ef) is really Wabei (ei,ab) */
      if(params.ref == 0) { /** RHF **/

	dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
	/*    dpd_contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0); */
	dpd_file2_mat_init(&newLIA);
	dpd_file2_mat_rd(&newLIA);
	for(Gam=0; Gam < moinfo.nirreps; Gam++) {
	  Gef = Gam; /* W is totally symmetric */
	  Gim = Gef ^ L_irr; 
	  dpd_buf4_mat_irrep_init(&L2, Gim);
	  dpd_buf4_mat_irrep_rd(&L2, Gim);
	  dpd_buf4_mat_irrep_shift13(&L2, Gim);

	  for(Gi=0; Gi < moinfo.nirreps; Gi++) {
	    Ga = Gi ^ L_irr;
	    Gm = Ga ^ Gam;
	    W.matrix[Gam] = dpd_block_matrix(moinfo.occpi[Gm],W.params->coltot[Gam]);

	    nrows = moinfo.occpi[Gi];
	    ncols = moinfo.occpi[Gm] * W.params->coltot[Gam];

	    for(A=0; A < moinfo.virtpi[Ga]; A++) {
	      a = moinfo.vir_off[Ga] + A;
	      am = W.row_offset[Gam][a];

	      dpd_buf4_mat_irrep_rd_block(&W, Gam, am, moinfo.occpi[Gm]);

	      if(nrows && ncols && moinfo.virtpi[Ga]) 
		C_DGEMV('n',nrows,ncols,1,L2.shift.matrix[Gim][Gi][0],ncols,W.matrix[Gam][0],1,
			1, &(newLIA.matrix[Gi][0][A]), moinfo.virtpi[Ga]);

	    }

	    dpd_free_block(W.matrix[Gam], moinfo.occpi[Gm], W.params->coltot[Gam]);
	  }
	  dpd_buf4_mat_irrep_close(&L2, Gim);
	}
	dpd_file2_mat_wrt(&newLIA);
	dpd_file2_mat_close(&newLIA);
	dpd_buf4_close(&L2);
	dpd_buf4_close(&W);

      }
      else if(params.ref == 1) { /** ROHF **/

	dpd_buf4_init(&W, CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
	dpd_contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&W);
	dpd_buf4_close(&L2);
	dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	dpd_contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&W);
	dpd_buf4_close(&L2);

	dpd_buf4_init(&W, CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "Lijab");
	dpd_contract442(&L2, &W, &newLia, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&W);
	dpd_buf4_close(&L2);
	dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
	dpd_contract442(&L2, &W, &newLia, 0, 0, 1.0, 1.0);
	dpd_buf4_close(&W);
	dpd_buf4_close(&L2);
      }
      else if(params.ref == 2) {

	dpd_buf4_init(&W, CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
	dpd_contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
	dpd_buf4_close(&W);
	dpd_buf4_close(&L2);
	dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
	dpd_contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
	dpd_buf4_close(&W);
	dpd_buf4_close(&L2);

	dpd_buf4_init(&W, CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "Lijab");
	dpd_contract442(&L2, &W, &newLia, 0, 0, 1, 1);
	dpd_buf4_close(&W);
	dpd_buf4_close(&L2);
	dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
	dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
	dpd_contract442(&L2, &W, &newLia, 0, 0, 1, 1);
	dpd_buf4_close(&W);
	dpd_buf4_close(&L2);

      }

      /* L1 RHS += -1/2 Lmnae*Wiemn */
      if(params.ref == 0) {
	dpd_buf4_init(&WMbIj, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
	dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
	dpd_contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1.0, 1.0);
	dpd_buf4_close(&LIjAb);
	dpd_buf4_close(&WMbIj);
      }
      else if(params.ref == 1) {

	dpd_buf4_init(&WMBIJ, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
	dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
	dpd_contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1.0, 1.0);
	dpd_buf4_close(&LIJAB);
	dpd_buf4_close(&WMBIJ);

	dpd_buf4_init(&WMbIj, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
	dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	dpd_contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1.0, 1.0);
	dpd_buf4_close(&LIjAb);
	dpd_buf4_close(&WMbIj);

	dpd_buf4_init(&Wmbij, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
	dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "Lijab");
	dpd_contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1.0, 1.0);
	dpd_buf4_close(&Lijab);
	dpd_buf4_close(&Wmbij);

	dpd_buf4_init(&WmBiJ, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
	dpd_buf4_init(&LiJaB, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
	dpd_contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1.0, 1.0);
	dpd_buf4_close(&LiJaB);
	dpd_buf4_close(&WmBiJ);
      }
      else if(params.ref == 2) {

	dpd_buf4_init(&WMBIJ, CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
	dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
	dpd_contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1, 1);
	dpd_buf4_close(&LIJAB);
	dpd_buf4_close(&WMBIJ);

	dpd_buf4_init(&WMbIj, CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
	dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
	dpd_contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1, 1);
	dpd_buf4_close(&LIjAb);
	dpd_buf4_close(&WMbIj);

	dpd_buf4_init(&Wmbij, CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
	dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "Lijab");
	dpd_contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1, 1);
	dpd_buf4_close(&Lijab);
	dpd_buf4_close(&Wmbij);

	dpd_buf4_init(&WmBiJ, CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
	dpd_buf4_init(&LiJaB, CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
	dpd_contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1, 1);
	dpd_buf4_close(&LiJaB);
	dpd_buf4_close(&WmBiJ);
      }


      /* L1 RHS += -Gef*Weifa */
      if(params.ref == 0) {

	/*     dpd_file2_init(&GAE, CC_LAMBDA, L_irr, 1, 1, "GAE"); */
	/*     dpd_buf4_init(&WaMeF, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); */
	/*     dpd_dot13(&GAE,&WaMeF,&newLIA, 0, 0, -1.0, 1.0); */
	/*     dpd_buf4_close(&WaMeF); */
	/*     dpd_file2_close(&GAE); */

	/* Above code replaced to remove disk-space and memory bottlenecks 7/26/05, -TDC */
	dpd_file2_init(&GAE, CC_LAMBDA, L_irr, 1, 1, "GAE");
	dpd_file2_mat_init(&GAE);
	dpd_file2_mat_rd(&GAE);
	dpd_file2_mat_init(&newLIA);
	dpd_file2_mat_rd(&newLIA);
	dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
	for(Gei=0; Gei < moinfo.nirreps; Gei++) {
	  dpd_buf4_mat_irrep_row_init(&W, Gei);
	  X = init_array(W.params->coltot[Gei]);
	  for(ei=0; ei < W.params->rowtot[Gei]; ei++) {
	    dpd_buf4_mat_irrep_row_rd(&W, Gei, ei);
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
	  dpd_buf4_mat_irrep_row_close(&W, Gei);
	  free(X);
	}
	dpd_buf4_close(&W);
	dpd_file2_mat_wrt(&newLIA);
	dpd_file2_mat_close(&newLIA);
	dpd_file2_mat_close(&GAE);
	dpd_file2_close(&GAE);
      }
      else if(params.ref == 1) {

	dpd_file2_init(&GAE, CC_LAMBDA, L_irr, 1, 1, "GAE");
	dpd_file2_init(&Gae, CC_LAMBDA, L_irr, 1, 1, "Gae");

	dpd_buf4_init(&WAMEF, CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
	dpd_dot13(&GAE,&WAMEF,&newLIA, 0, 0, -1.0, 1.0);
	dpd_buf4_close(&WAMEF);

	dpd_buf4_init(&WaMeF, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
	dpd_dot13(&Gae,&WaMeF,&newLIA, 0, 0, -1.0, 1.0);
	dpd_buf4_close(&WaMeF);

	dpd_buf4_init(&Wamef, CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
	dpd_dot13(&Gae,&Wamef,&newLia, 0, 0, -1.0, 1.0);
	dpd_buf4_close(&Wamef);

	dpd_buf4_init(&WAmEf, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
	dpd_dot13(&GAE,&WAmEf,&newLia, 0, 0, -1.0, 1.0);
	dpd_buf4_close(&WAmEf);

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

	dpd_file2_close(&Gae);
	dpd_file2_close(&GAE);
      }
      else if(params.ref == 2) {

	dpd_file2_init(&GAE, CC_LAMBDA, L_irr, 1, 1, "GAE");
	dpd_file2_init(&Gae, CC_LAMBDA, L_irr, 3, 3, "Gae");

	dpd_buf4_init(&W, CC_HBAR, 0, 21, 5, 21, 7, 0, "WAMEF");
	dpd_dot13(&GAE,&W,&newLIA, 0, 0, -1, 1);
	dpd_buf4_close(&W);

	dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
	dpd_dot13(&Gae,&W,&newLIA, 0, 0, -1, 1);
	dpd_buf4_close(&W);

	dpd_buf4_init(&W, CC_HBAR, 0, 31, 15, 31, 17, 0, "Wamef");
	dpd_dot13(&Gae,&W,&newLia, 0, 0, -1, 1);
	dpd_buf4_close(&W);

	dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
	dpd_dot13(&GAE,&W,&newLia, 0, 0, -1, 1);
	dpd_buf4_close(&W);

	dpd_file2_close(&Gae);
	dpd_file2_close(&GAE);

      }

      /* L1 RHS += -Gmn*Wmina */
      if(params.ref == 0) {
	dpd_file2_init(&GMI, CC_LAMBDA, L_irr, 0, 0, "GMI");

	dpd_buf4_init(&WmNiE, CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
	dpd_dot14(&GMI, &WmNiE, &newLIA, 0, 0, -1.0, 1.0);
	dpd_buf4_close(&WmNiE);

	dpd_file2_close(&GMI);
      }
      else if(params.ref == 1) {

	dpd_file2_init(&GMI, CC_LAMBDA, L_irr, 0, 0, "GMI");
	dpd_file2_init(&Gmi, CC_LAMBDA, L_irr, 0, 0, "Gmi");

	dpd_buf4_init(&WMNIE, CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
	dpd_dot14(&GMI, &WMNIE, &newLIA, 0, 0, -1.0, 1.0); 
	dpd_buf4_close(&WMNIE);

	dpd_buf4_init(&WmNiE, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
	dpd_dot14(&Gmi, &WmNiE, &newLIA, 0, 0, -1.0, 1.0);
	dpd_buf4_close(&WmNiE);

	dpd_buf4_init(&Wmnie, CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
	dpd_dot14(&Gmi, &Wmnie, &newLia, 0, 0, -1.0, 1.0);
	dpd_buf4_close(&Wmnie);

	dpd_buf4_init(&WMnIe, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
	dpd_dot14(&GMI, &WMnIe, &newLia, 0, 0, -1.0, 1.0);
	dpd_buf4_close(&WMnIe);

	dpd_file2_close(&Gmi);
	dpd_file2_close(&GMI);

      }
      else if(params.ref == 2) {

	dpd_file2_init(&GMI, CC_LAMBDA, L_irr, 0, 0, "GMI");
	dpd_file2_init(&Gmi, CC_LAMBDA, L_irr, 2, 2, "Gmi");

	dpd_buf4_init(&W, CC_HBAR, 0, 0, 21, 2, 21, 0, "WMNIE (M>N,EI)");
	dpd_dot14(&GMI, &W, &newLIA, 0, 0, -1, 1); 
	dpd_buf4_close(&W);

	dpd_buf4_init(&W, CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
	dpd_dot14(&Gmi, &W, &newLIA, 0, 0, -1, 1);
	dpd_buf4_close(&W);

	dpd_buf4_init(&W, CC_HBAR, 0, 10, 31, 12, 31, 0, "Wmnie (m>n,ei)");
	dpd_dot14(&Gmi, &W, &newLia, 0, 0, -1, 1);
	dpd_buf4_close(&W);

	dpd_buf4_init(&W, CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
	dpd_dot14(&GMI, &W, &newLia, 0, 0, -1, 1);
	dpd_buf4_close(&W);

	dpd_file2_close(&Gmi);
	dpd_file2_close(&GMI);
      }

      /* CC3 T3->L1 */
      if(!strcmp(params.wfn, "CC3")) {
	if(params.ref == 0) { 

	  dpd_file2_init(&XLD, CC3_MISC, 0, 0, 1, "CC3 XLD");
	  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
	  dpd_dot24(&XLD, &D, &newLIA, 0, 0, 1, 1);
	  dpd_buf4_close(&D);
	  dpd_file2_close(&XLD);

	  dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
	  dpd_buf4_init(&Z, CC3_MISC, 0, 10, 0, 10, 0, 0, "CC3 ZIFLN");
	  dpd_contract442(&Z, &L2, &newLIA, 0, 2, -0.5, 1);
	  dpd_buf4_close(&Z);
	  dpd_buf4_close(&L2);

	  dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  dpd_buf4_init(&Z, CC3_MISC, 0, 10, 0, 10, 0, 0, "CC3 ZIfLn");
	  dpd_contract442(&Z, &L2, &newLIA, 0, 2, -1, 1);
	  dpd_buf4_close(&Z);
	  dpd_buf4_close(&L2);

	  dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
	  dpd_buf4_init(&Z, CC3_MISC, 0, 11, 5, 11, 5, 0, "CC3 ZDFAN (AN,DF)");
	  dpd_contract442(&L2, &Z, &newLIA, 0, 0, 0.5, 1);
	  dpd_buf4_close(&Z);
	  dpd_buf4_close(&L2);

	  dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  dpd_buf4_init(&Z, CC3_MISC, 0, 11, 5, 11, 5, 0, "CC3 ZDfAn (An,Df)");
	  dpd_contract442(&L2, &Z, &newLIA, 0, 0, 1.0, 1);
	  dpd_buf4_close(&Z);
	  dpd_buf4_close(&L2);
	}
      }

      dpd_file2_close(&newLIA);
      dpd_file2_close(&newLia);

      /* newLia * Dia */
      if(params.ref == 0) { /** RHF **/

	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	dpd_file2_copy(&newLIA, CC_LAMBDA, "New LIA Increment");
	dpd_file2_close(&newLIA);

	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA Increment");
	if(params.local && local.filter_singles) local_filter_T1(&newLIA);
	else {
	  dpd_file2_init(&dIA, CC_DENOM, L_irr, 0, 1, "dIA");
	  dpd_file2_dirprd(&dIA, &newLIA);
	  dpd_file2_close(&dIA);
	}
	dpd_file2_close(&newLIA);

	dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
	dpd_file2_copy(&LIA, CC_LAMBDA, "New LIA");
	dpd_file2_close(&LIA);
	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "New LIA Increment");
	dpd_file2_axpy(&LIA, &newLIA, 1, 0);
	dpd_file2_close(&LIA);

	dpd_file2_copy(&newLIA, CC_LAMBDA, "New Lia");  /* spin-adaptation for RHF */
	dpd_file2_close(&newLIA);
      }
      else if(params.ref == 1) { /** ROHF **/

	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	dpd_file2_copy(&newLIA, CC_LAMBDA, "New LIA Increment");
	dpd_file2_close(&newLIA);

	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA Increment");
	dpd_file2_init(&dIA, CC_DENOM, L_irr, 0, 1, "dIA");
	dpd_file2_dirprd(&dIA, &newLIA);
	dpd_file2_close(&dIA);
	dpd_file2_close(&newLIA);

	dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
	dpd_file2_copy(&LIA, CC_LAMBDA, "New LIA");
	dpd_file2_close(&LIA);
	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "New LIA Increment");
	dpd_file2_axpy(&LIA, &newLIA, 1, 0);
	dpd_file2_close(&LIA);
	dpd_file2_close(&newLIA);

	dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 0, 1, "New Lia");
	dpd_file2_copy(&newLia, CC_LAMBDA, "New Lia Increment");
	dpd_file2_close(&newLia);

	dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 0, 1, "New Lia Increment");
	dpd_file2_init(&dia, CC_DENOM, L_irr, 0, 1, "dia");
	dpd_file2_dirprd(&dia, &newLia);
	dpd_file2_close(&dia);
	dpd_file2_close(&newLia);

	dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");
	dpd_file2_copy(&Lia, CC_LAMBDA, "New Lia");
	dpd_file2_close(&Lia);
	dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 0, 1, "New Lia");
	dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "New Lia Increment");
	dpd_file2_axpy(&Lia, &newLia, 1, 0);
	dpd_file2_close(&Lia);
	dpd_file2_close(&newLia);
      }
      else if(params.ref == 2) {

	dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
	dpd_file2_init(&dIA, CC_DENOM, L_irr, 0, 1, "dIA");
	dpd_file2_dirprd(&dIA, &newLIA);
	dpd_file2_close(&dIA);
	dpd_file2_close(&newLIA);

	dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 2, 3, "New Lia");
	dpd_file2_init(&dia, CC_DENOM, L_irr, 2, 3, "dia");
	dpd_file2_dirprd(&dia, &newLia);
	dpd_file2_close(&dia);
	dpd_file2_close(&newLia);
      }

#ifdef EOM_DEBUG
      check_sum("after L1 build",L_irr);
#endif

      return;
    }



  }} // namespace psi::cclambda
