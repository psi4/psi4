/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
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
    }
  }
  /* excited state - no inhomogenous term, first term is -energy*L*/
  else if (!params.zeta) {
    if (params.ref == 0 || params.ref == 1) {
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");
      dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 0, 1, "New Lia");
      dpd_file2_axpy(&LIA, &newLIA, -1.0 * L_params.cceom_energy, 0);
      dpd_file2_axpy(&Lia, &newLia, -1.0 * L_params.cceom_energy, 0);
      dpd_file2_close(&LIA);
      dpd_file2_close(&newLIA);
      dpd_file2_close(&Lia);
      dpd_file2_close(&newLia);
    }
    else if (params.ref == 2) {
      /* do nothing - TDC did not change to increments for the UHF case */
    }
  }

  if(params.ref == 0) {
    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");

    /* L1 RHS += Lie*Fea */
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&LIA,&FAE,&newLIA, 0, 1, 1, 1);
    dpd_file2_close(&FAE);

    /* L1 RHS += -Lma*Fim */
    dpd_file2_init(&FMI,CC_OEI, 0, 0, 0, "FMI");
    dpd_contract222(&FMI,&LIA,&newLIA, 0, 1, -1, 1);
    dpd_file2_close(&FMI);

    dpd_file2_close(&LIA);
    dpd_file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 2, 3, "New Lia");

    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");

    /* L1 RHS += Lie*Fea */
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 3, 3, "Fae");
    dpd_contract222(&Lia,&Fae,&newLia, 0, 1, 1, 1);
    dpd_contract222(&LIA,&FAE,&newLIA, 0, 1, 1, 1);
    dpd_file2_close(&Fae);
    dpd_file2_close(&FAE);

    /* L1 RHS += -Lma*Fim */
    dpd_file2_init(&FMI,CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&Fmi,CC_OEI, 0, 2, 2, "Fmi");
    dpd_contract222(&Fmi,&Lia,&newLia, 0, 1, -1, 1);
    dpd_contract222(&FMI,&LIA,&newLIA, 0, 1, -1, 1);
    dpd_file2_close(&Fmi);
    dpd_file2_close(&FMI);

    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);

    dpd_file2_close(&newLIA);
    dpd_file2_close(&newLia);
  }

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");

    /* L1 RHS += Lme*Wieam */
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
    dpd_contract422(&W, &LIA, &newLIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    dpd_file2_close(&LIA);
    dpd_file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/
  
    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 2, 3, "New Lia");

    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");

    dpd_buf4_init(&WMBEJ, CC2_HET1, 0, 20, 21, 20, 21, 0, "CC2 WMBEJ (ME,BJ)");
    dpd_contract422(&WMBEJ, &LIA, &newLIA, 1, 0, 1, 1);
    dpd_buf4_close(&WMBEJ);

    dpd_buf4_init(&Wmbej, CC2_HET1, 0, 30, 31, 30, 31, 0, "CC2 Wmbej (me,bj)");
    dpd_contract422(&Wmbej, &Lia, &newLia, 1, 0, 1, 1);
    dpd_buf4_close(&Wmbej);

    dpd_buf4_init(&WMbEj, CC2_HET1, 0, 20, 31, 20, 31, 0, "CC2 WMbEj (ME,bj)");
    dpd_contract422(&WMbEj, &Lia, &newLIA, 1, 0, -1, 1);
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WmBeJ, CC2_HET1, 0, 30, 21, 30, 21, 0, "CC2 WmBeJ (me,BJ)");
    dpd_contract422(&WmBeJ, &LIA, &newLia, 1, 0, -1, 1);
    dpd_buf4_close(&WmBeJ);

    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);

    dpd_file2_close(&newLIA);
    dpd_file2_close(&newLia);
 
  }

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    dpd_buf4_sort(&L2, CC_TMP0, rspq, 5, 0, "Z (2 AbIj - AbjI)");
    dpd_buf4_close(&L2);

    /* L1 RHS += 1/2 Limef*Wefam */
    /* Out-of-core contract442 */

    dpd_buf4_init(&W, CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
    dpd_buf4_init(&L2, CC_TMP0, L_irr, 5, 0, 5, 0, 0, "Z (2 AbIj - AbjI)");

    /* dpd_contract442(&L2, &W, &newLIA, 0, 0, 1, 1); */

    GW = W.file.my_irrep;
    GL2 = L2.file.my_irrep;
    GL1 = newLIA.my_irrep;

    dpd_file2_mat_init(&newLIA);
    dpd_file2_mat_rd(&newLIA);

    for(Gab=0; Gab < moinfo.nirreps; Gab++) {

      dpd_buf4_mat_irrep_row_init(&L2, Gab);
      dpd_buf4_mat_irrep_row_init(&W, Gab);

      for(ab=0; ab < L2.params->rowtot[Gab]; ab++) {

	dpd_buf4_mat_irrep_row_zero(&L2, Gab, ab);
	dpd_buf4_mat_irrep_row_rd(&L2, Gab, ab);

	dpd_buf4_mat_irrep_row_zero(&W, Gab, ab);
	dpd_buf4_mat_irrep_row_rd(&W, Gab, ab);

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
      dpd_buf4_mat_irrep_row_close(&L2, Gab);
      dpd_buf4_mat_irrep_row_close(&W, Gab);
    }
    dpd_file2_mat_wrt(&newLIA);
    dpd_file2_mat_close(&newLIA);

    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);

    dpd_file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 2, 3, "New Lia");

    dpd_buf4_init(&W, CC2_HET1, 0, 21, 7, 21, 7, 0, "CC2 WABEI (EI,A>B)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&W, CC2_HET1, 0, 26, 28, 26, 28, 0, "CC2 WAbEi (Ei,Ab)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&W, CC2_HET1, 0, 31, 17, 31, 17, 0, "CC2 Wabei (ei,a>b)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "Lijab");
    dpd_contract442(&L2, &W, &newLia, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&W, CC2_HET1, 0, 25, 29, 25, 29, 0, "CC2 WaBeI (eI,aB)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&L2, &W, &newLia, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    dpd_file2_close(&newLIA);
    dpd_file2_close(&newLia);
  }

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");

    /* L1 RHS += -1/2 Lmnae*Wiemn */
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    dpd_contract442(&W, &LIjAb, &newLIA, 0, 2, -1, 1);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&W);

    dpd_file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 2, 3, "New Lia");

    /* L1 RHS += -1/2 Lmnae*Wiemn */
    dpd_buf4_init(&WMBIJ, CC2_HET1, 0, 20, 2, 20, 2, 0, "CC2 WMBIJ (MB,I>J)");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1, 1);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&WMBIJ);

    dpd_buf4_init(&WMbIj, CC2_HET1, 0, 24, 22, 24, 22, 0, "CC2 WMbIj (Mb,Ij)");
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1, 1);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&WMbIj);

    dpd_buf4_init(&Wmbij, CC2_HET1, 0, 30, 12, 30, 12, 0, "CC2 Wmbij (mb,i>j)");
    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1, 1);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&Wmbij);

    dpd_buf4_init(&WmBiJ, CC2_HET1, 0, 27, 23, 27, 23, 0, "CC2 WmBiJ (mB,iJ)");
    dpd_buf4_init(&LiJaB, CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1, 1);
    dpd_buf4_close(&LiJaB);
    dpd_buf4_close(&WmBiJ);

    dpd_file2_close(&newLIA);
    dpd_file2_close(&newLia);
  }

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");

    /* L1 RHS += Gbj*<ij|ab> */
    dpd_file2_init(&G, CC_TMP0, L_irr, 1, 0, "CC2 GAI");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    dpd_contract422(&D, &G, &newLIA, 1, 0, 1, 1);
    dpd_buf4_close(&D);
    dpd_file2_close(&G);

    dpd_file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&newLIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&newLia, CC_LAMBDA, L_irr, 2, 3, "New Lia");

    dpd_file2_init(&GAI, CC_LAMBDA, L_irr, 1, 0, "CC2 GAI");
    dpd_file2_init(&Gai, CC_LAMBDA, L_irr, 3, 2, "CC2 Gai");

    /* L1 RHS += Gbj*<ij|ab> */
    /** AA **/
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_contract422(&D, &GAI, &newLIA, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_contract422(&D, &Gai, &newLIA, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    /** BB**/
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_contract422(&D, &Gai, &newLia, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_contract422(&D, &GAI, &newLia, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_file2_close(&Gai);
    dpd_file2_close(&GAI);

    dpd_file2_close(&newLIA);
    dpd_file2_close(&newLia);
  }

  if(params.ref == 0) { /** RHF **/
    /* newLia * Dia */
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
    /*dpd_file2_print(&newLIA,outfile);*/
    dpd_file2_close(&LIA);

    dpd_file2_copy(&newLIA, CC_LAMBDA, "New Lia");  /* spin-adaptation for RHF */
    dpd_file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/

    /* newLia * Dia */
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
  else if(params.ref == 2) { /** UHF **/

    /* newLia * Dia */
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
