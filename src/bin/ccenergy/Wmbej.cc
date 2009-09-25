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

/* Wmbej_build(): Build the Wmbej intermediate.
** 
** Wmbej = <mb||ej> + t(j,f) <mb||ef> - t(n,b) <mn||ej>
**         - [1/2 t(jn,fb) + t(j,f) t(n,b)] <mn||ef>
**
** Spin cases for UHF and ROHF orbitals:
** -------------------------------------
**
**
** TDC
** May 2000
*/

void Wmbej_build(void)
{
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ, W; 
  dpdbuf4 C, D, E, F, X, tIAjb, tiaJB, t2, Y, Z;
  dpdfile2 tIA, tia;
  int Gmb, mb, Gj, Ge, Gf, nrows, ncols, nlinks;

#ifdef TIME_CCENERGY
  timer_on("C->Wmbej");
#endif

  /* W(mb,je) <-- <mb||ej> */

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMbeJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "WMbEj");
    dpd_buf4_close(&D);
  }
  else if(params.ref == 1) { /*** ROHF ***/

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMBEJ", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "Wmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WmBEj", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMbeJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "WMbEj");
    dpd_buf4_copy(&D, CC_TMP0, "WmBeJ");
    dpd_buf4_close(&D);

  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&C, CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMBEJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "Wmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_buf4_scmcopy(&D, CC_TMP0, "WMbEj", 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_buf4_scmcopy(&D, CC_TMP0, "WmBeJ", 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&C, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WmBEj", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMbeJ", -1);
    dpd_buf4_close(&C);

  }

#ifdef TIME_CCENERGY
  timer_off("C->Wmbej");
  timer_on("F->Wmbej");
#endif
  
  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&F, &tIA, &WMbEj, 3, 1, 0, 1, 1); /* should be OOC-capable in libdpd */
    dpd_buf4_close(&WMbEj);
    dpd_buf4_close(&F);

    /*
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");

    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_sort(&Z, CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_axpy(&Z, &WMbeJ, 1.0);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_close(&Z);

    dpd_buf4_close(&F);
    */

    /* W(Mb,Je) <-- t(J,F) <Mb|Fe> */
    /* OOC code added to replace above on 3/23/05, TDC */
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_file2_mat_init(&tIA);
    dpd_file2_mat_rd(&tIA);

    for(Gmb=0; Gmb < moinfo.nirreps; Gmb++) {
      dpd_buf4_mat_irrep_row_init(&W, Gmb);
      dpd_buf4_mat_irrep_row_init(&F, Gmb);

      for(mb=0; mb < F.params->rowtot[Gmb]; mb++) {
	dpd_buf4_mat_irrep_row_rd(&W, Gmb, mb);
	dpd_buf4_mat_irrep_row_rd(&F, Gmb, mb);

	for(Gj=0; Gj < moinfo.nirreps; Gj++) {
	  Gf = Gj;  /* T1 is totally symmetric */
	  Ge = Gmb ^ Gf; /* <mb|fe> is totally symmetric */

	  nrows = moinfo.occpi[Gj];
	  ncols = moinfo.virtpi[Ge];
	  nlinks = moinfo.virtpi[Gf];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('n','n',nrows,ncols,nlinks,-1.0,tIA.matrix[Gj][0],nlinks,
		    &F.matrix[Gmb][0][F.col_offset[Gmb][Gf]],ncols,1.0,
		    &W.matrix[Gmb][0][W.col_offset[Gmb][Gj]],ncols);
	}

	dpd_buf4_mat_irrep_row_wrt(&W, Gmb, mb);
      }

      dpd_buf4_mat_irrep_row_close(&F, Gmb);
      dpd_buf4_mat_irrep_row_close(&W, Gmb);
    }

    dpd_file2_mat_close(&tIA);
    dpd_buf4_close(&W);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_contract424(&F, &tIA, &WMBEJ, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WMBEJ);
    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    dpd_contract424(&F, &tia, &Wmbej, 3, 1, 0, 1, 1);
    dpd_buf4_close(&Wmbej);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&F, &tia, &WMbEj, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_contract424(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WmBeJ);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract244(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_contract244(&tia, &F, &WmBEj, 1, 2, 1, -1, 1);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract244(&tIA, &F, &W, 1, 2, 1, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract244(&tia, &F, &W, 1, 2, 1, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }

#ifdef TIME_CCENERGY
  timer_off("F->Wmbej");
  timer_on("E->Wmbej");
#endif

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&E, &tIA, &WMbEj, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract424(&E, &tIA, &WMbeJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA);  
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_contract424(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMBEJ);
    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    dpd_contract424(&E, &tia, &Wmbej, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Wmbej);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&E, &tia, &WMbEj, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_contract424(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WmBeJ);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract424(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_contract424(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    dpd_contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_init(&E, CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    dpd_contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract424(&E, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    dpd_contract424(&E, &tIA, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);
  }

#ifdef TIME_CCENERGY
  timer_off("E->Wmbej");
#endif

  /* Convert to (ME,JB) for remaining terms */

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_buf4_sort(&WMbEj, CC_HBAR, prsq, 10, 10, "WMbEj");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_sort(&WMbeJ, CC_HBAR, psrq, 10, 10, "WMbeJ");
    dpd_buf4_close(&WMbeJ);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_buf4_sort(&WMBEJ, CC_HBAR, prsq, 10, 10, "WMBEJ");
    dpd_buf4_close(&WMBEJ);

    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    dpd_buf4_sort(&Wmbej, CC_HBAR, prsq, 10, 10, "Wmbej");
    dpd_buf4_close(&Wmbej);

    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_buf4_sort(&WMbEj, CC_HBAR, prsq, 10, 10, "WMbEj");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_buf4_sort(&WmBeJ, CC_HBAR, prsq, 10, 10, "WmBeJ");
    dpd_buf4_close(&WmBeJ);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_sort(&WMbeJ, CC_HBAR, psrq, 10, 10, "WMbeJ");
    dpd_buf4_close(&WMbeJ);

    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_buf4_sort(&WmBEj, CC_HBAR, psrq, 10, 10, "WmBEj");
    dpd_buf4_close(&WmBEj);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_sort(&W, CC_HBAR, prsq, 20, 20, "WMBEJ");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_sort(&W, CC_HBAR, prsq, 30, 30, "Wmbej");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_sort(&W, CC_HBAR, prsq, 20, 30, "WMbEj");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_sort(&W, CC_HBAR, prsq, 30, 20, "WmBeJ");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_sort(&W, CC_HBAR, psrq, 24, 24, "WMbeJ");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_sort(&W, CC_HBAR, psrq, 27, 27, "WmBEj");
    dpd_buf4_close(&W);
  }

#ifdef TIME_CCENERGY
  timer_on("X->Wmbej");
#endif

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    /*** ABAB ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    dpd_contract444(&D, &t2, &W, 0, 1, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_contract444(&D, &t2, &W, 0, 1, -0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABBA ***/
  
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_contract444(&D, &t2, &W, 0, 1, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /*** AAAA ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BBBB ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABAB ***/
  
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1); 
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1); 
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BABA ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABBA ***/
  
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BAAB ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /*** AAAA ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 20, 0, 20, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 21, 20, 21, 0, "D <IJ||AB> (IA,BJ)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BBBB ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 30, 10, 30, 10, 0, "Y (me,jn)");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 31, 30, 31, 0, "D <ij||ab> (ia,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABAB ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1); 
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 20, 10, 20, 10, 0, "Y (ME,jn)");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 31, 20, 31, 0, "D <Ij|Ab> (IA,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BABA ***/

    dpd_buf4_init(&W, CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1); 
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 30, 0, 30, 0, 0, "Y (me,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 21, 30, 21, 0, "D <Ij|Ab> (ia,BJ)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABBA ***/
  
    dpd_buf4_init(&W, CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 24, 22, 24, 22, 0, "Y (Me,Jn)");
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BAAB ***/
  
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 27, 23, 27, 23, 0, "Y (mE,jN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  }

#ifdef TIME_CCENERGY
  timer_off("X->Wmbej");
#endif
}
}} // namespace psi::ccenergy
