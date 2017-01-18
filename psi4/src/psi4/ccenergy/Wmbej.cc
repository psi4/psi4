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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"

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

void CCEnergyWavefunction::Wmbej_build(void)
{
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ, W;
  dpdbuf4 C, D, E, F, X, tIAjb, tiaJB, t2, Y, Z;
  dpdfile2 tIA, tia;
  int Gmb, mb, Gj, Ge, Gf, nrows, ncols, nlinks;

#ifdef TIME_CCENERGY
  timer_on("C->Wmbej");
#endif

  /* W(mb,je) <-- <mb||ej> */

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WMbEj");
    global_dpd_->buf4_close(&D);
  }
  else if(params_.ref == 1) { /*** ROHF ***/

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMBEJ", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "Wmbej", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WmBEj", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WMbEj");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WmBeJ");
    global_dpd_->buf4_close(&D);

  }
  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMBEJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "Wmbej", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "WMbEj", 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "WmBeJ", 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WmBEj", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);

  }

#ifdef TIME_CCENERGY
  timer_off("C->Wmbej");
  timer_on("F->Wmbej");
#endif

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&F, &tIA, &WMbEj, 3, 1, 0, 1, 1); /* should be OOC-capable in libdpd */
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_close(&F);

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
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->file2_mat_init(&tIA);
    global_dpd_->file2_mat_rd(&tIA);

    for(Gmb=0; Gmb < moinfo_.nirreps; Gmb++) {
      global_dpd_->buf4_mat_irrep_row_init(&W, Gmb);
      global_dpd_->buf4_mat_irrep_row_init(&F, Gmb);

      for(mb=0; mb < F.params->rowtot[Gmb]; mb++) {
	global_dpd_->buf4_mat_irrep_row_rd(&W, Gmb, mb);
	global_dpd_->buf4_mat_irrep_row_rd(&F, Gmb, mb);

    for(Gj=0; Gj < moinfo_.nirreps; Gj++) {
	  Gf = Gj;  /* T1 is totally symmetric */
	  Ge = Gmb ^ Gf; /* <mb|fe> is totally symmetric */

      nrows = moinfo_.occpi[Gj];
      ncols = moinfo_.virtpi[Ge];
      nlinks = moinfo_.virtpi[Gf];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('n','n',nrows,ncols,nlinks,-1.0,tIA.matrix[Gj][0],nlinks,
		    &F.matrix[Gmb][0][F.col_offset[Gmb][Gf]],ncols,1.0,
		    &W.matrix[Gmb][0][W.col_offset[Gmb][Gj]],ncols);
	}

	global_dpd_->buf4_mat_irrep_row_wrt(&W, Gmb, mb);
      }

      global_dpd_->buf4_mat_irrep_row_close(&F, Gmb);
      global_dpd_->buf4_mat_irrep_row_close(&W, Gmb);
    }

    global_dpd_->file2_mat_close(&tIA);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&F);

    global_dpd_->file2_close(&tIA);
  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->contract424(&F, &tIA, &WMBEJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->contract424(&F, &tia, &Wmbej, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&F, &tia, &WMbEj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->contract424(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract244(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract244(&tia, &F, &WmBEj, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&F);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params_.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract244(&tIA, &F, &W, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract244(&tia, &F, &W, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }

#ifdef TIME_CCENERGY
  timer_off("F->Wmbej");
  timer_on("E->Wmbej");
#endif

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&E, &tIA, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&E, &tIA, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->file2_close(&tIA);
  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->contract424(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->contract424(&E, &tia, &Wmbej, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&E, &tia, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->contract424(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WmBeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract424(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params_.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }

#ifdef TIME_CCENERGY
  timer_off("E->Wmbej");
#endif

  /* Convert to (ME,JB) for remaining terms */

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC_HBAR, prsq, 10, 10, "WMbEj");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC_HBAR, psrq, 10, 10, "WMbeJ");
    global_dpd_->buf4_close(&WMbeJ);

  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->buf4_sort(&WMBEJ, PSIF_CC_HBAR, prsq, 10, 10, "WMBEJ");
    global_dpd_->buf4_close(&WMBEJ);

    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->buf4_sort(&Wmbej, PSIF_CC_HBAR, prsq, 10, 10, "Wmbej");
    global_dpd_->buf4_close(&Wmbej);

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC_HBAR, prsq, 10, 10, "WMbEj");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->buf4_sort(&WmBeJ, PSIF_CC_HBAR, prsq, 10, 10, "WmBeJ");
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC_HBAR, psrq, 10, 10, "WMbeJ");
    global_dpd_->buf4_close(&WMbeJ);

    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->buf4_sort(&WmBEj, PSIF_CC_HBAR, psrq, 10, 10, "WmBEj");
    global_dpd_->buf4_close(&WmBEj);
  }
  else if(params_.ref == 2) { /** UHF **/

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prsq, 20, 20, "WMBEJ");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prsq, 30, 30, "Wmbej");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prsq, 20, 30, "WMbEj");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prsq, 30, 20, "WmBeJ");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psrq, 24, 24, "WMbeJ");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psrq, 27, 27, "WmBEj");
    global_dpd_->buf4_close(&W);
  }

#ifdef TIME_CCENERGY
  timer_on("X->Wmbej");
#endif

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /*** ABAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract444(&D, &t2, &W, 0, 1, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->contract444(&D, &t2, &W, 0, 1, -0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABBA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 1, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);
  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /*** AAAA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BBBB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BABA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABBA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BAAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params_.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /*** AAAA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 21, 20, 21, 0, "D <IJ||AB> (IA,BJ)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BBBB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Y (me,jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 31, 30, 31, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 20, 10, 20, 10, 0, "Y (ME,jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 31, 20, 31, 0, "D <Ij|Ab> (IA,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BABA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 30, 0, 30, 0, 0, "Y (me,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 21, 30, 21, 0, "D <Ij|Ab> (ia,BJ)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABBA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Y (Me,Jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BAAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Y (mE,jN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  }

#ifdef TIME_CCENERGY
  timer_off("X->Wmbej");
#endif
}
}} // namespace psi::ccenergy
