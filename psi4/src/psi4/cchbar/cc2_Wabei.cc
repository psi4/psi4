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
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <string>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* cc2_Wabei(): Compute the Wabei matrix from CC2 theory, which is
** given in spin orbitals as:
**
** Wabei = <ab||ei> - P(ab) t_m^a <mb||ei> + t_i^f <ab||ef>
**         - P(ab) t_i^f t_m^b <am||ef> + t_m^a t_n^b <mn||ei>
**         + t_m^a t_i^f t_n^b <mn||ef>
**
** The basic strategy for this code is to generate two intermediate
** quantities, Z1(Ab,EI) and Z2(Ei,Ab), which are summed in the final
** step to give the complete W(Ei,Ab) intermediate.  This is sorted
** to W(iE,bA) storage for use in the triples equations.
**
** TDC, Feb 2004
*/

void purge_cc2_Wabei(void);

void cc2_Wabei_build(void)
{
  int omit = 0;
  int e, E;
  int Gef, Gab, Gei, Ge, Gf, Gi;
  int nrows, ncols, nlinks;
  dpdfile2 t1, tIA, tia;
  dpdbuf4 Z, Z1, Z2, Z3;
  dpdbuf4 B, C, D, F, W;

  timer_on("F->Wabei");
  if(params.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC_TMP0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->buf4_close(&F);

  }

  else if (params.ref == 1) { /* ROHF */
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /* W(A>B,EI) <--- <AB||EI> */
    /* W(a>b,ei) <--- <ab||ei> */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 7, 11, 5, 1, "F <ai|bc>");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP2, rspq, 7, 11, "CC2 WABEI (A>B,EI)");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP2, rspq, 7, 11, "CC2 Wabei (a>b,ei)");
    global_dpd_->buf4_close(&F);

    /* W(Ab,Ei) <--- <Ab|Ei> */
    /* W(aB,eI) <--- <aB|eI> */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP2, rspq, 5, 11, "CC2 WAbEi (Ab,Ei)");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP2, rspq, 5, 11, "CC2 WaBeI (aB,eI)");
    global_dpd_->buf4_close(&F);
  }

  else if (params.ref == 2) { /* UHF */
    /* W(A>B,EI) <--- <AB||EI> */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, rspq, 7, 21, "CC2 WABEI (A>B,EI)");
    global_dpd_->buf4_close(&F);

    /* W(a>b,ei) <--- <ab||ei> */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 17, 31,15, 1, "F <ai|bc>");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, rspq, 17, 31, "CC2 Wabei (a>b,ei)");
    global_dpd_->buf4_close(&F);

    /* W(Ab,Ei) <--- <Ab|Ei> */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    global_dpd_->buf4_copy(&F, PSIF_CC_TMP0, "CC2 WAbEi (Ab,Ei)");
    global_dpd_->buf4_close(&F);

    /* W(aB,eI) <--- <aB|eI> */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, psrq, 29, 25, "CC2 WaBeI (aB,eI)");
    global_dpd_->buf4_close(&F);
  }
  timer_off("F->Wabei");

  timer_on("B->Wabei");
  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&t1);
    global_dpd_->file2_mat_rd(&t1);

    /* WEbEi <-- <Ab|Ef> * t(i,f) */

    /* Plus Combination */
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 8, 11, 8, 0, "Z1(ei,a>=b)");

    for(Gef=0; Gef < moinfo.nirreps; Gef++) {
      Gei = Gab = Gef; /* W and B are totally symmetric */
      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
	Gf = Ge ^ Gef; Gi = Gf;  /* t1 is totally symmetric */
	B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
	Z1.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi],Z1.params->coltot[Gei]);
	nrows = moinfo.occpi[Gi];
	ncols = Z1.params->coltot[Gef];
	nlinks = moinfo.virtpi[Gf];
	if(nrows && ncols && nlinks) {
	  for(E=0; E < moinfo.virtpi[Ge]; E++) {
	    e = moinfo.vir_off[Ge] + E;
	    global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
	    C_DGEMM('n','n',nrows,ncols,nlinks,0.5,t1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
		    0.0,Z1.matrix[Gei][0],ncols);
	    global_dpd_->buf4_mat_irrep_wrt_block(&Z1, Gei, Z1.row_offset[Gei][e], moinfo.occpi[Gi]);
	  }
	}
	global_dpd_->free_dpd_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
	global_dpd_->free_dpd_block(Z1.matrix[Gef], moinfo.occpi[Gi], Z1.params->coltot[Gei]);
      }
    }

    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&B);

    /* Minus Combination */
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 9, 11, 9, 0, "Z2(ei,a>=b)");

    for(Gef=0; Gef < moinfo.nirreps; Gef++) {
      Gei = Gab = Gef; /* W and B are totally symmetric */
      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
	Gf = Ge ^ Gef; Gi = Gf;  /* t1 is totally symmetric */
	B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
	Z2.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi],Z2.params->coltot[Gei]);
	nrows = moinfo.occpi[Gi];
	ncols = Z2.params->coltot[Gef];
	nlinks = moinfo.virtpi[Gf];
	if(nrows && ncols && nlinks) {
	  for(E=0; E < moinfo.virtpi[Ge]; E++) {
	    e = moinfo.vir_off[Ge] + E;
	    global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
	    C_DGEMM('n','n',nrows,ncols,nlinks,0.5,t1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
		    0.0,Z2.matrix[Gei][0],ncols);
	    global_dpd_->buf4_mat_irrep_wrt_block(&Z2, Gei, Z2.row_offset[Gei][e], moinfo.occpi[Gi]);
	  }
	}
	global_dpd_->free_dpd_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
	global_dpd_->free_dpd_block(Z2.matrix[Gef], moinfo.occpi[Gi], Z2.params->coltot[Gei]);
      }
    }

    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&B);

    global_dpd_->file2_mat_close(&t1);
    global_dpd_->file2_close(&t1);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 8, 0, "Z1(ei,a>=b)");
    global_dpd_->buf4_axpy(&Z1, &W, 1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 9, 0, "Z2(ei,a>=b)");
    global_dpd_->buf4_axpy(&Z2, &W, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);
  }
  else if (params.ref == 1) { /* ROHF */

    /** W(A>B,EI) <--- <AB||EF> * t1[I][F] **/
    /** W(a>b,ei) <--- <ab||ef> * t1[i][f] **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, 0, 7, 11, 7, 11, 0, "CC2 WABEI (A>B,EI)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <ab|cd>");
    global_dpd_->contract424(&B, &tIA, &W, 3, 1, 0, 0.5, 1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, 0, 7, 11, 7, 11, 0, "CC2 Wabei (a>b,ei)");
    global_dpd_->contract424(&B, &tia, &W, 3, 1, 0, 0.5, 1);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&W);

    /** W(Ab,Ei) <--- <Ab|Ef> * t1[i][f] **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "CC2 WAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract424(&B, &tia, &W, 3, 1, 0, 0.5, 1);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&W);

    /** W(aB,eI) <--- t1[I][F] * <aB|eF>  **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "CC2 WaBeI (aB,eI)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract424(&B, &tIA, &W, 3, 1, 0, 0.5, 1.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }

  else if (params.ref == 2) { /* UHF */
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /** W(A>B,EI) <--- <AB||EF> * t1[I][F] **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "CC2 WABEI (A>B,EI)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
    global_dpd_->contract424(&B, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&W);

    /** W(a>b,ei) <--- <ab||ef> * t1[i][f] **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 17, 31, 17, 31, 0, "CC2 Wabei (a>b,ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 17, 15, 15, 15, 1, "B <ab|cd>");
    global_dpd_->contract424(&B, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&W);

    /** W(Ab,Ei) <--- <Ab|Ef> * t1[i][f] **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "CC2 WAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->contract424(&B, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&W);

    /** W(aB,eI) <--- t1[I][F] * <aB|eF>  **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "CC2 ZBaIe (Ba,Ie)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->contract244(&tIA, &B, &Z, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, qpsr, 29, 25, "CC2 WaBeI (aB,eI)", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  timer_off("B->Wabei");

  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "CC2 ZMbEj (Ej,Mb)");
    global_dpd_->contract244(&t1, &Z, &W, 0, 2, 1, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 11, 11, 11, 0, "CC2 ZMbeJ (eJ,bM)");
    global_dpd_->contract424(&Z, &t1, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&t1);
  }
  else if (params.ref == 2) { /* UHF */
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /** AAA ***/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "CC2 ZABEI (AB,EI)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 ZMBEJ");
    global_dpd_->contract244(&tIA, &Z, &Z1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, qprs, 5, 21, "CC2 ZABEI (BA,EI)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "CC2 ZABEI (BA,EI)");
    global_dpd_->buf4_axpy(&Z, &Z1, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 21, 7, 21, 0, "CC2 WABEI (A>B,EI)");
    global_dpd_->buf4_axpy(&Z1, &W, 1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);

    /** ABAB **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "CC2 ZAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 ZMbEj");
    global_dpd_->contract244(&tIA, &Z, &Z1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 29, 26, 29, 26, 0, "CC2 ZAbEi (bA,Ei)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 26, 27, 26, 0, "CC2 ZmBEj");
    global_dpd_->contract244(&tia, &Z, &Z1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, qprs, 28, 26, "CC2 ZAbEi (Ab,Ei)", 1);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "CC2 ZAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "CC2 WAbEi (Ab,Ei)");
    global_dpd_->buf4_axpy(&Z, &W, -1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    /** BABA **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "CC2 ZaBeI (aB,eI)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 ZmBeJ");
    global_dpd_->contract244(&tia, &Z, &Z1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 28, 25, 28, 25, 0, "CC2 ZaBeI (Ba,eI)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 25, 24, 25, 0, "CC2 ZMbeJ");
    global_dpd_->contract244(&tIA, &Z, &Z1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, qprs, 29, 25, "CC2 ZaBeI (aB,eI)", 1);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "CC2 ZaBeI (aB,eI)");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "CC2 WaBeI (aB,eI)");
    global_dpd_->buf4_axpy(&Z, &W, -1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    /** BBBB **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "CC2 Zabei (ab,ei)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Zmbej");
    global_dpd_->contract244(&tia, &Z, &Z1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, qprs, 15, 31, "CC2 Zabei (ba,ei)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "CC2 Zabei (ba,ei)");
    global_dpd_->buf4_axpy(&Z, &Z1, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 15, 31, 17, 31, 0, "CC2 Wabei (a>b,ei)");
    global_dpd_->buf4_axpy(&Z1, &W, 1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }

  timer_on("Wabei_sort");
  if (params.ref == 0) { /* RHF */

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
    global_dpd_->buf4_scm(&W, 0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC2_HET1, rspq, 5, 11, "CC2 WAbEi", 1);
    global_dpd_->buf4_close(&W);

    if(params.wfn == "EOM_CC2") {
      global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
      global_dpd_->buf4_copy(&W, PSIF_CC2_HET1, "CC2 WAbEi (Ei,Ab)");
      global_dpd_->buf4_close(&W);
    }

  }
  if (params.ref == 1) { /* ROHF */

    /* sort to Wabei (ei,ab) */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, 0, 7, 11, 7, 11, 0, "CC2 WABEI (A>B,EI)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 11, 7, "CC2 WABEI (EI,A>B)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, 0, 7, 11, 7, 11, 0, "CC2 Wabei (a>b,ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 11, 7, "CC2 Wabei (ei,a>b)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "CC2 WAbEi (Ab,Ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 11, 5, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "CC2 WaBeI (aB,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 11, 5, "CC2 WaBeI (eI,aB)");
    global_dpd_->buf4_close(&W);

    /* purge before final sort */

    /*     purge_cc2_Wabei(); */
  }
  else if (params.ref == 2) { /* UHF */

    /* sort to Wabei (ei,ab) */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "CC2 WABEI (A>B,EI)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 21, 7, "CC2 WABEI (EI,A>B)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 17, 31, 17, 31, 0, "CC2 Wabei (a>b,ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 31, 17, "CC2 Wabei (ei,a>b)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "CC2 WAbEi (Ab,Ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 26, 28, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "CC2 WaBeI (aB,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 25, 29, "CC2 WaBeI (eI,aB)");
    global_dpd_->buf4_close(&W);
  }
  timer_off("Wabei_sort");

}

void purge_cc2_Wabei(void) {
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  /* Purge Wabei matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC_TMP2, 0, 11, 7,"CC2 WABEI (EI,A>B)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      e = W.params->roworb[h][ei][0];
      esym = W.params->psym[e];
      E = e - vir_off[esym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
        a = W.params->colorb[h][ab][0];
        b = W.params->colorb[h][ab][1];
        asym = W.params->rsym[a];
        bsym = W.params->ssym[b];
        A = a - vir_off[asym];
        B = b - vir_off[bsym];
        if ((E >= (virtpi[esym] - openpi[esym])) ||
            (A >= (virtpi[asym] - openpi[asym])) ||
            (B >= (virtpi[bsym] - openpi[bsym])) )
          W.matrix[h][ei][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC_TMP2, 0, 11, 7,"CC2 Wabei (ei,a>b)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      i = W.params->roworb[h][ei][1];
      isym = W.params->qsym[i];
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
        if (I >= (occpi[isym] - openpi[isym]))
          W.matrix[h][ei][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC_TMP2, 0, 11, 5,"CC2 WAbEi (Ei,Ab)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      e = W.params->roworb[h][ei][0];
      i = W.params->roworb[h][ei][1];
      esym = W.params->psym[e];
      isym = W.params->qsym[i];
      E = e - vir_off[esym];
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
        a = W.params->colorb[h][ab][0];
        asym = W.params->rsym[a];
        bsym = W.params->ssym[b];
        A = a - vir_off[asym];
        if ((E >= (virtpi[esym] - openpi[esym])) ||
            (I >= (occpi[isym] - openpi[isym])) ||
            (A >= (virtpi[asym] - openpi[asym])) )
          W.matrix[h][ei][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC_TMP2, 0, 11, 5,"CC2 WaBeI (eI,aB)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      for(ab=0; ab<W.params->coltot[h]; ab++) {
        b = W.params->colorb[h][ab][1];
        bsym = W.params->ssym[b];
        B = b - vir_off[bsym];
        if (B >= (virtpi[bsym] - openpi[bsym]))
          W.matrix[h][ei][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

}


}} // namespace psi::cchbar
