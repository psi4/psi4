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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* compute [ [H,eT1] , C1 ]
 *
 * [H,eT^1] matrix elements are in CC3_HET1
 *
 * C1 has symmetry C_irr and is read from EOM_CME with name "CME i"
 *
 * illustrative code that computes [H,C1] is in cc3_HC1.c
 *
 * output matrix elements are written to CC_HC1ET1
 *
*/
void HC1ET1_Wmbij(int i, int C_irr);
void HC1ET1_Wabei(int i, int C_irr);
void HC1ET1_Wmbij_rhf(int i, int C_irr);
void HC1ET1_Wabei_rhf(int i, int C_irr);

void cc3_HC1ET1 (int i, int C_irr) {

  /* only need new Wmbij and Wabei for CC3 EOM energies */

  if (params.ref == 0) {
    HC1ET1_Wmbij_rhf(i, C_irr);
    HC1ET1_Wabei_rhf(i, C_irr);
  }
  else {
    HC1ET1_Wmbij(i, C_irr);
    HC1ET1_Wabei(i, C_irr);
  }

  return;
}

void HC1ET1_Wmbij(int i, int C_irr)
{
  double dot;
  dpdbuf4 C, D, E, F, Ht, W, W1, X, Z;
  dpdfile2 CME, Cme;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if (params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }
  else if (params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /**** Term I ****/

    /***** Ht (MB,I>J) <--- -WMNIJ * CNB *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 20, 2, 20, 2, 0, "Ht_WMBIJ (MB,I>J)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 2, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    global_dpd_->contract424(&W, &CME, &Ht, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Ht);

    /***** Ht (mb,i>j) <--- -Wmnij * Cnb *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 30, 12, 30, 12, 0, "Ht_Wmbij (mb,i>j)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 12, 12, 12, 0, "CC3 Wmnij (m>n,i>j)");
    global_dpd_->contract424(&W, &Cme, &Ht, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Ht);

    /***** Ht (Mb,Ij) <--- -WMnIj * Cnb *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 24, 22, 24, 22, 0, "Ht_WMbIj (Mb,Ij)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->contract424(&W, &Cme, &Ht, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Ht);

    /***** Ht (mB,iJ) <--- ZBmJi <--- CNB * WNmJi *****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 26, 22, 26, 22, 0, "Z (Bm,Ji)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->contract244(&CME, &W, &Z, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_sort(&Z, PSIF_CC3_HC1ET1, qpsr, 27, 23, "Ht_WmBiJ (mB,iJ)");
    global_dpd_->buf4_close(&Z);

    /**** Term II ****/

    /***** Ht (MB,I>J) <--- -P(I/J) X (MB,JI) <--- WMBEJ(ME,JB) * CIE *****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 0, 20, 0, 20, 0, "Z (MI,JB)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 20, 20, 20, 0, "CC3 WMBEJ (ME,JB)");
    global_dpd_->contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 20, 0, "X (MB,JI)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, C_irr, 20, 0, 20, 0, 0, "X (MB,JI)");
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC3_HC1ET1, pqsr, 20, 2, "Ht_WMBIJ (MB,I>J)", 1);

    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 20, 0, 20, 2, 0, "Ht_WMBIJ (MB,I>J)");
    global_dpd_->buf4_axpy(&X, &Ht, -1);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_close(&Ht);

    /***** Ht (mb,i>j) <--- -P(i/j) X (mb,ij) <--- Wmbej (me,jb) * Cie *****/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 10, 30, 10, 30, 0, "Z (mi,jb)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 30, 30, 30, 0, "CC3 Wmbej (me,jb)");
    global_dpd_->contract424(&W, &Cme, &Z, 1, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 30, 10, "X (mb,ji)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, C_irr, 30, 10, 30, 10, 0, "X (mb,ji)");
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC3_HC1ET1, pqsr, 30, 12, "Ht_Wmbij (mb,i>j)", 1);

    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 30, 10, 30, 12, 0, "Ht_Wmbij (mb,i>j)");
    global_dpd_->buf4_axpy(&X, &Ht, -1);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_close(&Ht);

    /***** Ht (Mb,Ij) <--- CIE * WMbEj *****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 0, 30, 0, 30, 0, "Z (MI,jb)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 30, 20, 30, 0, "CC3 WMbEj (ME,jb)");
    global_dpd_->contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HC1ET1, psqr, 24, 22, "Ht_WMbIj (Mb,Ij)", 1);
    global_dpd_->buf4_close(&Z);

    /***** Ht (Mb,Ij) <--- Cje * WMbeI *****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 22, 24, 22, 24, 0, "Z (Mj,Ib)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 24, 24, 24, 0, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->contract424(&W, &Cme, &Z, 1, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HC1ET1, psrq, 24, 22, "Ht_WMbIj (Mb,Ij)", -1);
    global_dpd_->buf4_close(&Z);

    /***** Ht (mB,iJ) <--- Cie * WmBiJ *****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 10, 20, 10, 20, 0, "Z (mi,JB)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 20, 30, 20, 0, "CC3 WmBeJ (me,JB)");
    global_dpd_->contract424(&W, &Cme, &Z, 1, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HC1ET1, psqr, 27, 23, "Ht_WmBiJ (mB,iJ)", 1);
    global_dpd_->buf4_close(&Z);

    /***** Ht (mB,iJ) <--- CJE * WmEiB *****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 23, 27, 23, 27, 0, "Z (mJ,iB)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 27, 27, 27, 0, "CC3 WmBEj (mE,jB)");
    global_dpd_->contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HC1ET1, psrq, 27, 23, "Ht_WmBiJ (mB,iJ)", -1);
    global_dpd_->buf4_close(&Z);


    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 20, 2, 20, 2, 0, "Ht_WMBIJ (MB,I>J)");
    global_dpd_->buf4_sort(&Ht, PSIF_CC3_HC1ET1, rspq, 2, 20, "Ht_WMBIJ (I>J,MB)");
    global_dpd_->buf4_close(&Ht);

    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 30, 12, 30, 12, 0, "Ht_Wmbij (mb,i>j)");
    global_dpd_->buf4_sort(&Ht, PSIF_CC3_HC1ET1, rspq, 12, 30, "Ht_Wmbij (i>j,mb)");
    global_dpd_->buf4_close(&Ht);

    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 24, 22, 24, 22, 0, "Ht_WMbIj (Mb,Ij)");
    global_dpd_->buf4_sort(&Ht, PSIF_CC3_HC1ET1, rspq, 22, 24, "Ht_WMbIj (Ij,Mb)");
    global_dpd_->buf4_close(&Ht);

    global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 27, 23, 27, 23, 0, "Ht_WmBiJ (mB,iJ)");
    global_dpd_->buf4_sort(&Ht, PSIF_CC3_HC1ET1, rspq, 23, 27, "Ht_WmBiJ (iJ,mB)");
    global_dpd_->buf4_close(&Ht);


    /************ TEST *************/

    /*
    dpd_buf4_init(&W, CC3_HC1ET1, 0, 0, 20, 2, 20, 0, "Ht_WMBIJ (I>J,MB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMBIJ (I>J,MB)|WMBIJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 10, 30, 12, 30, 0, "Ht_Wmbij (i>j,mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmbij (i>j,mb)|Wmbij> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 22, 24, 22, 24, 0, "Ht_WMbIj (Ij,Mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbIj (Ij,Mb)|WMbIj> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 23, 27, 23, 27, 0, "Ht_WmBiJ (iJ,mB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBiJ (iJ,mB)|WmBiJ> = %15.10lf\n", dot);
    */

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }

  return;
}

void HC1ET1_Wabei(int i, int C_irr)
{
  double dot;
  dpdfile2 CME, Cme, tIA, tia;
  dpdbuf4 Ht, Z, Z1, Z2, Z3, B, C, D, E, F, W, X;

  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if (params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }
  else if (params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");


    /**** Term I ****/

    /***** Ht_WABEI <--- -P(A/B) CMA * WMBEI *****/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 5, 21, 5, 21, 0, "Z (AB,EI)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 21, 20, 21, 0, "CC3 WMBEJ (MB,EJ)");
    global_dpd_->contract244(&CME, &W, &Z, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qprs, 7, 21, "Ht_WABEI (A>B,EI)");
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 5, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    global_dpd_->buf4_axpy(&Z, &Ht, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Ht);

    /***** Ht_Wabei <--- Xbaei <--- -P(a/b) Zabei <-- Cma * Wmbei *****/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 15, 31, 15, 31, 0, "Z (ab,ei)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 31, 30, 31, 0, "CC3 Wmbej (mb,ej)");
    global_dpd_->contract244(&Cme, &W, &Z, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qprs, 17, 31, "Ht_Wabei (a>b,ei)");
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 15, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    global_dpd_->buf4_axpy(&Z, &Ht, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Ht);

    /***** Ht_WAbEi <--- -CMA * WMbEi *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 26, 24, 26, 0, "CC3 WMbEj (Mb,Ej)");
    global_dpd_->contract244(&CME, &W, &Ht, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Ht);

    /***** Ht_WAbEi <--- WAmEi * Cmb *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 26, 28, 26, 28, 0, "Ht_WAbEi (Ei,Ab)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 26, 26, 26, 26, 0, "CC3 WmBEj (Bm,Ej)");
    global_dpd_->contract424(&W, &Cme, &Ht, 1, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Ht);

    /***** Ht_WaBeI <--- -Cma * WmBeI *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 29, 25, 29, 25, 0, "Ht_WaBeI (aB,eI)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 25, 27, 25, 0, "CC3 WmBeJ (mB,eJ)");
    global_dpd_->contract244(&Cme, &W, &Ht, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Ht);

    /***** Ht_WaBeI <--- WaMeI * CMB *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 25, 29, 25, 29, 0, "Ht_WaBeI (eI,aB)");
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 25, 25, 25, 25, 0, "CC3 WMbeJ (bM,eJ)");
    global_dpd_->contract424(&W, &CME, &Ht, 1, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Ht);

    /**** Term II ****/

    /***** Ht_WABEI <--- <AB||EF> * CIF *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 7, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
    global_dpd_->contract424(&B, &CME, &Ht, 3, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Ht);


    /***** Ht_Wabei <--- <ab||ef> * Cif *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 17, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 17, 15, 15, 15, 1, "B <ab|cd>");
    global_dpd_->contract424(&B, &Cme, &Ht, 3, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Ht);

    /***** HAbEi <--- <Ab|Ef> * Cif *****/
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->contract424(&B, &Cme, &Ht, 3, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Ht);

    /***** HaBeI <--- C[I][F] * <aB|eF>  *****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 24, 28, 24, 28, 0, "Z (Ie,Ba)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->contract244(&CME, &B, &Z, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, qpsr, 25, 29, "Ht_WaBeI (eI,aB)", 1);
    global_dpd_->buf4_close(&Z);


    /** term III **/

    /** H( A>B,EI) <--- -P(A/B)( <AM||EF>*C[I][F]*t1[M][B] ) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 20, 21, 20, 21, 0, "Z (MA,EI)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 5, 21, 5, 21, 0, "Z1(BA,EI)");
    global_dpd_->contract244(&tIA, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, qprs, 7, 21, "Ht_WABEI (A>B,EI)", 1);
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 5, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    global_dpd_->buf4_axpy(&Z1, &Ht, -1);
    global_dpd_->buf4_close(&Ht);
    global_dpd_->buf4_close(&Z1);


    /** H(a>b,ei) <--- -P(a/b)( <am||ef>*C[i][f]*t1[m][b] ) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 30, 31, 30, 31, 0, "Z (ma,ei)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &Cme, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 15, 31, 15, 31, 0, "Z1(ba,ei)");
    global_dpd_->contract244(&tia, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, qprs, 17, 31, "Ht_Wabei (a>b,ei)", 1);
    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 15, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    global_dpd_->buf4_axpy(&Z1, &Ht, -1);
    global_dpd_->buf4_close(&Ht);
    global_dpd_->buf4_close(&Z1);

    /** H(Ab,Ei) <--- -<Am|Ef>*C[i][f]*t1[m][b] **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 26, 26, 26, 26, 0, "Z(Am,Ei)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    global_dpd_->contract424(&F, &Cme, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 26, 28, 26, 28, 0, "Ht_WAbEi (Ei,Ab)");
    global_dpd_->contract424(&Z, &tia, &Ht, 1, 0, 0, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Ht);

    /** Ht(Ab,Ei) <--- -<Mb|Ef>*C[i][f]*t1[M][A] **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract424(&F, &Cme, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    global_dpd_->contract244(&tIA, &Z, &Ht, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Ht);

    /** Ht(aB,eI) <--- -<aM|eF>*C[I][F]*t1[M][B] **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 25, 25, 25, 25, 0, "Z(aM,eI)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    global_dpd_->contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 25, 29, 25, 29, 0, "Ht_WaBeI (eI,aB)");
    global_dpd_->contract424(&Z, &tIA, &Ht, 1, 0, 0, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Ht);

    /** Ht(aB,eI) <--- -C[I][F] * <mB|eF> * t1[m][a] **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 27, 25, 27, 25, 0, "Z(mB,eI)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 29, 25, 29, 25, 0, "Ht_WaBeI (aB,eI)");
    global_dpd_->contract244(&tia, &Z, &Ht, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Ht);


    /** term 6 **/

    /** Ht(A>B,EI) <---  0.5 * P(A/B) <MN||EF> * t1[M][A]*C[I][F]*t1[N][B] **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 2, 21, 2, 21, 0, "Z (M>N,EI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->contract424(&D, &CME, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 0, 21, 2, 21, 0, "Z (M>N,EI)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 21, 20, 21, 20, 0, "Z1 (EI,MB)");
    global_dpd_->contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 5, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    global_dpd_->contract244(&tIA, &Z1, &Ht, 0, 2, 0, 1, 1);
    global_dpd_->buf4_close(&Ht);
    global_dpd_->buf4_close(&Z1);

    /** W(a>b,ei) <---  0.5 * P(a/b) <nm||ef> * t1[m][b]*C[i][f]*t1[n][a] **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 12, 31, 12, 31, 0, "Z (m>n,ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D, &Cme, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 10, 31, 12, 31, 0, "Z (m>n,ei)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 31, 30, 31, 30, 0, "Z1 (ei,mb)");
    global_dpd_->contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 15, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    global_dpd_->contract244(&tia, &Z1, &Ht, 0, 2, 0, 1, 1);
    global_dpd_->buf4_close(&Ht);
    global_dpd_->buf4_close(&Z1);


    /** W(Ab,Ei) <---  <Nm|Ef> * t1[m][b]*C[i][f]*t1[N][A] **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 22, 26, 22, 26, 0, "Z (Nm,Ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &Cme, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 26, 24, 26, 24, 0, "Z1 (Ei,Nb)");
    global_dpd_->contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    global_dpd_->contract244(&tIA, &Z1, &Ht, 0, 2, 0, 1, 1);
    global_dpd_->buf4_close(&Ht);

    /** W(aB,eI) <---  <nM|eF> * t1[M][B]*C[I][F]*t1[n][a] **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 23, 25, 23, 25, 0, "Z (nM,eI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract424(&D, &CME, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 25, 27, 25, 27, 0, "Z1 (eI,nB)");
    global_dpd_->contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 29, 25, 29, 25, 0, "Ht_WaBeI (aB,eI)");
    global_dpd_->contract244(&tia, &Z1, &Ht, 0, 2, 0, 1, 1);
    global_dpd_->buf4_close(&Ht);


    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);


    /* combine W(Ab,Ei) and W(Ei,Ab)  */

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    global_dpd_->buf4_sort_axpy(&Ht, PSIF_CC_TMP0, rspq, 26, 28, "Ht_WAbEi (Ei,Ab)", 1);
    global_dpd_->buf4_close(&Ht);

    global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 29, 25, 29, 25, 0, "Ht_WaBeI (aB,eI)");
    global_dpd_->buf4_sort_axpy(&Ht, PSIF_CC_TMP0, rspq, 25, 29, "Ht_WaBeI (eI,aB)", 1);
    global_dpd_->buf4_close(&Ht);

    /* sort to Wabei (ei,ab) */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 7, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP2, rspq, 21, 7, "Ht_WABEI (EI,A>B)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 17, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP2, rspq, 31, 17, "Ht_Wabei (ei,a>b)");
    global_dpd_->buf4_close(&W);

    /* sort to Wabei (ie,ba) */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, C_irr, 21, 7, 21, 7, 0, "Ht_WABEI (EI,A>B)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1ET1, qprs, 20, 7, "Ht_WABEI (IE,B>A)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1ET1, C_irr, 20, 7, 20, 7, 0, "Ht_WABEI (IE,B>A)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP2, C_irr, 31, 17, 31, 17, 0, "Ht_Wabei (ei,a>b)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1ET1, qprs, 30, 17, "Ht_Wabei (ie,b>a)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1ET1, C_irr, 30, 17, 30, 17, 0, "Ht_Wabei (ie,b>a)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 26, 28, 26, 28, 0, "Ht_WAbEi (Ei,Ab)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1ET1, qpsr, 27, 29, "Ht_WAbEi (iE,bA)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 25, 29, 25, 29, 0, "Ht_WaBeI (eI,aB)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1ET1, qpsr, 24, 28, "Ht_WaBeI (Ie,Ba)");
    global_dpd_->buf4_close(&W);


    /************ TEST *************/
    /*
    dpd_buf4_init(&W, CC3_HC1ET1, 0, 20, 5, 20, 7, 0, "Ht_WABEI (IE,B>A)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WABEI(IE,B>A)|WABEI> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 30, 15, 30, 17, 0, "Ht_Wabei (ie,b>a)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wabei (ie,b>a)|Wabei> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 27, 29, 27, 29, 0, "Ht_WAbEi (iE,bA)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WAbEi (iE,bA)|WAbEi> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 24, 28, 24, 28, 0, "Ht_WaBeI (Ie,Ba)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WaBeI (Ie,Ba)|WaBeI> = %15.10lf\n", dot);
    */

  }

  return;
}

void HC1ET1_Wmbij_rhf(int i, int C_irr)
{
  double dot;
  dpdbuf4 C, D, E, F, Ht, W, W1, X, Z;
  dpdfile2 CME;
  char CME_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);

  global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);

  /***** Ht (Mb,Ij) <--- -WMnIj * Cnb *****/
  global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 10, 0, 10, 0, 0, "Ht_WMbIj (Mb,Ij)");
  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
  global_dpd_->contract424(&W, &CME, &Ht, 1, 0, 1, -1.0, 0.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Ht);

  /***** Ht (Mb,Ij) <--- CIE * WMbEj *****/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 0, 10, 0, 10, 0, "Z (MI,jb)");
  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
  global_dpd_->contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HC1ET1, psqr, 10, 0, "Ht_WMbIj (Mb,Ij)", 1);
  global_dpd_->buf4_close(&Z);

  /***** Ht (Mb,Ij) <--- Cje * WMbeI *****/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 0, 10, 0, 10, 0, "Z (Mj,Ib)");
  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
  global_dpd_->contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HC1ET1, psrq, 10, 0, "Ht_WMbIj (Mb,Ij)", -1);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Ht, PSIF_CC3_HC1ET1, C_irr, 10, 0, 10, 0, 0, "Ht_WMbIj (Mb,Ij)");
  global_dpd_->buf4_sort(&Ht, PSIF_CC3_HC1ET1, rspq, 0, 10, "Ht_WMbIj (Ij,Mb)");
  global_dpd_->buf4_close(&Ht);

  global_dpd_->file2_close(&CME);

  return;
}

void HC1ET1_Wabei_rhf(int i, int C_irr)
{
  double dot;
  dpdfile2 CME, tIA, T1;
  dpdbuf4 Ht, Z, Z1, Z2, Z3, B, C, D, E, F, W, X;
  char CME_lbl[32];
  int Gef, Gei, Gab, Ge, Gf, Gi;
  int EE, e;
  int nrows, ncols, nlinks;

  sprintf(CME_lbl, "%s %d", "CME", i);

  global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
  global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

  /***** Ht_WAbEi <--- -CMA * WMbEi *****/
  global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 11, 10, 11, 0, "CC3 WMbEj (Mb,Ej)");
  global_dpd_->contract244(&CME, &W, &Ht, 0, 0, 0, -1.0, 0.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Ht);

  /***** Ht_WAbEi <--- WAmEi * Cmb *****/
  global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "Ht_WAbEi (Ei,Ab)");
  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 11, 11, 11, 0, "CC3 WMbeJ (bM,eJ)");
  global_dpd_->contract424(&W, &CME, &Ht, 1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Ht);

  /***** HAbEi <--- <Ab|Ef> * Cif *****/
/*
  dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  dpd_contract424(&B, &CME, &Ht, 3, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&B);
  dpd_buf4_close(&Ht);
*/

  // Added new B(+)/B(-) code from cchbar 12/29/09, -TDC
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
  global_dpd_->file2_init(&T1, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 11, 8, 11, 8, 0, "Z1(ei,a>=b)");
  global_dpd_->buf4_scm(&Z1, 0.0); /* needed for empty occpi or virtpi irreps */
  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gab = Gef; /* B is totally symmetric */
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gf = Ge ^ Gef;
      Gi = Gf ^ C_irr;  /* T1 is not necessarily totally symmetric */
      Gei = Ge ^ Gi;
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gab]);
      Z1.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi],Z1.params->coltot[Gab]);

      nrows = moinfo.occpi[Gi];
      ncols = Z1.params->coltot[Gab];
      nlinks = moinfo.virtpi[Gf];
      if(nrows && ncols && nlinks) {
        for(EE=0; EE < moinfo.virtpi[Ge]; EE++) {
          e = moinfo.vir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
          C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,
                  B.matrix[Gef][0],ncols,0.0,Z1.matrix[Gei][0],ncols);
          global_dpd_->buf4_mat_irrep_wrt_block(&Z1,Gei,Z1.row_offset[Gei][e],moinfo.occpi[Gi]);
        }
      }
      global_dpd_->free_dpd_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gab]);
      global_dpd_->free_dpd_block(Z1.matrix[Gei], moinfo.occpi[Gi], Z1.params->coltot[Gab]);
    }
  }
  global_dpd_->buf4_close(&Z1);
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);

  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
  global_dpd_->file2_init(&T1, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 11, 9, 11, 9, 0, "Z2(ei,a>=b)");
  global_dpd_->buf4_scm(&Z2, 0.0); /* needed for empty occpi or virtpi irreps */
  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gab = Gef; /* W and B are totally symmetric */
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gf = Ge ^ Gef;
      Gi = Gf ^ C_irr;  /* T1 is not necessarily totally symmetric */
      Gei = Ge ^ Gi;
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gab]);
      Z2.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi],Z2.params->coltot[Gab]);

      nrows = moinfo.occpi[Gi];
      ncols = Z2.params->coltot[Gab];
      nlinks = moinfo.virtpi[Gf];
      if(nrows && ncols && nlinks) {
        for(EE=0; EE < moinfo.virtpi[Ge]; EE++) {
          e = moinfo.vir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
          C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,
                  B.matrix[Gef][0],ncols,0.0,Z2.matrix[Gei][0],ncols);
          global_dpd_->buf4_mat_irrep_wrt_block(&Z2, Gei, Z2.row_offset[Gei][e], moinfo.occpi[Gi]);
        }
      }
      global_dpd_->free_dpd_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gab]);
      global_dpd_->free_dpd_block(Z2.matrix[Gei], moinfo.occpi[Gi], Z2.params->coltot[Gab]);
    }
  }
  global_dpd_->buf4_close(&Z2);
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 11, 5, 11, 8, 0, "Z1(ei,a>=b)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 11, 5, 11, 9, 0, "Z2(ei,a>=b)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "Ht_WAbEi (Ei,Ab)");
  global_dpd_->buf4_axpy(&Z1, &W, 1);
  global_dpd_->buf4_axpy(&Z2, &W, 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  /** H(Ab,Ei) <--- -<Am|Ef>*C[i][f]*t1[m][b] **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 11, 11, 11, 11, 0, "Z(Am,Ei)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
  global_dpd_->contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
  global_dpd_->buf4_close(&F);

  global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "Ht_WAbEi (Ei,Ab)");
  global_dpd_->contract424(&Z, &tIA, &Ht, 1, 0, 0, -1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&Ht);

  /** Ht(Ab,Ei) <--- -<Mb|Ef>*C[i][f]*t1[M][A] **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
  global_dpd_->buf4_close(&F);

  global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  global_dpd_->contract244(&tIA, &Z, &Ht, 0, 0, 0, -1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&Ht);

  /** W(Ab,Ei) <---  <Nm|Ef> * t1[m][b]*C[i][f]*t1[N][A] **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 0, 11, 0, 11, 0, "Z (Nm,Ei)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract424(&D, &CME, &Z, 3, 1, 0, 1, 0);
  global_dpd_->buf4_close(&D);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 11, 10, 11, 10, 0, "Z1 (Ei,Nb)");
  global_dpd_->contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  global_dpd_->contract244(&tIA, &Z1, &Ht, 0, 2, 0, 1, 1);
  global_dpd_->buf4_close(&Ht);

  global_dpd_->file2_close(&CME);
  global_dpd_->file2_close(&tIA);

  /* combine W(Ab,Ei) and W(Ei,Ab)  */
  global_dpd_->buf4_init(&Ht, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  global_dpd_->buf4_sort_axpy(&Ht, PSIF_CC_TMP0, rspq, 11, 5, "Ht_WAbEi (Ei,Ab)", 1);
  global_dpd_->buf4_close(&Ht);

  /* sort to Wabei (ie,ba) */
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "Ht_WAbEi (Ei,Ab)");
  global_dpd_->buf4_sort(&W, PSIF_CC3_HC1ET1, qpsr, 10, 5, "Ht_WAbEi (iE,bA)");
  global_dpd_->buf4_close(&W);

  return;
}


}} // namespace psi::cceom
