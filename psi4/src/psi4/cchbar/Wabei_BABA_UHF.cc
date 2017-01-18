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
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* WaBeI_UHF(): Computes all contributions to the aBeI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (eI,aB) ordering and is referred to on disk as "WaBeI".
**
** The spin-orbital expression for the Wabei elements is:
**
** Wabei = <ab||ei> - Fme t_mi^ab + t_i^f <ab||ef>
**         - P(ab) t_m^b <am||ef>t_i^f + 1/2 tau_mn^ab <mn||ef> t_i^f
**         + 1/2 <mn||ei> tau_mn^ab - P(ab) <mb||ef> t_mi^af
**         - P(ab) t_m^a { <mb||ei> - t_ni^bf <mn||ef> }
**
** (cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).)
**
** For the aBeI spin case, we evaluate these contractions with two
** target orderings, (aB,eI) and (eI,aB), depending on the term.
** After all terms have been evaluated, the (aB,eI) terms are sorted
** into (eI,aB) ordering and both groups arer added together.
**
** TDC, June 2002
*/
void build_Z1A_BABA();

void WaBeI_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /** W(eI,aB) <--- <eI|aB> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WeIaB");
  global_dpd_->buf4_close(&F);

  /** W(eI,aB) <--- - F_me t_mI^aB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1, 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);


  /** Z(Ie,Ba) <--- t_I^F <Fe|Ba> **/
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 29,29,29,29, 0, "B <aB|cD>");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  for(int Gef=0; Gef < moinfo.nirreps; Gef++){
    int Gei = Gef;
    int Gab = Gef;
    for(int Ge=0; Ge<moinfo.nirreps; Ge++){
      int Gf = Ge ^ Gef;
      int Gi = Gf;
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.avirtpi[Gf], B.params->coltot[Gef]);
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.aoccpi[Gi],W.params->coltot[Gei]);
      int nrows = moinfo.aoccpi[Gi];
      int ncols = W.params->coltot[Gei];
      int nlinks = moinfo.avirtpi[Gf];
      if(nrows && ncols){
        for(int EE=0; EE < moinfo.bvirtpi[Ge]; EE++){
          int e = moinfo.bvir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&B,Gef,B.row_offset[Gef][e],moinfo.avirtpi[Gf]);
          global_dpd_->buf4_mat_irrep_rd_block(
              &W,Gei,W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
          C_DGEMM(
              'n','n',nrows,ncols,nlinks,
              1.0,T1.matrix[Gi][0],nlinks,
              B.matrix[Gef][0],ncols,
              1.0,W.matrix[Gei][0],ncols);
          global_dpd_->buf4_mat_irrep_wrt_block(
              &W,Gei,W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
        }
      }
      global_dpd_->free_dpd_block(B.matrix[Gef],moinfo.avirtpi[Gf],W.params->coltot[Gef]);
      global_dpd_->free_dpd_block(W.matrix[Gei],moinfo.aoccpi[Gi],W.params->coltot[Gei]);
    }
  }
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);

  /*
   * 4 terms can be expressed as - (Tau_Mn^Ab W_mNeI)
   * Notes:
   *      1. W_MnIe intermediate is read from disk W(Mn,eI)
   *      3. TauiJaB (nM,aB) is read from disk.
   *      4. tauiJaB is sorted (Mn,aB) order
   *      5. contract W(Mn,eI)tau(Mn,aB) --> W(eI,aB)
   * --AMJ 6/16
   */
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 22, 25, 22,25, 0, "WMnIe (Mn,eI)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR,  0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0,  23, 29,  23, 29, 0, "tauiJaB");
  global_dpd_->buf4_sort(&T, PSIF_CC_TMP0, qprs, 22,29, "tauJiaB");
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&T, PSIF_CC_TMP0, 0,22,29,22,29, 0, "tauJiaB");
  global_dpd_->contract444(&Z,&T,&W,1, 1, 1, 1);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**
   * WAbEi <-- (T2+T1*T1) *F
   *
   * Z1a(Ia,mF) = t(I,F)t(m,a) - T(ma,IF)
   * <mB|eF>(Be,mF)Z1a(Ia,mF) = W1(Be,Ia)            [Contract 444]
   * Z(aM,eI)<--  -<aM|eF> t_I^F                     [Contract 424]
   * W(eI,aB)<--  Z(aM,eI) t_M^B                     [Contract 424]
   * W'(ae,IB)<-- <am||ef>(ae,mf)T2(IB,mf)           [Contract 444]
   * W'(ae,IB)<-- <aM|eF>(ae,MF) t_IM^BF(MF,IB)      [Contract 444]
   * Z(Be,Ia) sort axpy(qrsp) WaBeI (eI,aB)
   *
   */

  build_Z1A_BABA();

    /** Z(Be,Ia)<--  <mB|eF>Z1a(Ia,mF) **/
  global_dpd_->buf4_init(&F,PSIF_CC_FINTS, 0, 28, 27, 28, 27, 0, "F <iA|bC> (Ab,iC)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 27, 24 ,27, 0, "Z1a(Ia,mF)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "Z(Be,Ia)");
  global_dpd_->contract444(&F, &Z, &W, 0, 0, 1, 0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&F);

    /** Z(aM,eI)<--  -<aM|eF> t_I^F) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25,25,25,25, 0, "Z'(aM,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0,0,1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, -1, 0);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);
  global_dpd_->file2_close(&T1);

    /** W(eI,aB)<-- Z(aM,eI) t_M^B **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z'(aM,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);
  global_dpd_->file2_close(&T1);

    /** W'(ae,IB)<-- t_Im^Bf<am||ef> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 30, 15, 30, 0, "F <ai||bc> (ab,ic)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0,  "tIAjb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0,  0, 15, 20, 15, 20, 0, "W'(ae,IB)");
  global_dpd_->contract444(&F, &T, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);

  /** W'(IB,ae) <-- <aM|eF>t_IM^BF **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0, "W'(ae,IB)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 20, 15, 20,0, "F <aI|bC> (ab,IC)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20,20,20,20,0, "tIAJB");
  global_dpd_->contract444(&F, &T, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);

  /** Add Z(Be,Ia) to target WaBeI **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "Z(Be,Ia)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, qrsp, 25, 29, "WeIaB",1 );
  global_dpd_->buf4_close(&W);

  /** WaBeI <-- - t_m^a { <mB|eI> + t_In^Bf <mn||ef> + t_IN^BF <mN|eF> }
                + t_M^B {-<Ma|Ie> + t_In^Fa <Mn|Fe> }
      Evaluate in three steps:
         (1) Z_mBeI =  <mB|eI> + t_In^Bf <mn||ef> + tIN^BF <mN|eF>
            stored (me,IB)
         (2) Z_MaeI = -<Ma|Ie> + t_In^Fa <Mn|Fe>
            [stored (aM,eI)]
         (3) WaBeI <-- - t_m^a Z_mBeI + t_M^B Z_MaeI
      Store targets in     W'(ae,IB) and  W(eI,aB)
  **/

  /** Z(mB,eI) <-- <mB|eI> **/
  global_dpd_->buf4_init(
    &D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
  global_dpd_->buf4_sort(&D, PSIF_CC_TMP0, prsq, 30,20, "Z(me,IB)");
  global_dpd_->buf4_close(&D);

  /** <mn||ef> t_In^Bf --> Z(me,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_init(
    &D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(
    &T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <mN|eF> t_IN^BF --> Z(me,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_init(
    &D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** W'(ae,IB) <-- - t_m^a Z(me,IB) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0, "W'(ae,IB)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(aM,eI) <-- - <Ma|Ie> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_TMP0, qpsr, 25, 25, "Z(aM,eI)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** Z(Me,Ia) <-- t_In^Fa <Mn|Fe> **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(Me,Ia) --> Z(aM,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, spqr, 25, 25, "Z(aM,eI)", 1);
  global_dpd_->buf4_close(&Z);

  /** W(eI,aB) <-- t_M^B Z_MaeI **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Combine accumulated W'(ae,IB) and W(eI,aB) terms into WeIaB ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0,"W'(ae,IB)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, qrps, 25, 29, "WeIaB",1 );
  global_dpd_->buf4_close(&W);
}

void build_Z1A_BABA()
{
  dpdfile2 TIF, Tma;
  dpdbuf4 Z, T2;
  int I, F, m, a, Iorb, Forb, aorb, morb, GI, Gm, Ga, GF;
  int h, row, col;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA" );
  global_dpd_->buf4_scmcopy(&T2, PSIF_CC_TMP0, "Z1a(Ia,mF)",-1);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 27, 24, 27, 0, "Z1a(Ia,mF)");
  global_dpd_->file2_init(&TIF, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&Tma, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&TIF);
  global_dpd_->file2_mat_init(&Tma);
  global_dpd_->file2_mat_rd(&Tma);
  global_dpd_->file2_mat_rd(&TIF);

  for(h = 0; h<moinfo.nirreps; h++){
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row = 0; row< Z.params->rowtot[h]; row++){
      Iorb = Z.params->roworb[h][row][0];
      aorb = Z.params->roworb[h][row][1];
      I = TIF.params->rowidx[Iorb];
      a = Tma.params->colidx[aorb];
      GI = TIF.params->psym[Iorb];
      Ga = Tma.params->qsym[aorb];
      for(col = 0; col < Z.params->coltot[h]; col++ ){
        morb = Z.params->colorb[h][col][0];
        Forb = Z.params->colorb[h][col][1];
        m = Tma.params->rowidx[morb];
        F = TIF.params->colidx[Forb];
        Gm = Tma.params->psym[morb];
        GF = TIF.params->qsym[Forb];

        if(GI == GF && Ga == Gm){
          Z.matrix[h][row][col] -= TIF.matrix[GI][I][F] * Tma.matrix[Gm][m][a];
        }

      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->file2_mat_close(&TIF);
  global_dpd_->file2_mat_close(&Tma);
  global_dpd_->file2_close(&TIF);
  global_dpd_->file2_close(&Tma);
  global_dpd_->buf4_close(&Z);
}
}} // namespace psi::cchbar
