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

/* WAbEi_UHF(): Computes all contributions to the AbEi spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (Ei,Ab) ordering and is referred to on disk as "WEiAb".
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
** For the AbEi spin case, we evaluate these contractions with two
** target orderings, (Ab,Ei) and (Ei,Ab), depending on the term.
** After all terms have been evaluated, the (Ab,Ei) terms are sorted
** into (Ei,Ab) ordering and both groups arer added together.
**
** TDC, June 2002
*/
void build_Z1A_ABAB();
void build_Z1B_ABAB();

void WAbEi_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(Ei,Ab) <--- <Ei|Ab> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WEiAb");
  global_dpd_->buf4_close(&F);


  /** W(Ei,Ab) <--- - F_ME t_Mi^Ab **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1, 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /** W(Ei,Ab) <-- t(i,f) B(Ef,Ab) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  for(int Gef=0; Gef < moinfo.nirreps; Gef++) {
    int Gei = Gef;
    int Gab = Gef; /* W and B are totally symmetric */
    for(int Ge=0; Ge<moinfo.nirreps; Ge++){
      int Gf = Ge ^ Gef;
      int Gi = Gf;
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.bvirtpi[Gf],B.params->coltot[Gef]);
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.boccpi[Gi],W.params->coltot[Gei]);
      int nrows = moinfo.boccpi[Gi];
      int ncols = W.params->coltot[Gei];
      int nlinks = moinfo.bvirtpi[Gf];
      if(nrows && ncols){
        for(int EE=0; EE< moinfo.avirtpi[Ge]; EE++){
          int e = moinfo.avir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e],moinfo.bvirtpi[Gf]);
          global_dpd_->buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][e],moinfo.boccpi[Gi]);
          C_DGEMM('n','n',nrows,ncols,nlinks,1.0,T1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
              1.0,W.matrix[Gei][0],ncols);
          global_dpd_->buf4_mat_irrep_wrt_block(&W, Gei,W.row_offset[Gei][e],moinfo.boccpi[Gi]);
        }
      }
      global_dpd_->free_dpd_block(B.matrix[Gef],moinfo.bvirtpi[Gf],W.params->coltot[Gef]);
      global_dpd_->free_dpd_block(W.matrix[Gei],moinfo.boccpi[Gi],W.params->coltot[Gei]);
    }
  }
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);

  /* W(Ei,Ab) <-- T1(i,f)D(Ef,Mn)Tau(Mn,Ab)
   * T1(i,f)D(Ef,Mn)Tau(Mn,Ab)  == [W(Mn,Ei) Tau(Mn,Ab)]
   * Notes:
   *      1. W_mNiE intermediate is read from disk (mN,Ei)order to temp buffer Z
   *      2. W_mNiE is sorted to (Ei,Mn) order, Saved to disk, Re-Read into buffer Z
   *            in (Ei, Mn) order
   *      3. Tau_IJAB (Mn,Ab) is read from disk.
   *      5. Read W_AbEi (Ei, Ab) into buffer W.
   *      4. Loop over Ei(row index) of W_EiAb target:
   * --AMJ 1/16
   */
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 23, 26, 23,26, 0, "WmNiE (mN,Ei)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rsqp, 26, 22, "WmNiE (Ei,Mn)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR,  0, 26, 28, 26, 28, 0, "WEiAb");
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR,  0, 26, 22, 26, 22, 0, "WmNiE (Ei,Mn)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0,  22, 28,  22, 28, 0, "tauIjAb");
  for(int Gei=0; Gei< moinfo.nirreps; Gei++) {
    int Gab = Gei;
    int Gnm = Gei; /* Everything is totally symmetric */
    int nrows = T.params->rowtot[Gnm];
    int ncols = T.params->coltot[Gab];
    if (nrows && ncols) {
      global_dpd_->buf4_mat_irrep_init(&Z, Gei);
      global_dpd_->buf4_mat_irrep_rd(&Z, Gei);
      global_dpd_->buf4_mat_irrep_init(&T, Gnm);
      global_dpd_->buf4_mat_irrep_rd(&T, Gnm);
      global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
      for(int ei=0; ei< W.params->rowtot[Gei]; ei++){
         global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
         C_DGEMV('t',nrows,ncols,1,T.matrix[Gei][0],ncols,Z.matrix[Gei][ei],1,
             1.0,W.matrix[Gei][0],1);
         global_dpd_->buf4_mat_irrep_row_wrt(&W, Gei, ei);
      }
      global_dpd_->buf4_mat_irrep_row_close(&W, Gei);
      global_dpd_->buf4_mat_irrep_close(&T, Gnm);
      global_dpd_->buf4_mat_irrep_close(&Z, Gei);
    }
  }
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**
   * WAbEi <-- (T2+T1*T1) *F
   *
   * Z1a(iA,Mf) = t(i,f)t(M,A) - T(MA,if)
   * <Mb|Ef>(bE,Mf)Z1b(iA,Mf) = W1(bE,iA)           [Contract 444]
   * Z1b(ib,mf) = t(mb,if) - t(i,f)t(m,b)
   * <Am|Ef>(AE,mf)Z1b(ib,mf) = W2(AE,ib)           [Contract 444]
   * <AM||EF>(AE,MF)T2(ib,MF) =>W2(AE,ib)           [Contract 444]
   * W2(AE,ib) sort axpy(qrps) WAbEi (Ei,Ab)
   * W1(bE,iA) sort axpy(qrsp) WAbEi (Ei,Ab)
   *
   */

  /** Z1a = t_i^f t_M^A + T_Mi^Af **/
  build_Z1A_ABAB();

  /** W1(bE,iA)<--  <Mb|Ef>Z1a(iA,Mf) **/
  global_dpd_->buf4_init(&F,PSIF_CC_FINTS, 0, 29, 24, 29, 24, 0, "F <Ia|Bc> (aB,Ic)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 24, 27 ,24, 0, "Z1a(iA,Mf)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 27, 29, 27, 0, "W1(bE,iA)");
  global_dpd_->contract444(&F, &Z, &W, 0, 0, -1, 0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&F);

  /** Z1b = +t_mi^bf - t_i^f t_m^b **/
  build_Z1B_ABAB();

  /** W2(ib,AE)<-- <Am|Ef>Z1b(ib,mf) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0,  5, 30,  5, 30, 0, "F <Ai|Bc> (AB,ic)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0 , 0, 30, 30, 30, 30, 0, "Z1b(ib,mf)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0,  0,  5, 30,  5, 30, 0, "W2(AE,ib)");
  global_dpd_->contract444(&F, &Z, &W, 0, 0, 1, 0);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** W2(ib,AE)<--t_iM^bF<AM||EF> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0,  5, 20,  5, 20, 0, "F <AI||BC> (AB,IC)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0,  "tiaJB");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0,  0,  5, 30,  5, 30, 0, "W2(AE,ib)");
  global_dpd_->contract444(&F, &T, &W, 0,0, 1, 1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&W);

  /** Add W2 and W1 to target **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 27, 29, 27, 0, "W1(bE,iA)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, qrsp, 26, 28, "WEiAb",1 );
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0,  5, 30,  5, 30, 0,"W2(AE,ib)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, qrps, 26, 28, "WEiAb",1 );
  global_dpd_->buf4_close(&W);

  /** WAbEi <-- - t_M^A { <Mb|Ei> + t_iN^bF <MN||EF> + t_in^bf <Mn|Ef> }
                + t_m^b {-<mA|iE> + t_iN^fA <mN|fE> }
      Evaluate in three steps:
         (1) Z_MbEi =  <Mb|Ei> + t_iN^bF <MN||EF> + tin^bf <Mn|Ef>  [stored (Mb,Ei)]
         (2) Z_mAEi = -<mA|iE> + t_iN^fA <mN|fE>                    [stored (Am,Ei)]
         (3) WAbEi <-- - t_M^A Z_MbEi + t_m^b Z_mAEi
      Store targets in     W'(Ab,Ei) and  W(Ei,AB)
  **/
  /** Z(Mb,Ei) <-- <Mb|Ei> **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
  global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "Z(Mb,Ei)");
  global_dpd_->buf4_close(&D);

  /** <MN||EF> t_iN^bF --> Z(ME,ib) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <Mn|Ef> t_in^bf --> Z(ME,ib) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(ME,ib) --> Z(Mb,Ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, psqr, 24, 26, "Z(Mb,Ei)", 1);
  global_dpd_->buf4_close(&Z);

  /** W'(Ab,Ei) <-- - t_M^A Z(Mb,Ei) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(Am,Ei) <-- - <mA|iE> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
  global_dpd_->buf4_copy(&C, PSIF_CC_TMP0, "Z(Am,Ei)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** Z(mE,iA) <-- t_iN^fA <mN|fE> **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mE,iA)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(mE,iA) --> Z(Am,Ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mE,iA)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, spqr, 26, 26, "Z(Am,Ei)", 1);
  global_dpd_->buf4_close(&Z);

  /** W(Ei,AB) <-- t_m^b Z_mAEi **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Combine accumulated W'(Ab,Ei) and W(Ei,Ab) terms into WEiAb ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rspq, 26, 28, "WEiAb", 1);
  global_dpd_->buf4_close(&W);
}

void build_Z1A_ABAB(){
  dpdfile2 Tif, TMA;
  dpdbuf4 Z, T2;
  int i, f, M, A, iorb, forb, Aorb, Morb, Gi, GM, GA, Gf;
  int h, row, col;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb" );
  global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "Z1a(iA,Mf)");
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 24, 27, 24, 0, "Z1a(iA,Mf)");
  global_dpd_->file2_init(&TMA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&Tif, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&Tif);
  global_dpd_->file2_mat_init(&TMA);
  global_dpd_->file2_mat_rd(&TMA);
  global_dpd_->file2_mat_rd(&Tif);

  for(h = 0; h<moinfo.nirreps; h++){
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row = 0; row< Z.params->rowtot[h]; row++){
      iorb = Z.params->roworb[h][row][0];
      Aorb = Z.params->roworb[h][row][1];
      i = Tif.params->rowidx[iorb];
      A = TMA.params->colidx[Aorb];
      Gi = Tif.params->psym[iorb];
      GA = TMA.params->qsym[Aorb];
      for(col = 0; col < Z.params->coltot[h]; col++ ){
        Morb = Z.params->colorb[h][col][0];
        forb = Z.params->colorb[h][col][1];
        M = TMA.params->rowidx[Morb];
        f = Tif.params->colidx[forb];
        GM = TMA.params->psym[Morb];
        Gf = Tif.params->qsym[forb];

        if(Gi == Gf && GA == GM){
          Z.matrix[h][row][col] += Tif.matrix[Gi][i][f] * TMA.matrix[GM][M][A];
        }

      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->file2_mat_close(&Tif);
  global_dpd_->file2_mat_close(&TMA);
  global_dpd_->file2_close(&Tif);
  global_dpd_->file2_close(&TMA);
  global_dpd_->buf4_close(&Z);
}

void build_Z1B_ABAB(){
  dpdfile2 T1;
  dpdbuf4 Z, T2;
  int i, f, b, m, iorb, forb, borb, morb, Gi, Gm, Gb, Gf;
  int h, row, col;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb" );
  global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "Z1b(ib,mf)");
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z1b(ib,mf)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  for(h = 0; h < moinfo.nirreps; h++){
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row = 0; row < Z.params->rowtot[h]; row++){
      iorb = Z.params->roworb[h][row][0];
      borb = Z.params->roworb[h][row][1];
      i = T1.params->rowidx[iorb];
      b = T1.params->colidx[borb];
      Gi = T1.params->psym[iorb];
      Gb = T1.params->qsym[borb];
      for(col = 0; col < Z.params->coltot[h]; col++){
        morb = Z.params->colorb[h][col][0];
        forb = Z.params->colorb[h][col][1];
        m = T1.params->rowidx[morb];
        f = T1.params->colidx[forb];
        Gm = T1.params->psym[morb];
        Gf = T1.params->qsym[forb];

        if(Gi == Gf && Gb == Gm) {
           Z.matrix[h][row][col] -= T1.matrix[Gi][i][f] *T1.matrix[Gm][m][b];
        }
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
}
}} // namespace psi::cchbar
