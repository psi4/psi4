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

/* WABEI_UHF(): Computes all contributions to the ABEI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (EI,AB) ordering and is referred to on disk as "WEIAB".
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
** For the ABEI spin case, we evaluate these contractions with two
** target orderings, (AB,EI) and (EI,AB), depending on the term.
** After all terms have been evaluated, the (AB,EI) terms are sorted
** into (EI,AB) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void build_Z1_AAAA(void);

void WABEI_UHF(void)
{
  dpdfile2 FME, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /** W(EI,AB) <--- <EI||AB> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WEIAB");
  global_dpd_->buf4_close(&F);

  /** W(EI,AB) <--- - F_ME t_MI^AB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  global_dpd_->file2_mat_init(&FME);
  global_dpd_->file2_mat_rd(&FME);
  for(int Gei=0; Gei < moinfo.nirreps; Gei++) {
    int Gmi = Gei;
    int Gab = Gei;
    global_dpd_->buf4_mat_irrep_init(&T2,Gmi);
    global_dpd_->buf4_mat_irrep_rd(&T2,Gmi);
    int row=0;
    for(int Ge=0; Ge<moinfo.nirreps; Ge++){
      int Gm= Ge;
      int Gi= Gm ^ Gmi;
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.aoccpi[Gi],W.params->coltot[Gei]);
      int nrows = moinfo.aoccpi[Gm];
      int ncols = moinfo.aoccpi[Gi] * W.params->coltot[Gei];
      if(nrows && ncols){
        for(int EE=0; EE< moinfo.avirtpi[Ge]; EE++){
          int e = moinfo.avir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
          C_DGEMV('t',nrows,ncols, -1.0,&T2.matrix[Gmi][row][0],ncols,
              &FME.matrix[Gm][0][EE],moinfo.avirtpi[Ge], 1.0,W.matrix[Gei][0],1);
          global_dpd_->buf4_mat_irrep_wrt_block(&W,Gei,W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
        }
      }
      row+= moinfo.aoccpi[Gm]* moinfo.aoccpi[Gi];
      global_dpd_->free_dpd_block(W.matrix[Gei],moinfo.aoccpi[Gi],W.params->coltot[Gei]);
    }
    global_dpd_->buf4_mat_irrep_close(&T2,Gmi);

  }
  global_dpd_->file2_mat_close(&FME);
  global_dpd_->file2_close(&FME);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);

  /** W'(AB,EI) <--- <AB||EF> t_I^F **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 7, 5, 5, 1, "B <AB|CD>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  for(int Gef=0; Gef < moinfo.nirreps; Gef++) {
    int Gei = Gef;
    int Gab = Gef; /* W and B are totally symmetric */
    for(int Ge=0; Ge<moinfo.nirreps; Ge++){
      int Gf= Ge ^ Gef;
      int Gi= Gf;
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.avirtpi[Gf],B.params->coltot[Gef]);
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.aoccpi[Gi],W.params->coltot[Gei]);
      int nrows = moinfo.aoccpi[Gi];
      int ncols = W.params->coltot[Gei];
      int nlinks = moinfo.avirtpi[Gf];
      if(nrows && ncols){
        for(int EE=0; EE< moinfo.avirtpi[Ge]; EE++){
          int e = moinfo.avir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e],moinfo.avirtpi[Gf]);
          global_dpd_->buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
          C_DGEMM('n','n',nrows,ncols,nlinks,1.0,T1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
              1.0,W.matrix[Gei][0],ncols);
          global_dpd_->buf4_mat_irrep_wrt_block(&W, Gei,W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
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
   * W(EI,AB) <-- T1(I,F)D(EF,MN)Tau(MN,AB)
   * T1(I,F)D(EF,MN)Tau(MN,AB) ==  (W(MN,EI) Tau(MN,AB)
   * Notes:
   *      1. W_MNIE intermediate is read from disk (M>N-,EI)order to temp buffer Z
   *      2. W_MNIE is sorted to (EI,M>N-) order, Saved to disk, Re-Read into buffer Z
   *            in (EI, M>N-) order
   *      3. Tau_IJAB (MN,AB) is read from disk.
   *      5. Read W_ABEI (EI, A>B-) into buffer W.
   *      4. Loop over EI(row index) of W_EIAB target:
   */
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 2, 21, 2,21, 0, "WMNIE (M>N,EI)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rspq, 21, 2, "WMNIE (EI,M>N)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR,  0, 21, 7, 21, 7, 0, "WEIAB");
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR,  0, 21, 2, 21, 2, 0, "WMNIE (EI,M>N)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0,  2, 7,  2, 7, 0, "tauIJAB");
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
         C_DGEMV('t',nrows,ncols,-1,T.matrix[Gei][0],ncols,Z.matrix[Gei][ei],1,
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

  /* (T2+T1*T1)*F
   * <AI||BC> sorted F(AB,IC)
   * Build Z1(IA,MF) = t_{IM}^{AF} -t_I^Ft_M^A
   * contract F(AB,IC)Z1(IA,MF) => W1(BE,IA)
   * <Ai|Bc> sorted F(AB,ic)
   * t_Im^Af = t(IA,mf)
   * contract F(AB,ic)t(IA,mf) => W1(BE,IA)
   * sort W1(BE,IA) W2(EI,AB)
   * Read W2 with anti=1
   * Axpy W2=> WEIAB
   *
   **/
  build_Z1_AAAA();//Z(IA,MF)

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z1(IA,MF)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 5, 20, 5, 20, 0, "F <AI||BC> (AB,IC)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 20, 5, 20,  0, "W1(BE,IA)");
  global_dpd_->contract444(&F,&Z, &W, 0, 0, -1, 0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 5 ,30,"F <Ai|Bc> (AB,ic)");
  global_dpd_->buf4_close(&F);                                  //tIAmf
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0,"tIAjb");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS,  0, 5, 30, 5, 30,  0,"F <Ai|Bc> (AB,ic)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "W1(BE,IA)");
  global_dpd_->contract444(&F, &T, &W, 0, 0, -1, 1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "W1(BE,IA)");
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, qrsp, 21, 5, "W2(EI,AB)");
  global_dpd_->buf4_close(&Z1);

  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB" );
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 7, 21, 5, 1, "W2(EI,AB)");
  global_dpd_->buf4_axpy(&Z,&W,1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /* WEIAB <--- Pab t_M^A {<MB||EI> +t_IN^BF<MN||EF> -t_nI^fB<Mn|Ef>}
   * 1.   t_IN^BF * <MN||EF> + t_nI^fB * <Mn|Ef> --> Z(ME,IB)
   * 2.   Z(ME,IB) --sort--> Z(EI,MB)
   * 3. - t_M^A( <MB||EI> + Z(EI,MB) ) --> W'(EI,AB)
   * 4. WABEI <-- W'(EI,AB)- W'(EI,AB)
  **/
  /** Z(EI,BM) <-- -<IE||BM> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
  global_dpd_->buf4_sort_axpy(&C, PSIF_CC_TMP0, qprs, 21, 21, "Z(EI,BM)",-1);
  global_dpd_->buf4_close(&C);

  /** <MN||EF> t_IN^BF --> Z(ME,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <Mn|Ef> t_In^Bf --> Z(ME,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(ME,IB) --> Z(EI,BM) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, qrsp, 21, 21, "Z(EI,BM)",1);
  global_dpd_->buf4_close(&Z);

  /** Z(AB,EI) <-- -t_M^A Z(MB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 5, 21,  5, 0, "Z(EI,BA)");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 21, 21, 21, 21, 0, "Z(EI,BM)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z1, &T1, &Z, 3, 0, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&Z);

  /** WEIAB <--- Z(EI,AB) - Z(EI,BA)**/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 7, 21, 5, 1, "Z(EI,BA)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21,7, 21, 7, 0, "WEIAB" );
  global_dpd_->buf4_axpy(&Z,&W,-1.0);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_close(&W);

}

void build_Z1_AAAA(void)
{
  dpdbuf4 T2, Z1, Fint;
  dpdfile2 T1;
  int row,col,h;
  int GI,GB,GM,GF,GA;
  int i,b,m,f,a;
  int I,B,M,F,A;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "Z1(IA,MF)");
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z1(IA,MF)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  for(h = 0; h < moinfo.nirreps; h++){
    global_dpd_->buf4_mat_irrep_init(&Z1, h);
    global_dpd_->buf4_mat_irrep_rd(&Z1, h);
    for(row = 0; row< Z1.params->rowtot[h]; row ++){
      i  = Z1.params->roworb[h][row][0];
      a  = Z1.params->roworb[h][row][1];
      I  = T1.params->rowidx[i];
      A  = T1.params->colidx[a];
      GI = T1.params->psym[i];
      GA = T1.params->qsym[a];
      for(col =0; col < Z1.params->coltot[h]; col ++){
        m = Z1.params->colorb[h][col][0];
        f = Z1.params->colorb[h][col][1];
        M = T1.params->rowidx[m];
        F = T1.params->colidx[f];
        GM = T1.params->psym[m];
        GF = T1.params->qsym[f];

        if( GI == GF && GA == GM ){
          Z1.matrix[h][row][col] -= (T1.matrix[GI][I][F] * T1.matrix[GM][M][A]);
        }
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z1, h);
    global_dpd_->buf4_mat_irrep_close(&Z1, h);
  }
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
}//build_Z1_AAAA
}} // namespace psi::cchbar
