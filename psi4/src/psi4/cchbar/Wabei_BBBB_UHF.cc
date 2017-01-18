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

/* Wabei_UHF(): Computes all contributions to the abei spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (ei,ab) ordering and is referred to on disk as "Wabei".
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
** For the abei spin case, we evaluate these contractions with two
** target orderings, (ab,ei) and (ei,ab), depending on the term.
** After all terms have been evaluated, the (ab,ei) terms are sorted
** into (ei,ab) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void build_Z1_BBBB(void);

void Wabei_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /** W(ei,ab) <--- <ei||ab> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "Weiab");
  global_dpd_->buf4_close(&F);

  /** W(ei,ab) <--- - F_me t_mi^ab **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab");
  global_dpd_->file2_mat_init(&Fme);
  global_dpd_->file2_mat_rd(&Fme);
  for(int Gei=0; Gei < moinfo.nirreps; Gei++) {
    int Gmi = Gei;
    int Gab = Gei;
    global_dpd_->buf4_mat_irrep_init(&T2,Gmi);
    global_dpd_->buf4_mat_irrep_rd(&T2,Gmi);
    int row=0;
    for(int Ge=0; Ge<moinfo.nirreps; Ge++){
      int Gm= Ge;
      int Gi= Gm ^ Gmi;
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.boccpi[Gi],W.params->coltot[Gei]);
      int nrows = moinfo.boccpi[Gm];
      int ncols = moinfo.boccpi[Gi] * W.params->coltot[Gei];
      if(nrows && ncols){
        for(int EE=0; EE< moinfo.bvirtpi[Ge]; EE++){
          int e = moinfo.bvir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][e],moinfo.boccpi[Gi]);
          C_DGEMV('t',nrows,ncols, -1.0,&T2.matrix[Gmi][row][0],ncols,
              &Fme.matrix[Gm][0][EE],moinfo.bvirtpi[Ge], 1.0,W.matrix[Gei][0],1);
          global_dpd_->buf4_mat_irrep_wrt_block(&W,Gei,W.row_offset[Gei][e],moinfo.boccpi[Gi]);
        }
      }
      row+= moinfo.boccpi[Gm]* moinfo.boccpi[Gi];
      global_dpd_->free_dpd_block(W.matrix[Gei],moinfo.boccpi[Gi],W.params->coltot[Gei]);
    }
    global_dpd_->buf4_mat_irrep_close(&T2,Gmi);

  }
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);


  /** W(ab,ei) <-- T(i,f)B(ef,ab) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0,  31, 17, 31, 17, 0, "Weiab");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 15, 17, 15, 15, 1, "B <ab|cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  for(int Gef=0; Gef < moinfo.nirreps; Gef++) {
    int Gei = Gef;
    int Gab = Gef; /* W and B are totally symmetric */
    for(int Ge=0; Ge<moinfo.nirreps; Ge++){
      int Gf= Ge ^ Gef;
      int Gi= Gf;
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.bvirtpi[Gf],B.params->coltot[Gef]);
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.boccpi[Gi],W.params->coltot[Gei]);
      int nrows = moinfo.boccpi[Gi];
      int ncols = W.params->coltot[Gei];
      int nlinks = moinfo.bvirtpi[Gf];
      if(nrows && ncols){
        for(int EE=0; EE< moinfo.bvirtpi[Ge]; EE++){
          int e = moinfo.bvir_off[Ge] + EE;
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

  /*
   * T1(i,f)D(ef,mn)Tau(mn,ab)  == [W(mn,ei) Tau(mn,ab)]
   * 4 terms can be expressed as - (Tau_mn^ab W_mnei)
   * Notes:
   *      1. W_mnie intermediate is read from disk (m>n-,ei)order to temp buffer Z
   *      2. W_mnie is sorted to (EI,M>N-) order, Saved to disk, Re-Read into buffer Z
   *            in (EI, M>N-) order
   *      3. Tau_ijab (MN,AB) is read from disk.
   *      5. Read W_abei (ei, a>b-) into buffer W.
   *      4. Loop over ei(row index) of W_eiab target:
   * --AMJ 1/16
   */
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 12, 31, 12,31, 0, "Wmnie (m>n,ei)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rspq, 31, 12, "Wmnie (ei,m>n)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR,  0, 31, 17, 31, 17, 0, "Weiab");
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR,  0, 31, 12, 31, 12, 0, "Wmnie (ei,m>n)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0,  12, 17,  12, 17, 0, "tauijab");
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

  /*
   * <ai||bc> sorted F(ab,ic)
   * Build Z1(ia,mf) = t_{im}^{af} -t_i^ft_m^a
   * contract F(ab,ic)Z1(ia,mf) => W1(be,ia)
   * <aI|bC> sorted F(ab,IC)
   * t_iM^aF = t(ia,MF)
   * contract F(ab,IC)t(ia,MF) => W1(be,ia)
   * sort W1(be,ia) W2(ei,ab)
   * Read W2 with anti=1
   * Axpy W2=> Weiab
   *
   * --AMJ 02/19/2016
   **/
  build_Z1_BBBB();//Z(ia,mf)
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z1(ia,mf)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 30, 15, 30, 0, "F <ai||bc> (ab,ic)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 15, 30, 15, 30,  0, "W1(be,ia)");
  global_dpd_->contract444(&F,&Z, &W, 0, 0, -1, 0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 15 ,20,"F <aI|bC> (ab,IC)");
  global_dpd_->buf4_close(&F);                                  //tiaMF
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0,"tiaJB");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS,  0, 15, 20, 15, 20,  0,"F <aI|bC> (ab,IC)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 15, 30, 15, 30, 0, "W1(be,ia)");
  global_dpd_->contract444(&F, &T, &W, 0, 0, -1, 1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 15, 30, 15, 30, 0, "W1(be,ia)");
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, qrsp, 31, 15, "W2(ei,ab)");
  global_dpd_->buf4_close(&Z1);

  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab" );
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 31, 17, 31, 15, 1, "W2(ei,ab)");
  global_dpd_->buf4_axpy(&Z,&W,1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);


  /*- Pab t_m^a {<mb||ei> +t_in^bf<mn||ef> -t_Ni^Fb<mN|eF>}
   * 1.   t_in^bf * <mn||ef> + t_Ni^Fb * <mN|eF> --> Z(me,ib)
   * 2.   Z(me,ib) --sort--> Z(ei,mb)
   * 3. - t_m^a( <mb||ei> + Z(ei,mb) ) --> W'(ei,ab)
   * 4. Wabei <-- W'(ei,ab)- W'(ei,ab)
   */
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
  /* * a more efficient way to do scmcopy + sort
   * * sort_axpy into an empty buffer */
  global_dpd_->buf4_sort_axpy(&C, PSIF_CC_TMP0, qprs, 31, 31, "Z(ei,bm)",-1);
  global_dpd_->buf4_close(&C);

  /** <mn||ef> t_in^bf --> Z(me,ib) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(me,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <mN|eF> t_iN^bF --> Z(me,ib) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(me,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(me,ib) --> Z(mb,ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(me,ib)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, qrsp, 31, 31, "Z(ei,bm)",1);
  global_dpd_->buf4_close(&Z);

  /** Z(ei,ba) <-- -t_m^a Z(ei,bm) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 31, 15, 31,  15, 0, "Z(ei,ba)");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 31, 31, 31, 31, 0, "Z(ei,bm)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z1, &T1, &Z, 3, 0, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&Z);

  /** W(ei,ab) <-- Z(ei,ba) -  Z'(ei,ab) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 31, 17, 31, 15, 1, "Z(ei,ba)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 31,17, 31, 17, 0, "Weiab" );
  global_dpd_->buf4_axpy(&Z,&W,-1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

}

void build_Z1_BBBB(void)
{
  dpdbuf4 T2, Z1, Fint;
  dpdfile2 T1;
  int row,col,h;
  int GI,GB,GM,GF,GA;
  int i,b,m,f,a;
  int I,B,M,F,A;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "Z1(ia,mf)");
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z1(ia,mf)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
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
}
}} // namespace psi::cchbar
