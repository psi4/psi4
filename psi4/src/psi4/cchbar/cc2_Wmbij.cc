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
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

void purge_cc2_Wmbij(void);

void cc2_Wmbij_build(void)
{

  dpdfile2 t1, tIA, tia;
  dpdbuf4 E, D, C, F;
  dpdbuf4 W, W1;
  dpdbuf4 Z, Z1;
  dpdbuf4 X;

  if(params.ref == 0) { /** RHF **/
    /* W(Mb,Ij) <-- <Mb|Ij> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
    global_dpd_->buf4_copy(&E, PSIF_CC2_HET1, "CC2 WMbIj");
    global_dpd_->buf4_close(&E);
  }
  else if(params.ref == 1) { /** ROHF **/
    /** W(MB,I>J) <--- <MB||IJ> **/
    /** W(mb,i>j) <--- <mb||ij> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC2_HET1, rspq, 10, 2, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_sort(&E, PSIF_CC2_HET1, rspq, 10, 2, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_close(&E);

    /** W(Mb,Ij) <--- <Mb|Ij> **/
    /** W(mB,iJ) <--- <mB|iJ> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC2_HET1, rspq, 10, 0, "CC2 WMbIj");
    global_dpd_->buf4_sort(&E, PSIF_CC2_HET1, rspq, 10, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_close(&E);
  }
  else if(params.ref == 2) { /*** UHF ***/
    /** W(MB,I>J) <--- <MB||IJ> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_sort(&E, PSIF_CC2_HET1, rspq, 20, 2, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_close(&E);

    /** W(mb,i>j) <--- <mb||ij> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC2_HET1, rspq, 30, 12, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_close(&E);

    /** W(Mb,Ij) <--- <Mb|Ij> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC2_HET1, rspq, 24, 22, "CC2 WMbIj");
    global_dpd_->buf4_close(&E);

    /** W(mB,iJ) <--- <mB|iJ> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->buf4_sort(&E, PSIF_CC2_HET1, rspq, 27, 23, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_close(&E);
  }

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /* W(Mb,Ij) <-- -t(n,b) * W(Mn,Ij) */
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&W1, PSIF_CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
    global_dpd_->contract424(&W1, &t1, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&W1);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&t1);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /**** W(MB,I>J) <-- -ZMBJI <-- P(I/J)( -<JE||MB> * t1[I][E] ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z (MB,JI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tIA, &Z, 1, 1, 0, -1, 0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 10, 0, "X (MB,IJ)");
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "X (MB,IJ)");
    global_dpd_->buf4_axpy(&Z, &X, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_axpy(&X, &W, 1);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_close(&W);

    /**** W(mb,i>j) <-- -Zmbji <-- P(i/j)( -<je||mb> * t1[i][e] ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z (mb,ji)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tia, &Z, 1, 1, 0, -1, 0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 10, 0, "X (mb,ij)");
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "X (mb,ij)");
    global_dpd_->buf4_axpy(&Z, &X, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 2, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_axpy(&X, &W, 1);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_close(&W);

    /**** W(Mb,Ij) <-- ( <Mb|Ej> * t1[I][E] ) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &W, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /**** W(Mb,Ij) <-- ( <Mb|Ie> * t1[j][e] ) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_close(&W);

    /**** W(mB,iJ) <-- ( <mB|eJ> * t1[i][e] ) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tia, &D, &W, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /**** W(mB,iJ) <-- ( <mB|iE> * t1[J][E] ) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params.ref == 2) { /*** UHF ***/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /**** W(MB,I>J) <-- -ZMBJI <-- P(I/J)( -<JE||MB> * t1[I][E] ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Z (MB,JI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->contract424(&C, &tIA, &Z, 1, 1, 0, -1, 0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 20, 0, "X (MB,IJ)");
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "X (MB,IJ)");
    global_dpd_->buf4_axpy(&Z, &X, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 20, 0, 20, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_axpy(&X, &W, 1);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_close(&W);

    /**** W(mb,i>j) <-- -Zmbji <-- P(i/j)( -<je||mb> * t1[i][e] ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Z (mb,ji)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tia, &Z, 1, 1, 0, -1, 0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 30, 10, "X (mb,ij)");
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "X (mb,ij)");
    global_dpd_->buf4_axpy(&Z, &X, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 30, 10, 30, 12, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_axpy(&X, &W, 1);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_close(&W);

    /**** W(Mb,Ij) <-- ( <Mb|Ej> * t1[I][E] ) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 24, 22, 24, 22, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->contract244(&tIA, &D, &W, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /**** W(Mb,Ij) <-- ( <Mb|Ie> * t1[j][e] ) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 24, 22, 24, 22, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->contract424(&C, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_close(&W);

    /**** W(mB,iJ) <-- ( <mB|eJ> * t1[i][e] ) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 27, 23, 27, 23, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->contract244(&tia, &D, &W, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /**** W(mB,iJ) <-- ( <mB|iE> * t1[J][E] ) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 27, 23, 27, 23, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    global_dpd_->contract424(&C, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /* W(Mb,Ij) <-- +P(ij) t(I,E) * <Mb||Ej> */
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &t1, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "CC2 ZMbIj (jM,bI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract424(&D, &t1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC2_HET1, qrsp, 10, 0, "CC2 WMbIj", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&t1);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /**** W(MB,I>J) <-- ( t1[N][B] * W(MN,I>J) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 2, 10, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_init(&Z, PSIF_CC2_HET1, 0, 0, 2, 2, 2, 0, "CC2 WMNIJ (M>N,I>J)");
    global_dpd_->contract424(&Z, &tIA, &W, 1, 0, 1, -0.5, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(mb,i>j) <-- ( t1[n][b] * W(mn,i>j) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 2, 10, 2, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_init(&Z, PSIF_CC2_HET1, 0, 0, 2, 2, 2, 0, "CC2 Wmnij (m>n,i>j)");
    global_dpd_->contract424(&Z, &tia, &W, 1, 0, 1, -0.5, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(Mb,Ij) <-- ( t1[n][b] * W(Mn,Ij) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&Z, PSIF_CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
    global_dpd_->contract424(&Z, &tia, &W, 1, 0, 1, -0.5, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(mB,iJ) <-- ( t1[N][B] * W(mN,iJ) ****/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 0, 11, 0, 0, "Z (Bm,Ji)");
    global_dpd_->buf4_init(&Z, PSIF_CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
    global_dpd_->contract244(&tIA, &Z, &Z1, 0, 0, 0, -0.5, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC2_HET1, qpsr, 10, 0, "CC2 WmBiJ (mB,iJ)", 1);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params.ref == 2) { /*** UHF ***/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /**** W(MB,I>J) <-- ( t1[N][B] * W(MN,I>J) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 20, 2, 20, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_init(&Z, PSIF_CC2_HET1, 0, 0, 2, 2, 2, 0, "CC2 WMNIJ (M>N,I>J)");
    global_dpd_->contract424(&Z, &tIA, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(mb,i>j) <-- ( t1[n][b] * W(mn,i>j) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 30, 12, 30, 12, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_init(&Z, PSIF_CC2_HET1, 0, 10, 12, 12, 12, 0, "CC2 Wmnij (m>n,i>j)");
    global_dpd_->contract424(&Z, &tia, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(Mb,Ij) <-- ( t1[n][b] * W(Mn,Ij) ****/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 24, 22, 24, 22, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&Z, PSIF_CC2_HET1, 0, 22, 22, 22, 22, 0, "CC2 WMnIj");
    global_dpd_->contract424(&Z, &tia, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(mB,iJ) <-- ( t1[N][B] * W(mN,iJ) ****/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 26, 22, 26, 22, 0, "Z (Bm,Ji)");
    global_dpd_->buf4_init(&Z, PSIF_CC2_HET1, 0, 22, 22, 22, 22, 0, "CC2 WMnIj");
    global_dpd_->contract244(&tIA, &Z, &Z1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC2_HET1, qpsr, 27, 23, "CC2 WmBiJ (mB,iJ)", 1);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /* Z(Mb,Ej) = <Mb|Ej> + t(j,f) * <Mb|Ef> */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj (Mb,Ej)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract424(&F, &t1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    /* W(Mb,Ij) <-- t(I,E) * Z(Mb,Ej) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "CC2 ZMbIj");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj (Mb,Ej)");
    global_dpd_->contract244(&t1, &Z1, &Z, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&t1);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /**** W(MB,I>J) <-- 0.5*P(I/J)XMBIJ <--- ( t1[I][E] * ZMBEJ ) <--  <MB||EF> * t1[J][F] ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z (MB,EJ)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->contract244(&tIA, &Z, &W, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(mb,i>j) <-- P(i/j) (Zmbif * t1[j][f]) <-- 0.5*( t1[i][e] * <mb||ef> ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z (mb,ej)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 2, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->contract244(&tia, &Z, &W, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(Mb,Ij) <-- (ZIfMb * t1[j][f]) <--  t1[I][E] * <Mb|Ef> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (If,Mb)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract244(&tIA, &F, &Z, 1, 2, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->contract424(&Z, &tia, &W, 1, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(mB,iJ) <-- ZiFmB * t1[J][F] <-- t1[i][e] * <mB|eF> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (iF,mB)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract244(&tia, &F, &Z, 1, 2, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->contract424(&Z, &tIA, &W, 1, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params.ref == 2) { /*** UHF ***/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /**** W(MB,I>J) <-- 0.5*P(I/J)XMBIJ <--- ( t1[I][E] * ZMBEJ ) <--  <MB||EF> * t1[J][F] ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z (MB,EJ)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 0, 20, 2, 0, "CC2 ZMBIJ (MB,I>J)");
    global_dpd_->contract244(&tIA, &Z, &W, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 0, 20, 2, 0, "CC2 ZMBIJ (MB,I>J)");
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 20, 0, 20, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(mb,i>j) <-- P(i/j) (Zmbif * t1[j][f]) <-- 0.5*( t1[i][e] * <mb||ef> ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Z (mb,ej)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 10, 30, 12, 0, "CC2 Zmbij (mb,i>j)");
    global_dpd_->contract244(&tia, &Z, &W, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 10, 30, 12, 0, "CC2 Zmbij (mb,i>j)");
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 30, 10, 30, 12, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(Mb,Ij) <-- (ZIfMb * t1[j][f]) <--  t1[I][E] * <Mb|Ef> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z (If,Mb)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract244(&tIA, &F, &Z, 1, 2, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "CC2 ZMbIj");
    global_dpd_->contract424(&Z, &tia, &W, 1, 1, 0, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 24, 22, 24, 22, 0, "CC2 WMbIj");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "CC2 ZMbIj");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(mB,iJ) <-- ZiFmB * t1[J][F] <-- t1[i][e] * <mB|eF> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "Z (iF,mB)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract244(&tia, &F, &Z, 1, 2, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "CC2 ZmBiJ (mB,iJ)");
    global_dpd_->contract424(&Z, &tIA, &W, 1, 1, 0, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 27, 23, 27, 23, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "CC2 ZmBiJ (mB,iJ)");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }

  if(params.ref == 1) { /** ROHF **/

    /* do purge before sort */
    purge_cc2_Wmbij();

    /* do final sort to get (Ij,Mb) */
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 2, 10, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 2, 10, "CC2 WMBIJ (I>J,MB)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 2, 10, 2, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 2, 10, "CC2 Wmbij (i>j,mb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 0, 10, "CC2 WMbIj (Ij,Mb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 0, 10, "CC2 WmBiJ (iJ,mB)");
    global_dpd_->buf4_close(&W);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 20, 2, 20, 2, 0, "CC2 WMBIJ (MB,I>J)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 2, 20, "CC2 WMBIJ (I>J,MB)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 30, 12, 30, 12, 0, "CC2 Wmbij (mb,i>j)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 12, 30, "CC2 Wmbij (i>j,mb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 24, 22, 24, 22, 0, "CC2 WMbIj");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 22, 24, "CC2 WMbIj (Ij,Mb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 27, 23, 27, 23, 0, "CC2 WmBiJ (mB,iJ)");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, rspq, 23, 27, "CC2 WmBiJ (iJ,mB)");
    global_dpd_->buf4_close(&W);
  }
}

void purge_cc2_Wmbij(void) {
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

  global_dpd_->file4_init(&W, PSIF_CC2_HET1, 0, 10, 2,"CC2 WMBIJ (MB,I>J)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      b = W.params->roworb[h][mb][1];
      bsym = W.params->qsym[b];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
        if (B >= (virtpi[bsym] - openpi[bsym]))
          W.matrix[h][mb][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC2_HET1, 0, 10, 2,"CC2 Wmbij (mb,i>j)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      m = W.params->roworb[h][mb][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
        i = W.params->colorb[h][ij][0];
        j = W.params->colorb[h][ij][1];
        isym = W.params->rsym[i];
        jsym = W.params->ssym[j];
        I = i - occ_off[isym];
        J = j - occ_off[jsym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
            (I >= (occpi[isym] - openpi[isym])) ||
            (J >= (occpi[jsym] - openpi[jsym])) )
          W.matrix[h][mb][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC2_HET1, 0, 10, 0,"CC2 WMbIj");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      for(ij=0; ij<W.params->coltot[h]; ij++) {
        j = W.params->colorb[h][ij][1];
        jsym = W.params->ssym[j];
        J = j - occ_off[jsym];
        if (J >= (occpi[jsym] - openpi[jsym]))
          W.matrix[h][mb][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC2_HET1, 0, 10, 0,"CC2 WmBiJ (mB,iJ)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      m = W.params->roworb[h][mb][0];
      b = W.params->roworb[h][mb][1];
      msym = W.params->psym[m];
      bsym = W.params->qsym[b];
      M = m - occ_off[msym];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
        i = W.params->colorb[h][ij][0];
        isym = W.params->rsym[i];
        I = i - occ_off[isym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
            (B >= (virtpi[bsym] - openpi[bsym])) ||
            (I >= (occpi[isym] - openpi[isym])) )
          W.matrix[h][mb][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);
}

}} // namespace psi::cchbar
