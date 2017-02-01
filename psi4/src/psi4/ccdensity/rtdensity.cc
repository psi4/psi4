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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void rtdensity(struct TD_Params S)
{
  dpdfile2 D, T1, L1, Z;
  dpdbuf4 T2, L2;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 0, 0, "RTDIJ");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 0, 0, "RTDij");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "Lia");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 1, 1, "RTDAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 1, 1, "RTDab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 0, 1, "RTDIA");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 0, 0, "Z(I,M)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);

    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */

    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 0, 0, "Z(I,M)");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);

    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */

    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 1, 1, "Z(A,E)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 0, 1, "RTDia");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 0, 0, "Z(i,m)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 0, 0, "Z(i,m)");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 1, 1, "Z(a,e)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_TMP, "RTDAI");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_TMP, "RTDai");
    global_dpd_->file2_close(&L1);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 0, 0, "RTDIJ");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 2, 2, "RTDij");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 2, 3, "Lia");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 1, 1, "RTDAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 3, 3, "RTDab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 0, 1, "RTDIA");

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 0, 0, "Z(I,M)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);

    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */

    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 0, 0, "Z(I,M)");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);

    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */

    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 1, 1, "Z(A,E)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, S.irrep, 2, 3, "RTDia");

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 2, 2, "Z(i,m)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 2, 2, "Z(i,m)");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, S.irrep, 3, 3, "Z(a,e)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_TMP, "RTDAI");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_GL, S.irrep, 2, 3, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_TMP, "RTDai");
    global_dpd_->file2_close(&L1);
  }

  return;
}

}} // namespace psi::ccdensity
