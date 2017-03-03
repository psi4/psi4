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
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* build excited state -> excited state transition moment for oscillator strengths 
   and rotational strengths */

/* LR1_OO  =  LIE * RJE */
/* LR1_oo  = Lie * Rje */
/* LR1_VV = LMA * RMB */
/* LR1_vv = Lma * Rmb */
/* LT1_OO  =  LIE * TJE */
/* LT1_oo  = Lie * Tje */
/* LT1_VV = LMA * TMB */
/* LT1_vv = Lma * Tmb */
/* L2R1_OV = RME * LIMAE + Rme + LImAe */
/* L2R1_ov = Rme * Limae + RME + LiMaE */
/* L1R2_OV = LME * RIMAE + Lme * RImAe */
/* L1R2_ov = Lme * Rimae + LME * RiMaE */
/* L1T2_OV = LME * TIMAE + Lme * TImAe */
/* L1T2_ov = Lme * Timae + LME * TiMaE */
/* LR2_OO  = 0.5 * LIMEF * RJMEF + LImEf * RJmEf */
/* LR2_oo  = 0.5 * Limef * Rjmef + LiMeF * RjMeF */
/* LR2_VV = 0.5 * LMNEA * RMNEB + LmNeA * RmNeB */
/* LR2_vv = 0.5 * Lmnea * Rmneb + LMnEa * RMnEb */
/* LT2_OO = 0.5 * LIMEF * TJMEF + LImEf * TJmEf */
/* LT2_oo = 0.5 * Limef * Tjmef + LiMeF * TjMeF */
/* LT2_VV = 0.5 * LMNEA * TMNEB + LmNeA * TmNeB */
/* LT2_vv = 0.5 * Lmnea * Tmneb + LMnEa * TMnEb */
/* LR_OO = LR1_OO + LR2_OO */
/* LR_oo = LR1_oo + LR2_oo */
/* LR_VV = LR1_VV + LR2_VV */
/* LR_vv = LR1_vv + LR2_vv */
/* LT_OO = LT1_OO + LT2_OO */
/* LT_oo = LT1_oo + LT2_oo */

void ex_tdensity_intermediates(struct TD_Params S, struct TD_Params U)
{
  dpdfile2 L1, R1, T1, I, LR1, LR2, LT1, LT2;
  dpdbuf4 L2, T2, R2;
  int rohf = 0;

  /* Set LHS to Excited State */
  int LHS    = PSIF_CC_GL;
  /* Convenient Irrep Handles */
  int Lirrep = S.irrep;
  int Rirrep = U.irrep;
  int Tirrep = Lirrep^Rirrep;

  /* Auxiliary Irrep Handle */
  // -> T's (ground) with LHS (excited)
  // -> Should always be Lirrep (based on code structure)
  int Airrep = Lirrep;

  if ( (params.ref == 0) || (params.ref == 1) ) rohf = 1;

  /* LR1_OO(I,J)  =  LIE * RJE */
  global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 0, "LR1_OO");
  global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
  global_dpd_->contract222(&L1, &R1, &I, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&I);

  /* LR1_oo(i,j)  = Lia * Rje */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 0, "LR1_oo");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "Ria");
    global_dpd_->contract222(&L1, &R1, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 2, 2, "LR1_oo");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 2, 3, "Ria");
    global_dpd_->contract222(&L1, &R1, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&I);
  }

  /* LR1_VV(A,B) = LMA * RMB */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 1, 1, "LR1_VV");
  global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
  global_dpd_->contract222(&L1, &R1, &I, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&I);

  /* LR1_vv(a,b) = Lma * Rmb */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 1, 1, "LR1_vv");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "Ria");
    global_dpd_->contract222(&L1, &R1, &I, 1, 1, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 3, 3, "LR1_vv");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 2, 3, "Ria");
    global_dpd_->contract222(&L1, &R1, &I, 1, 1, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&I);
  }
  
  /* LT1_OO(I,J)  =  LIE * TJE */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 0, "LT1_OO");
  global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract222(&L1, &T1, &I, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&I);

  /* LT1_oo(i,j)  = Lia * Rje */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 0, "LT1_oo");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&L1, &T1, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 2, 2, "LT1_oo");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&L1, &T1, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&I);
  }

  /* LT1_VV(A,B) = LMA * TMB */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 1, 1, "LT1_VV");
  global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract222(&L1, &T1, &I, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&I);

  /* LT1_vv(a,b) = Lma * Tmb */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 1, 1, "LT1_vv");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&L1, &T1, &I, 1, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 3, 3, "LT1_vv");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&L1, &T1, &I, 1, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&I);
  }

  /* L2R1_OV(I,A) = RME * LIMAE + Rme + LImAe */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "Ria");
    global_dpd_->dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 2, 3, "Ria");
    global_dpd_->dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* L2R1_OV(i,a) = Rme * Limae + RME + LiMaE */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_ov");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 2, 7, 0, "Lijab");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "Ria");
    global_dpd_->dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 2, 3, "L2R1_ov");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 10, 15, 12, 17, 0, "Lijab");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 2, 3, "Ria");
    global_dpd_->dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* L1R2_OV(I,A) = LME * RIMAE + Lme * RImAe */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 1, "L1R2_OV");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 2, 7, 0, "RIJAB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 1, "L1R2_OV");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 2, 7, 0, "RIJAB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&I);
  }

  /* L1R2_ov(i,a) = Lme * Rimae + LME * RiMaE */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 1, "L1R2_ov");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 2, 7, 0, "Rijab");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 0, 5, 0, "RiJaB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 2, 3, "L1R2_ov");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 10, 15, 12, 17, 0, "Rijab");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 23, 29, 23, 29, 0, "RiJaB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&I);
  }

  /* L1T2_OV = LME * TIMAE + Lme * TImAe */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 1, "L1T2_OV");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &T2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 1, "L1T2_OV");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &T2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&I);
  }

  /* L1T2_ov = Lme * Timae + LME * TiMaE */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 1, "L1T2_ov");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &T2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 2, 3, "L1T2_ov");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &T2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&I);
  }

  /* LR2_OO(I,J)  = 0.5 * LIMEF * RJMEF + LImEf * RJmEf */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 0, "LR2_OO");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 7, 2, 7, 0, "RIJAB");
    global_dpd_->contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 0, "LR2_OO");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 7, 2, 7, 0, "RIJAB");
    global_dpd_->contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* LR2_oo(i,j)  = 0.5 * Limef * Rjmef + LiMeF * RjMeF */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 0, "LR2_oo");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 7, 2, 7, 0, "Rijab");
    global_dpd_->contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 0, 5, 0, "RiJaB");
    global_dpd_->contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 2, 2, "LR2_oo");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 10, 17, 12, 17, 0, "Rijab");
    global_dpd_->contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 23, 29, 23, 29, 0, "RiJaB");
    global_dpd_->contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* LR2_VV(A,B) = 0.5 * LMNEA * RMNEB + LmNeA * RmNeB */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 1, 1, "LR2_VV");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 2, 5, 2, 7, 0, "RIJAB");
    global_dpd_->contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 0, 5, 0, "RiJaB");
    global_dpd_->contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 1, 1, "LR2_VV");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 2, 5, 2, 7, 0, "RIJAB");
    global_dpd_->contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 23, 29, 23, 29, 0, "RiJaB");
    global_dpd_->contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* LR2_vv(a,b) = 0.5 * Lmnea * Rmneb + LMnEa * RMnEb */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 1, 1, "LR2_vv");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 2, 5, 2, 7, 0, "Rijab");
    global_dpd_->contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 3, 3, "LR2_vv");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 12, 15, 12, 17, 0, "Rijab");
    global_dpd_->contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, Rirrep, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* LT2_OO(I,J) = 0.5 * LIMEF * TJMEF + LImEf * TJmEf */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 0, "LT2_OO");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 0, "LT2_OO");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  
  /* LT2_oo(i,j) = 0.5 * Limef * Tjmef + LiMeF * TjMeF */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 0, "LT2_oo");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 2, 2, "LT2_oo");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* LT2_VV(A,B) = 0.5 * LMNEA * TMNEB + LmNeA * TmNeB */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 1, 1, "LT2_VV");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 1, 1, "LT2_VV");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* LT2_vv(a,b) = 0.5 * Lmnea * Tmneb + LMnEa * TMnEb */

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 1, 1, "LT2_vv");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 3, 3, "LT2_vv");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I);
  }

  /* LR_OO = LR1_OO + LR2_OO */
  /* LR_oo = LR1_oo + LR2_oo */
  /* LR_VV = LR1_VV + LR2_VV */
  /* LR_vv = LR1_vv + LR2_vv */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 0, "LR_OO");
  global_dpd_->file2_init(&LR1, PSIF_EOM_TMP, Tirrep, 0, 0, "LR1_OO");
  global_dpd_->file2_init(&LR2, PSIF_EOM_TMP, Tirrep, 0, 0, "LR2_OO");
  global_dpd_->file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
  global_dpd_->file2_close(&LR2);
  global_dpd_->file2_close(&LR1);
  global_dpd_->file2_close(&I);

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 0, 0, "LR_oo");
    global_dpd_->file2_init(&LR1, PSIF_EOM_TMP, Tirrep, 0, 0, "LR1_oo");
    global_dpd_->file2_init(&LR2, PSIF_EOM_TMP, Tirrep, 0, 0, "LR2_oo");
    global_dpd_->file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&LR2);
    global_dpd_->file2_close(&LR1);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 2, 2, "LR_oo");
    global_dpd_->file2_init(&LR1, PSIF_EOM_TMP, Tirrep, 2, 2, "LR1_oo");
    global_dpd_->file2_init(&LR2, PSIF_EOM_TMP, Tirrep, 2, 2, "LR2_oo");
    global_dpd_->file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&LR2);
    global_dpd_->file2_close(&LR1);
    global_dpd_->file2_close(&I);
  }

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 1, 1, "LR_VV");
  global_dpd_->file2_init(&LR1, PSIF_EOM_TMP, Tirrep, 1, 1, "LR1_VV");
  global_dpd_->file2_init(&LR2, PSIF_EOM_TMP, Tirrep, 1, 1, "LR2_VV");
  global_dpd_->file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
  global_dpd_->file2_close(&LR2);
  global_dpd_->file2_close(&LR1);
  global_dpd_->file2_close(&I);

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 1, 1, "LR_vv");
    global_dpd_->file2_init(&LR1, PSIF_EOM_TMP, Tirrep, 1, 1, "LR1_vv");
    global_dpd_->file2_init(&LR2, PSIF_EOM_TMP, Tirrep, 1, 1, "LR2_vv");
    global_dpd_->file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&LR2);
    global_dpd_->file2_close(&LR1);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Tirrep, 3, 3, "LR_vv");
    global_dpd_->file2_init(&LR1, PSIF_EOM_TMP, Tirrep, 3, 3, "LR1_vv");
    global_dpd_->file2_init(&LR2, PSIF_EOM_TMP, Tirrep, 3, 3, "LR2_vv");
    global_dpd_->file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&LR2);
    global_dpd_->file2_close(&LR1);
    global_dpd_->file2_close(&I);
  }

  /* LT_OO = LT1_OO + LT2_OO */
  /* LT_oo = LT1_oo + LT2_oo */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 0, "LT_OO");
  global_dpd_->file2_init(&LT1, PSIF_EOM_TMP, Airrep, 0, 0, "LT1_OO");
  global_dpd_->file2_init(&LT2, PSIF_EOM_TMP, Airrep, 0, 0, "LT2_OO");
  global_dpd_->file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
  global_dpd_->file2_close(&LT2);
  global_dpd_->file2_close(&LT1);
  global_dpd_->file2_close(&I);

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 0, 0, "LT_oo");
    global_dpd_->file2_init(&LT1, PSIF_EOM_TMP, Airrep, 0, 0, "LT1_oo");
    global_dpd_->file2_init(&LT2, PSIF_EOM_TMP, Airrep, 0, 0, "LT2_oo");
    global_dpd_->file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&LT2);
    global_dpd_->file2_close(&LT1);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 2, 2, "LT_oo");
    global_dpd_->file2_init(&LT1, PSIF_EOM_TMP, Airrep, 2, 2, "LT1_oo");
    global_dpd_->file2_init(&LT2, PSIF_EOM_TMP, Airrep, 2, 2, "LT2_oo");
    global_dpd_->file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&LT2);
    global_dpd_->file2_close(&LT1);
    global_dpd_->file2_close(&I);
  }

  /* LT_VV = LT1_VV + LT2_VV */
  /* LT_vv = LT1_vv + LT2_vv */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 1, 1, "LT_VV");
  global_dpd_->file2_init(&LT1, PSIF_EOM_TMP, Airrep, 1, 1, "LT1_VV");
  global_dpd_->file2_init(&LT2, PSIF_EOM_TMP, Airrep, 1, 1, "LT2_VV");
  global_dpd_->file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
  global_dpd_->file2_close(&LT2);
  global_dpd_->file2_close(&LT1);
  global_dpd_->file2_close(&I);

  if (rohf) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 1, 1, "LT_vv");
    global_dpd_->file2_init(&LT1, PSIF_EOM_TMP, Airrep, 1, 1, "LT1_vv");
    global_dpd_->file2_init(&LT2, PSIF_EOM_TMP, Airrep, 1, 1, "LT2_vv");
    global_dpd_->file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&LT2);
    global_dpd_->file2_close(&LT1);
    global_dpd_->file2_close(&I);
  }
  else {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, Airrep, 3, 3, "LT_vv");
    global_dpd_->file2_init(&LT1, PSIF_EOM_TMP, Airrep, 3, 3, "LT1_vv");
    global_dpd_->file2_init(&LT2, PSIF_EOM_TMP, Airrep, 3, 3, "LT2_vv");
    global_dpd_->file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&LT2);
    global_dpd_->file2_close(&LT1);
    global_dpd_->file2_close(&I);
  }

  return;
}

}} // namespace psi::ccdensity
