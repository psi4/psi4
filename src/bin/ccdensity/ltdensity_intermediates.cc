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

/* build left transition moment for oscillator strengths 
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

void ltdensity_intermediates(struct TD_Params S)
{
  dpdfile2 L1, R1, T1, I, LR1, LR2, LT1, LT2;
  dpdbuf4 L2, T2, R2;
  int rohf = 0;

  if ( (params.ref == 0) || (params.ref == 1) ) rohf = 1;

  /* LR1_OO(I,J)  =  LIE * RJE */
  dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 0, "LR1_OO");
  dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
  dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "RIA");
  dpd_contract222(&L1, &R1, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);
  dpd_file2_close(&I);

  /* LR1_oo(i,j)  = Lia * Rje */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 0, "LR1_oo");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "Ria");
    dpd_contract222(&L1, &R1, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&L1);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 2, 2, "LR1_oo");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&R1, CC_GR, S.irrep, 2, 3, "Ria");
    dpd_contract222(&L1, &R1, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&L1);
    dpd_file2_close(&I);
  }

  /* LR1_VV(A,B) = LMA * RMB */

  dpd_file2_init(&I, EOM_TMP, S.irrep, 1, 1, "LR1_VV");
  dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
  dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "RIA");
  dpd_contract222(&L1, &R1, &I, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);
  dpd_file2_close(&I);

  /* LR1_vv(a,b) = Lma * Rmb */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 1, 1, "LR1_vv");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "Ria");
    dpd_contract222(&L1, &R1, &I, 1, 1, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&L1);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 3, 3, "LR1_vv");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&R1, CC_GR, S.irrep, 2, 3, "Ria");
    dpd_contract222(&L1, &R1, &I, 1, 1, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&L1);
    dpd_file2_close(&I);
  }
  
  /* LT1_OO(I,J)  =  LIE * TJE */

  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT1_OO");
  dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&L1, &T1, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&L1);
  dpd_file2_close(&I);

  /* LT1_oo(i,j)  = Lia * Rje */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT1_oo");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&L1, &T1, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&L1);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 2, 2, "LT1_oo");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&L1, &T1, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&L1);
    dpd_file2_close(&I);
  }

  /* LT1_VV(A,B) = LMA * TMB */

  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT1_VV");
  dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&L1, &T1, &I, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&L1);
  dpd_file2_close(&I);

  /* LT1_vv(a,b) = Lma * Tmb */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT1_vv");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&L1, &T1, &I, 1, 1, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&L1);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 3, 3, "LT1_vv");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&L1, &T1, &I, 1, 1, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&L1);
    dpd_file2_close(&I);
  }

  /* L2R1_OV(I,A) = RME * LIMAE + Rme + LImAe */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "Ria");
    dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_file2_init(&R1, CC_GR, S.irrep, 2, 3, "Ria");
    dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* L2R1_OV(i,a) = Rme * Limae + RME + LiMaE */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 2, 7, 0, "Lijab");
    dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "Ria");
    dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 15, 12, 17, 0, "Lijab");
    dpd_file2_init(&R1, CC_GR, S.irrep, 2, 3, "Ria");
    dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* L1R2_OV(I,A) = LME * RIMAE + Lme * RImAe */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 1, "L1R2_OV");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 2, 7, 0, "RIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 1, "L1R2_OV");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 2, 7, 0, "RIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, S.irrep, 22, 28, 22, 28, 0, "RIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  }

  /* L1R2_ov(i,a) = Lme * Rimae + LME * RiMaE */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 1, "L1R2_ov");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 2, 7, 0, "Rijab");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 2, 3, "L1R2_ov");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 10, 15, 12, 17, 0, "Rijab");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, S.irrep, 23, 29, 23, 29, 0, "RiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  }

  /* L1T2_OV = LME * TIMAE + Lme * TImAe */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L1T2_OV");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L1T2_OV");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_close(&I);
  }

  /* L1T2_ov = Lme * Timae + LME * TiMaE */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L1T2_ov");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 2, 3, "L1T2_ov");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_close(&I);
  }

  /* LR2_OO(I,J)  = 0.5 * LIMEF * RJMEF + LImEf * RJmEf */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 0, "LR2_OO");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 7, 2, 7, 0, "RIJAB");
    dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RIjAb");
    dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 0, "LR2_OO");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 7, 2, 7, 0, "RIJAB");
    dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 22, 28, 22, 28, 0, "RIjAb");
    dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* LR2_oo(i,j)  = 0.5 * Limef * Rjmef + LiMeF * RjMeF */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 0, "LR2_oo");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 7, 2, 7, 0, "Rijab");
    dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RiJaB");
    dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 2, 2, "LR2_oo");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 10, 17, 12, 17, 0, "Rijab");
    dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 23, 29, 23, 29, 0, "RiJaB");
    dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* LR2_VV(A,B) = 0.5 * LMNEA * RMNEB + LmNeA * RmNeB */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 1, 1, "LR2_VV");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 2, 5, 2, 7, 0, "RIJAB");
    dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RiJaB");
    dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 1, 1, "LR2_VV");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 2, 5, 2, 7, 0, "RIJAB");
    dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 23, 29, 23, 29, 0, "RiJaB");
    dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* LR2_vv(a,b) = 0.5 * Lmnea * Rmneb + LMnEa * RMnEb */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 1, 1, "LR2_vv");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 2, 5, 2, 7, 0, "Rijab");
    dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RIjAb");
    dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 3, 3, "LR2_vv");
    dpd_buf4_init(&L2, CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 12, 15, 12, 17, 0, "Rijab");
    dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&R2, CC_GR, S.irrep, 22, 28, 22, 28, 0, "RIjAb");
    dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* LT2_OO(I,J) = 0.5 * LIMEF * TJMEF + LImEf * TJmEf */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT2_OO");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT2_OO");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  
  /* LT2_oo(i,j) = 0.5 * Limef * Tjmef + LiMeF * TjMeF */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT2_oo");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 2, 2, "LT2_oo");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* LT2_VV(A,B) = 0.5 * LMNEA * TMNEB + LmNeA * TmNeB */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT2_VV");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT2_VV");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* LT2_vv(a,b) = 0.5 * Lmnea * Tmneb + LMnEa * TMnEb */

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT2_vv");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 3, 3, "LT2_vv");
    dpd_buf4_init(&L2, CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I);
  }

  /* LR_OO = LR1_OO + LR2_OO */
  /* LR_oo = LR1_oo + LR2_oo */
  /* LR_VV = LR1_VV + LR2_VV */
  /* LR_vv = LR1_vv + LR2_vv */

  dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_file2_init(&LR1, EOM_TMP, S.irrep, 0, 0, "LR1_OO");
  dpd_file2_init(&LR2, EOM_TMP, S.irrep, 0, 0, "LR2_OO");
  dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
  dpd_file2_close(&LR2);
  dpd_file2_close(&LR1);
  dpd_file2_close(&I);

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 0, 0, "LR_oo");
    dpd_file2_init(&LR1, EOM_TMP, S.irrep, 0, 0, "LR1_oo");
    dpd_file2_init(&LR2, EOM_TMP, S.irrep, 0, 0, "LR2_oo");
    dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
    dpd_file2_close(&LR2);
    dpd_file2_close(&LR1);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 2, 2, "LR_oo");
    dpd_file2_init(&LR1, EOM_TMP, S.irrep, 2, 2, "LR1_oo");
    dpd_file2_init(&LR2, EOM_TMP, S.irrep, 2, 2, "LR2_oo");
    dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
    dpd_file2_close(&LR2);
    dpd_file2_close(&LR1);
    dpd_file2_close(&I);
  }

  dpd_file2_init(&I, EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_file2_init(&LR1, EOM_TMP, S.irrep, 1, 1, "LR1_VV");
  dpd_file2_init(&LR2, EOM_TMP, S.irrep, 1, 1, "LR2_VV");
  dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
  dpd_file2_close(&LR2);
  dpd_file2_close(&LR1);
  dpd_file2_close(&I);

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 1, 1, "LR_vv");
    dpd_file2_init(&LR1, EOM_TMP, S.irrep, 1, 1, "LR1_vv");
    dpd_file2_init(&LR2, EOM_TMP, S.irrep, 1, 1, "LR2_vv");
    dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
    dpd_file2_close(&LR2);
    dpd_file2_close(&LR1);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, S.irrep, 3, 3, "LR_vv");
    dpd_file2_init(&LR1, EOM_TMP, S.irrep, 3, 3, "LR1_vv");
    dpd_file2_init(&LR2, EOM_TMP, S.irrep, 3, 3, "LR2_vv");
    dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
    dpd_file2_close(&LR2);
    dpd_file2_close(&LR1);
    dpd_file2_close(&I);
  }

  /* LT_OO = LT1_OO + LT2_OO */
  /* LT_oo = LT1_oo + LT2_oo */

  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT_OO");
  dpd_file2_init(&LT1, EOM_TMP, 0, 0, 0, "LT1_OO");
  dpd_file2_init(&LT2, EOM_TMP, 0, 0, 0, "LT2_OO");
  dpd_file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
  dpd_file2_close(&LT2);
  dpd_file2_close(&LT1);
  dpd_file2_close(&I);

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT_oo");
    dpd_file2_init(&LT1, EOM_TMP, 0, 0, 0, "LT1_oo");
    dpd_file2_init(&LT2, EOM_TMP, 0, 0, 0, "LT2_oo");
    dpd_file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
    dpd_file2_close(&LT2);
    dpd_file2_close(&LT1);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 2, 2, "LT_oo");
    dpd_file2_init(&LT1, EOM_TMP, 0, 2, 2, "LT1_oo");
    dpd_file2_init(&LT2, EOM_TMP, 0, 2, 2, "LT2_oo");
    dpd_file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
    dpd_file2_close(&LT2);
    dpd_file2_close(&LT1);
    dpd_file2_close(&I);
  }

  /* LT_VV = LT1_VV + LT2_VV */
  /* LT_vv = LT1_vv + LT2_vv */

  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT_VV");
  dpd_file2_init(&LT1, EOM_TMP, 0, 1, 1, "LT1_VV");
  dpd_file2_init(&LT2, EOM_TMP, 0, 1, 1, "LT2_VV");
  dpd_file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
  dpd_file2_close(&LT2);
  dpd_file2_close(&LT1);
  dpd_file2_close(&I);

  if (rohf) {
    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT_vv");
    dpd_file2_init(&LT1, EOM_TMP, 0, 1, 1, "LT1_vv");
    dpd_file2_init(&LT2, EOM_TMP, 0, 1, 1, "LT2_vv");
    dpd_file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
    dpd_file2_close(&LT2);
    dpd_file2_close(&LT1);
    dpd_file2_close(&I);
  }
  else {
    dpd_file2_init(&I, EOM_TMP, 0, 3, 3, "LT_vv");
    dpd_file2_init(&LT1, EOM_TMP, 0, 3, 3, "LT1_vv");
    dpd_file2_init(&LT2, EOM_TMP, 0, 3, 3, "LT2_vv");
    dpd_file2_axpbycz(&LT1, &LT2, &I, 1.0, 1.0, 0.0);
    dpd_file2_close(&LT2);
    dpd_file2_close(&LT1);
    dpd_file2_close(&I);
  }

  return;
}

}} // namespace psi::ccdensity
