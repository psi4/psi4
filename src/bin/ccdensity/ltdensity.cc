/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void ltdensity_rohf(struct TD_Params S);
void ltdensity_uhf(struct TD_Params S);
void ltdensity_intermediates(struct TD_Params S);
void sort_ltd_rohf(struct TD_Params S);
void sort_ltd_uhf(struct TD_Params S);

void ltdensity_rohf(struct TD_Params S)
{
  dpdfile2 DAI, Dai, DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, Int, XIJ, Xij, R1;
  dpdbuf4 T2, L2, R2, I2;
  dpdfile2 D, T1, L1, Z;

  if(S.irrep == 0) {
    dpd_file2_init(&D, CC_TMP, 0, 0, 0, "LTDIJ");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, 0, 0, 0, "LTDij");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, 0, 1, 1, "LTDAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, 0, 1, 1, "LTDab");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);
  }

  /* R_I^A */
  dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "RIA");
  dpd_file2_copy(&R1, CC_TMP, "LTDIA");
  dpd_file2_close(&R1);

  if(S.irrep == 0) {
    dpd_file2_init(&D, CC_TMP, 0, 0, 1, "LTDIA");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);

    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */

    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);

    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */

    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(A,E)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);
  }

  /* R_i^a */
  dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "Ria");
  dpd_file2_copy(&R1, CC_TMP, "LTDia");
  dpd_file2_close(&R1);

  if(S.irrep == 0) {
    dpd_file2_init(&D, CC_TMP, 0, 0, 1, "LTDia");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(i,m)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(i,m)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(a,e)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_TMP, "LTDAI");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_copy(&L1, CC_TMP, "LTDai");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_scm(&L1, (1/S.R0));
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_scm(&L1, (1/S.R0));
    dpd_file2_close(&L1);
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_scm(&L2, (1/S.R0));
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_scm(&L2, (1/S.R0));
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_scm(&L2, (1/S.R0));
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_scm(&L2, (1/S.R0));
    dpd_buf4_close(&L2);
  }

  ltdensity_intermediates(S);

  dpd_file2_init(&TIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&Tia, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_init(&RIA, CC_GR, S.irrep, 0, 1, "RIA");
  dpd_file2_init(&Ria, CC_GR, S.irrep, 0, 1, "Ria");
  dpd_file2_init(&LIA, CC_GLG, 0, 0, 1, "LIA");
  dpd_file2_init(&Lia, CC_GLG, 0, 0, 1, "Lia");

  /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */

  dpd_file2_init(&DIJ, CC_TMP, S.irrep, 0, 0, "LTDIJ");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_file2_axpy(&Int, &DIJ, -1.0, 1);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_contract222(&TIA, &Int, &DIJ, 0, 0, -1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_file2_close(&DIJ);

  dpd_file2_init(&Dij, CC_TMP, S.irrep, 0, 0, "LTDij");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 0, "LR_oo");
  dpd_file2_axpy(&Int, &Dij, -1.0, 1);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_contract222(&Tia, &Int, &Dij, 0, 0, -1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_file2_close(&Dij);

  /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */

  dpd_file2_init(&DAB, CC_TMP, S.irrep, 1, 1, "LTDAB");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_file2_axpy(&Int, &DAB, 1.0, 0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_contract222(&Int, &TIA, &DAB, 1, 1, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_file2_close(&DAB);

  dpd_file2_init(&Dab, CC_TMP, S.irrep, 1, 1, "LTDab");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 1, 1, "LR_vv");
  dpd_file2_axpy(&Int, &Dab, 1.0, 0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_contract222(&Int, &Tia, &Dab, 1, 1, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_file2_close(&Dab);

  /* D[a][i] = +L2R1_ov[i][a] */

  dpd_file2_init(&DAI, CC_TMP, S.irrep, 0, 1, "LTDAI");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_file2_axpy(&Int, &DAI, 1.0, 0);
  dpd_file2_close(&Int);
  dpd_file2_close(&DAI);

  dpd_file2_init(&Dai, CC_TMP, S.irrep, 0, 1, "LTDai");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_file2_axpy(&Int, &Dai, 1.0, 0);
  dpd_file2_close(&Int);
  dpd_file2_close(&Dai);

  dpd_file2_init(&DIA, CC_TMP, S.irrep, 0, 1, "LTDIA");
  dpd_file2_init(&Dia, CC_TMP, S.irrep, 0, 1, "LTDia");

  /* D[i][a] = L1R2_ov[i][a] */

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L1R2_OV");
  dpd_file2_axpy(&Int, &DIA, 1.0, 0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L1R2_ov");
  dpd_file2_axpy(&Int, &Dia, 1.0, 0);
  dpd_file2_close(&Int);

  /* - LR_OO[M][I] * t1[M][A] */

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_contract222(&Int, &TIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 0, "LR_oo");
  dpd_contract222(&Int, &Tia, &Dia, 1, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  /* - t1[I][E] * LR_vv[E][A] */

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_contract222(&TIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 1, 1, "LR_vv");
  dpd_contract222(&Tia, &Int, &Dia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  /* - LT2_OO[M][I] * r1[M][A] */

  dpd_file2_init(&Int, EOM_TMP, 0, 0, 0, "LT2_OO");
  dpd_contract222(&Int, &RIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, 0, 0, 0, "LT2_oo");
  dpd_contract222(&Int, &Ria, &Dia, 1, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  /* - r1[I][E] * LT2_VV[E][A] */

  dpd_file2_init(&Int, EOM_TMP, 0, 1, 1, "LT2_VV");
  dpd_contract222(&RIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, 0, 1, 1, "LT2_vv");
  dpd_contract222(&Ria, &Int, &Dia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  /* + L2R1_ov[M][E] * t2[i][m][a][e] */

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB"); 
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb"); 
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab"); 
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB"); 
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_buf4_close(&T2);
    
  /* - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_file2_init(&XIJ, EOM_TMP, S.irrep, 0, 0, "XIJ");
  dpd_contract222(&TIA, &Int, &XIJ, 0, 0, 1.0, 0.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&XIJ, EOM_TMP, S.irrep, 0, 0, "XIJ");
  dpd_contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&XIJ);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_file2_init(&Xij, EOM_TMP, S.irrep, 0, 0, "Xij");
  dpd_contract222(&Tia, &Int, &Xij, 0, 0, 1.0, 0.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Xij, EOM_TMP, S.irrep, 0, 0, "Xij");
  dpd_contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Xij);

  dpd_file2_close(&DIA);
  dpd_file2_close(&Dia);

  dpd_file2_close(&TIA);
  dpd_file2_close(&Tia);
  dpd_file2_close(&RIA);
  dpd_file2_close(&Ria);
  dpd_file2_close(&LIA);
  dpd_file2_close(&Lia);

  return;
}

void ltdensity_uhf(struct TD_Params S)
{
  dpdfile2 DAI, Dai, DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, Int, XIJ, Xij, R1;
  dpdbuf4 T2, L2, R2, I2;
  dpdfile2 D, T1, L1, Z;

  if(S.irrep == 0 ) { /* Symmetric Transition */

    dpd_file2_init(&D, CC_TMP, 0, 0, 0, "LTDIJ");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, 0, 2, 2, "LTDij");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, 0, 1, 1, "LTDAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, 0, 3, 3, "LTDab");
    dpd_buf4_init(&L2, CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);
  }

  /* R_I^A */
    
  dpd_file2_init(&R1, CC_GR, S.irrep, 0, 1, "RIA");
  dpd_file2_copy(&R1, CC_TMP, "LTDIA");
  dpd_file2_close(&R1); 

  if(S.irrep == 0) { // Symmetric Transitions

    dpd_file2_init(&D, CC_TMP, 0, 0, 1, "LTDIA");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);

    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */

    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(A,E)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);
  }

  /* R_i^a */
    
  dpd_file2_init(&R1, CC_GR, S.irrep, 2, 3, "Ria");
  dpd_file2_copy(&R1, CC_TMP, "LTDia"); 
  dpd_file2_close(&R1);

  if(S.irrep == 0) { /* Symmetric Transitions */

    dpd_file2_init(&D, CC_TMP, 0, 2, 3, "LTDia");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Z, CC_TMP0, 0, 2, 2, "Z(i,m)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_init(&Z, CC_TMP0, 0, 2, 2, "Z(i,m)");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    dpd_file2_init(&Z, CC_TMP0, 0, 3, 3, "Z(a,e)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_TMP, "LTDAI");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_copy(&L1, CC_TMP, "LTDai");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_scm(&L1, (1/S.R0));
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_scm(&L1, (1/S.R0));
    dpd_file2_close(&L1);
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_scm(&L2, (1/S.R0));
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_buf4_scm(&L2, (1/S.R0));
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_scm(&L2, (1/S.R0));
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_scm(&L2, (1/S.R0));
    dpd_buf4_close(&L2);
  }

  ltdensity_intermediates(S);

  dpd_file2_init(&TIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&Tia, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_init(&RIA, CC_GR, S.irrep, 0, 1, "RIA");
  dpd_file2_init(&Ria, CC_GR, S.irrep, 2, 3, "Ria");
  dpd_file2_init(&LIA, CC_GL, S.irrep, 0, 1, "LIA");
  dpd_file2_init(&Lia, CC_GL, S.irrep, 2, 3, "Lia");

  /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */

  dpd_file2_init(&DIJ, CC_TMP, S.irrep, 0, 0, "LTDIJ");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_file2_axpy(&Int, &DIJ, -1.0, 1);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_contract222(&TIA, &Int, &DIJ, 0, 0, -1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_file2_close(&DIJ);

  dpd_file2_init(&Dij, CC_TMP, S.irrep, 2, 2, "LTDij");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 2, "LR_oo");
  dpd_file2_axpy(&Int, &Dij, -1.0, 1);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_contract222(&Tia, &Int, &Dij, 0, 0, -1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_file2_close(&Dij);

  /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */

  dpd_file2_init(&DAB, CC_TMP, S.irrep, 1, 1, "LTDAB");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_file2_axpy(&Int, &DAB, 1.0, 0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_contract222(&Int, &TIA, &DAB, 1, 1, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_file2_close(&DAB);

  dpd_file2_init(&Dab, CC_TMP, S.irrep, 3, 3, "LTDab");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 3, 3, "LR_vv");
  dpd_file2_axpy(&Int, &Dab, 1.0, 0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_contract222(&Int, &Tia, &Dab, 1, 1, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_file2_close(&Dab);

  /* D[a][i] = +L2R1_ov[i][a] */

  dpd_file2_init(&DAI, CC_TMP, S.irrep, 0, 1, "LTDAI");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_file2_axpy(&Int, &DAI, 1.0, 0);
  dpd_file2_close(&Int);
  dpd_file2_close(&DAI);

  dpd_file2_init(&Dai, CC_TMP, S.irrep, 2, 3, "LTDai");
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_file2_axpy(&Int, &Dai, 1.0, 0);
  dpd_file2_close(&Int);
  dpd_file2_close(&Dai);

  dpd_file2_init(&DIA, CC_TMP, S.irrep, 0, 1, "LTDIA");

  dpd_file2_init(&Dia, CC_TMP, S.irrep, 2, 3, "LTDia");

  /* D[i][a] = L1R2_ov[i][a] */

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L1R2_OV");
  dpd_file2_axpy(&Int, &DIA, 1.0, 0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 3, "L1R2_ov");
  dpd_file2_axpy(&Int, &Dia, 1.0, 0);
  dpd_file2_close(&Int);

  /* - LR_OO[M][I] * t1[M][A] */

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_contract222(&Int, &TIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 2, "LR_oo");
  dpd_contract222(&Int, &Tia, &Dia, 1, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  /* - t1[I][E] * LR_vv[E][A] */

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_contract222(&TIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 3, 3, "LR_vv");
  dpd_contract222(&Tia, &Int, &Dia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  /* - LT2_OO[M][I] * r1[M][A] */

  dpd_file2_init(&Int, EOM_TMP, 0, 0, 0, "LT2_OO");
  dpd_contract222(&Int, &RIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, 0, 2, 2, "LT2_oo");
  dpd_contract222(&Int, &Ria, &Dia, 1, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  /* - r1[I][E] * LT2_vv[E][A] */

  dpd_file2_init(&Int, EOM_TMP, 0, 1, 1, "LT2_VV");
  dpd_contract222(&RIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Int, EOM_TMP, 0, 3, 3, "LT2_vv");
  dpd_contract222(&Ria, &Int, &Dia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Int);

  /* + L2R1_ov[M][E] * t2[i][m][a][e] */

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB"); 
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb"); 
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab"); 
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB"); 
  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  dpd_file2_close(&Int);
  dpd_buf4_close(&T2);
    
  /* - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_file2_init(&XIJ, EOM_TMP, S.irrep, 0, 0, "XIJ");
  dpd_contract222(&TIA, &Int, &XIJ, 0, 0, 1.0, 0.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&XIJ, EOM_TMP, S.irrep, 0, 0, "XIJ");
  dpd_contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&XIJ);

  dpd_file2_init(&Int, EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_file2_init(&Xij, EOM_TMP, S.irrep, 2, 2, "Xij");
  dpd_contract222(&Tia, &Int, &Xij, 0, 0, 1.0, 0.0);
  dpd_file2_close(&Int);

  dpd_file2_init(&Xij, EOM_TMP, S.irrep, 2, 2, "Xij");
  dpd_contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Xij);

  dpd_file2_close(&DIA);
  dpd_file2_close(&Dia);

  dpd_file2_close(&TIA);
  dpd_file2_close(&Tia);
  dpd_file2_close(&RIA);
  dpd_file2_close(&Ria);
  dpd_file2_close(&LIA);
  dpd_file2_close(&Lia);

  return;
}

}} // namespace psi::ccdensity
