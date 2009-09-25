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

void rtdensity(struct TD_Params S)
{
  dpdfile2 D, T1, L1, Z;
  dpdbuf4 T2, L2;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&D, CC_TMP, S.irrep, 0, 0, "RTDIJ");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 0, 0, "RTDij");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 7, 2, 7, 0, "Lijab");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "Lia");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 1, 1, "RTDAB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 1, 1, "RTDab");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "RTDIA");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Z, CC_TMP0, S.irrep, 0, 0, "Z(I,M)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);

    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */

    dpd_file2_init(&Z, CC_TMP0, S.irrep, 0, 0, "Z(I,M)");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);

    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */

    dpd_file2_init(&Z, CC_TMP0, S.irrep, 1, 1, "Z(A,E)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "RTDia");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Z, CC_TMP0, S.irrep, 0, 0, "Z(i,m)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_init(&Z, CC_TMP0, S.irrep, 0, 0, "Z(i,m)");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    dpd_file2_init(&Z, CC_TMP0, S.irrep, 1, 1, "Z(a,e)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 2, 5, 2, 7, 0, "Lijab");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_TMP, "RTDAI");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "Lia");
    dpd_file2_copy(&L1, CC_TMP, "RTDai");
    dpd_file2_close(&L1);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&D, CC_TMP, S.irrep, 0, 0, "RTDIJ");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 2, 2, "RTDij");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 10, 17, 12, 17, 0, "Lijab");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&L1, CC_GL, S.irrep, 2, 3, "Lia");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 1, 1, "RTDAB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 3, 3, "RTDab");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 12, 15, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GL, S.irrep, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "RTDIA");

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GL, S.irrep, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Z, CC_TMP0, S.irrep, 0, 0, "Z(I,M)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);

    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */

    dpd_file2_init(&Z, CC_TMP0, S.irrep, 0, 0, "Z(I,M)");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);

    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */

    dpd_file2_init(&Z, CC_TMP0, S.irrep, 1, 1, "Z(A,E)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_TMP, S.irrep, 2, 3, "RTDia");

    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_file2_init(&L1, CC_GL, S.irrep, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GL, S.irrep, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Z, CC_TMP0, S.irrep, 2, 2, "Z(i,m)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_init(&Z, CC_TMP0, S.irrep, 2, 2, "Z(i,m)");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 10, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    dpd_file2_init(&Z, CC_TMP0, S.irrep, 3, 3, "Z(a,e)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_file2_init(&L1, CC_GL, S.irrep, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_TMP, "RTDAI");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GL, S.irrep, 2, 3, "Lia");
    dpd_file2_copy(&L1, CC_TMP, "RTDai");
    dpd_file2_close(&L1);
  }

  return;
}

}} // namespace psi::ccdensity
