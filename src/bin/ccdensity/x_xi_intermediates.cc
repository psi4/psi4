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

/* build one-electron intermediates for the construction of the xi amplitudes */

void x_xi_intermediates(void)
{
  dpdfile2 L1, R1, T1, I, LR1, LR2, LT1, LT2;
  dpdbuf4 L2, T2, R2, D, H2, I2;
  int L_irr, R_irr, G_irr;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

  if (params.ref == 0) {

    /* RD_OO(I,J) = (2RImEf - RImFe) <jm|ef> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);

    /* RD_VV(E,A) = (2RMnFe - RMnEf) <mn|fa> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);

    /* R2Wamef_OV(n,g) = Rnmef Wgmef */
    dpd_file2_init(&I, EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);

    /* R2Wamef_OOOV = Rmnfg Weifg */
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(Mn,eI)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(Mn,eI)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 0, 10, "R2Wamef_OoOv");
    dpd_buf4_close(&I2);

    /* R1Wamef_VV (a,e) = Rmf Wamef */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I);

    /* R1Wmnie_OO (M,I) = Rne Wmnie */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I);

  }
  else if (params.ref == 1) {
    /* RD_OO(I,J)  =  RIMEF * <JM||EF> + RImEf <Jm|Ef> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  
    /* RD_oo(i,j)  =  Rimef * <jm||ef> + RiMeF <jM|eF> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);

    /* RD_VV(E,A)  =  RMNFE * <MN||FA> + RmNfE <mN|fA> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
    dpd_buf4_init(&R2, CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  
    /* RD_vv(e,a)  =  Rmnfe * <mn||fa> + RMnFe <Mn|Fa> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "RD_vv");
    dpd_buf4_init(&R2, CC_GR, R_irr, 2, 5, 2, 7, 0, "Rijab");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  
    /* R2Wamef_OOVO = Rmnfg Weifg -> OOOV (mn,ie) */
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 2, 11, 2, 11, 0, "Z(M>N,EI)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
    dpd_contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 2, 11, 2, 11, 0, "Z(M>N,EI)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 2, 10, "R2Wamef_OOOV"); 
    dpd_buf4_close(&I2);

    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(Mn,eI)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    dpd_contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(Mn,eI)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 0, 10, "R2Wamef_OoOv");
    dpd_buf4_close(&I2);
  
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(mN,Ei)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJAb");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(mN,Ei)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 0, 10, "R2Wamef_oOoV");
    dpd_buf4_close(&I2);

    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 2, 11, 2, 11, 0, "Z(m>n,ei)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
    dpd_contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 2, 11, 2, 11, 0, "Z(m>n,ei)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 2, 10, "R2Wamef_ooov");
    dpd_buf4_close(&I2);
  
    /* R2Wamef_OV(n,g) = Rnmef Wgmef */
    dpd_file2_init(&I, EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  
    dpd_file2_init(&I, EOM_TMP, R_irr, 0, 1, "R2Wamef_ov");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);

    /* R1Wmnie_OO (M,I) = Rne Wmnie */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 2, 10, 0, "WMNIE");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_oo");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 2, 10, 0, "Wmnie");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I);

    /* R1Wamef_VV (a,e) = Rmf Wamef */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_vv");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
  
    dpd_file2_close(&I);
  }
  else { /* UHF */
    /* RD_OO(I,J)  =  RIMEF * <JM||EF> + RImEf <Jm|Ef> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  
    /* RD_oo(i,j)  =  Rimef * <jm||ef> + RiMeF <jM|eF> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 2, 2, "RD_oo");
    dpd_buf4_init(&R2, CC_GR, R_irr, 10, 17, 12, 17, 0, "Rijab");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);

    /* RD_VV(E,A)  =  RMNFE * <MN||FA> + RmNfE <mN|fA> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
    dpd_buf4_init(&R2, CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  
    /* RD_vv(e,a)  =  Rmnfe * <mn||fa> + RMnFe <Mn|Fa> */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 3, 3, "RD_vv");
    dpd_buf4_init(&R2, CC_GR, R_irr, 12, 15, 12, 17, 0, "Rijab");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&R2);
    dpd_file2_close(&I);
  
    /* R2Wamef_OOVO = Rmnfg W(eifg) -> Z(mn,ei) -> OOOV (mn,ie) */
    /* make using (AM,EF) so we can avoid making the (MA,EF) at all*/
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 2, 21, 2, 21, 0, "Z(M>N,EI)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    dpd_contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 2, 21, 2, 21, 0, "Z(M>N,EI)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 2, 20, "R2Wamef_OOOV"); 
    dpd_buf4_close(&I2);
  
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 22, 25, 22, 25, 0, "Z(Mn,eI)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 22, 29, 22, 29, 0, "RIjaB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 22, 25, 22, 25, 0, "Z(Mn,eI)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 22, 24, "R2Wamef_OoOv");
    dpd_buf4_close(&I2);
  
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 23, 26, 23, 26, 0, "Z(mN,Ei)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 23, 28, 23, 28, 0, "RiJAb");
    dpd_buf4_init(&H2, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 23, 26, 23, 26, 0, "Z(mN,Ei)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 23, 27, "R2Wamef_oOoV");
    dpd_buf4_close(&I2);
  
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 12, 31, 12, 31, 0, "Z(m>n,ei)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 12, 17, 12, 17, 0, "Rijab");
    dpd_buf4_init(&H2, CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    dpd_contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP1, R_irr, 12, 31, 12, 31, 0, "Z(m>n,ei)");
    dpd_buf4_sort(&I2, EOM_TMP, pqsr, 12, 30, "R2Wamef_ooov");
    dpd_buf4_close(&I2);

    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);
  
    /* R2Wamef_OV(n,g) = Rnmef Wgmef */
    dpd_file2_init(&I, EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
  
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_buf4_init(&H2, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
  
    dpd_file2_close(&I);
  
    dpd_file2_init(&I, EOM_TMP, R_irr, 2, 3, "R2Wamef_ov");
  
    dpd_buf4_init(&R2, CC_GR, R_irr, 10, 17, 12, 17, 0, "Rijab");
    dpd_buf4_init(&H2, CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    dpd_buf4_init(&H2, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_close(&R2);

    dpd_file2_close(&I);

    /* R1Wmnie_OO (M,I) = Rne Wmnie */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");

    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 20, 2, 20, 0, "WMNIE");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 22, 24, 22, 24, 0, "WMnIe");
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);

    dpd_file2_close(&I);
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 2, 2, "R1Wmnie_oo");

    dpd_buf4_init(&H2, CC_HBAR, 0, 10, 30, 12, 30, 0, "Wmnie");
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 23, 27, 23, 27, 0, "WmNiE");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);

    dpd_file2_close(&I);

    /* R1Wamef_VV (a,e) = Rmf Wamef */
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");

    dpd_buf4_init(&H2, CC_HBAR, 0, 21, 5, 21, 7, 0, "WAMEF");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);

    dpd_file2_close(&I);
    dpd_file2_init(&I, EOM_TMP_XI, R_irr, 3, 3, "R1Wamef_vv");

    dpd_buf4_init(&H2, CC_HBAR, 0, 31, 15, 31, 17, 0, "Wamef");
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
  
    dpd_file2_close(&I);
  }
  return;
}

}} // namespace psi::ccdensity
