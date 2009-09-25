/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* compute remaining matrix elements for [ H, e^T1 ] which are not needed by
 * ccenergy for CC3 but are needed for CC3 EOM and CC3 Lambda and CC3 response.
 *
 * Wmbej = <mb||ej> + t_j^f <mb||ef> - t_n^b <mn||ej> - t_j^f t_n^b <mn||ef> 
 *
*/

void HET1_Wmbej(void);
void purge_HET1_Wmbej(void);
void HET1_Wabef(void);

void cc3_HET1(void) {

  HET1_Wmbej();

  if (params.ref == 1) purge_HET1_Wmbej();

  HET1_Wabef();

  return;
}


void HET1_Wmbej(void)
{
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ;
  dpdbuf4 D, C, F, E, X, Y, t2, W, Z;
  dpdfile2 tIA, tia;
  int i;

  if(params.ref == 0) { /** RHF **/

    /* <mb||ej> -> Wmbej */

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMbeJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "WMbEj");
    dpd_buf4_close(&D);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    /* F -> Wmbej */

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");

    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&F, &tIA, &WMbEj, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WMbEj);

    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");

    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_sort(&Z, CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_axpy(&Z, &WMbeJ, 1.0);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_close(&Z);

    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);

    /* E -> Wmbej */
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&E, &tIA, &WMbEj, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract424(&E, &tIA, &WMbeJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA);  


    /* Sort to (ME,JB) */

    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_buf4_sort(&WMbEj, CC3_HET1, prsq, 10, 10, "CC3 WMbEj (ME,jb)");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_sort(&WMbeJ, CC3_HET1, psrq, 10, 10, "CC3 WMbeJ (Me,Jb)");
    dpd_buf4_close(&WMbeJ);


    /* T1^2 -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);

    /* also store WMbEj (Mb,Ej)  WMbeJ (bM,eJ)*/
    dpd_buf4_init(&WMbEj, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
    dpd_buf4_sort(&WMbEj, CC3_HET1, psqr, 10, 11, "CC3 WMbEj (Mb,Ej)");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WMbEj, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
    dpd_buf4_sort(&WMbEj, CC3_HET1, spqr, 11, 11, "CC3 WMbeJ (bM,eJ)");
    dpd_buf4_close(&WMbEj);

  }
  else if(params.ref == 1) { /** ROHF **/

    /* W(mb,je) <-- <mb||ej> */

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMBEJ", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "Wmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WmBEj", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMbeJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "WMbEj");
    dpd_buf4_copy(&D, CC_TMP0, "WmBeJ");
    dpd_buf4_close(&D);

    /* F -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_contract424(&F, &tIA, &WMBEJ, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WMBEJ);
    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    dpd_contract424(&F, &tia, &Wmbej, 3, 1, 0, 1, 1);
    dpd_buf4_close(&Wmbej);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&F, &tia, &WMbEj, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_contract424(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WmBeJ);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract244(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_contract244(&tia, &F, &WmBEj, 1, 2, 1, -1, 1);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* E -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_contract424(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMBEJ);
    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    dpd_contract424(&E, &tia, &Wmbej, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Wmbej);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&E, &tia, &WMbEj, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_contract424(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WmBeJ);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract424(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_contract424(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);

    /* Convert to (ME,JB) for remaining terms */

    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_buf4_sort(&WMBEJ, CC3_HET1, prsq, 10, 10, "CC3 WMBEJ (ME,JB)");
    dpd_buf4_close(&WMBEJ);

    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    dpd_buf4_sort(&Wmbej, CC3_HET1, prsq, 10, 10, "CC3 Wmbej (me,jb)");
    dpd_buf4_close(&Wmbej);

    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    dpd_buf4_sort(&WMbEj, CC3_HET1, prsq, 10, 10, "CC3 WMbEj (ME,jb)");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_buf4_sort(&WmBeJ, CC3_HET1, prsq, 10, 10, "CC3 WmBeJ (me,JB)");
    dpd_buf4_close(&WmBeJ);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_sort(&WMbeJ, CC3_HET1, psrq, 10, 10, "CC3 WMbeJ (Me,Jb)");
    dpd_buf4_close(&WMbeJ);

    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_buf4_sort(&WmBEj, CC3_HET1, psrq, 10, 10, "CC3 WmBEj (mE,jB)");
    dpd_buf4_close(&WmBEj);

    /* X -> Wmbej */
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /*** AAAA ***/

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMBEJ (ME,JB)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_sort(&W, CC3_HET1, psqr, 10, 11, "CC3 WMBEJ (MB,EJ)");
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BBBB ***/

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 Wmbej (me,jb)");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_sort(&W, CC3_HET1, psqr, 10, 11, "CC3 Wmbej (mb,ej)");
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABAB ***/
  
    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_sort(&W, CC3_HET1, psqr, 10, 11, "CC3 WMbEj (Mb,Ej)");
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BABA ***/

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WmBeJ (me,JB)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_sort(&W, CC3_HET1, psqr, 10, 11, "CC3 WmBeJ (mB,eJ)");
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABBA ***/
  
    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    dpd_buf4_sort(&W, CC3_HET1, psqr, 10, 11, "CC3 WMbeJ (Mb,eJ)");
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BAAB ***/

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WmBEj (mE,jB)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_sort(&W, CC3_HET1, psqr, 10, 11, "CC3 WmBEj (mB,Ej)");
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* W(mb,je) <-- <mb||ej> */

    dpd_buf4_init(&C, CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMBEJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "Wmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_buf4_scmcopy(&D, CC_TMP0, "WMbEj", 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_buf4_scmcopy(&D, CC_TMP0, "WmBeJ", 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&C, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WmBEj", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMbeJ", -1);
    dpd_buf4_close(&C);

    /* F -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract244(&tIA, &F, &W, 1, 2, 1, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract244(&tia, &F, &W, 1, 2, 1, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* E -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    dpd_contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_init(&E, CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    dpd_contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract424(&E, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    dpd_contract424(&E, &tIA, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);

    /* Convert to (ME,JB) for remaining terms */

    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_sort(&W, CC3_HET1, prsq, 20, 20, "CC3 WMBEJ (ME,JB)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_sort(&W, CC3_HET1, prsq, 30, 30, "CC3 Wmbej (me,jb)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_sort(&W, CC3_HET1, prsq, 20, 30, "CC3 WMbEj (ME,jb)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_sort(&W, CC3_HET1, prsq, 30, 20, "CC3 WmBeJ (me,JB)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_sort(&W, CC3_HET1, psrq, 24, 24, "CC3 WMbeJ (Me,Jb)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_sort(&W, CC3_HET1, psrq, 27, 27, "CC3 WmBEj (mE,jB)");
    dpd_buf4_close(&W);

    /* X -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /*** AAAA ***/

    dpd_buf4_init(&Y, CC_TMP0, 0, 20, 0, 20, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 21, 20, 21, 0, "D <IJ||AB> (IA,BJ)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 20, 20, 20, 0, "CC3 WMBEJ (ME,JB)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BBBB ***/

    dpd_buf4_init(&Y, CC_TMP0, 0, 30, 10, 30, 10, 0, "Y (me,jn)");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 31, 30, 31, 0, "D <ij||ab> (ia,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 30, 30, 30, 0, "CC3 Wmbej (me,jb)");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABAB ***/

    dpd_buf4_init(&Y, CC_TMP0, 0, 20, 10, 20, 10, 0, "Y (ME,jn)");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 31, 20, 31, 0, "D <Ij|Ab> (IA,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 30, 20, 30, 0, "CC3 WMbEj (ME,jb)");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BABA ***/

    dpd_buf4_init(&Y, CC_TMP0, 0, 30, 0, 30, 0, 0, "Y (me,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 21, 30, 21, 0, "D <Ij|Ab> (ia,BJ)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 20, 30, 20, 0, "CC3 WmBeJ (me,JB)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABBA ***/
  
    dpd_buf4_init(&Y, CC_TMP0, 0, 24, 22, 24, 22, 0, "Y (Me,Jn)");
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 24, 24, 24, 24, 0, "CC3 WMbeJ (Me,Jb)");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BAAB ***/
  
    dpd_buf4_init(&Y, CC_TMP0, 0, 27, 23, 27, 23, 0, "Y (mE,jN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 27, 27, 27, 27, 0, "CC3 WmBEj (mE,jB)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* also store lists as Wmbej (mb,ej) */
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 20, 20, 20, 0, "CC3 WMBEJ (ME,JB)");
    dpd_buf4_sort(&W, CC3_HET1, psqr, 20, 21, "CC3 WMBEJ (MB,EJ)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 30, 30, 30, 0, "CC3 Wmbej (me,jb)");
    dpd_buf4_sort(&W, CC3_HET1, psqr, 30, 31, "CC3 Wmbej (mb,ej)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 30, 20, 30, 0, "CC3 WMbEj (ME,jb)");
    dpd_buf4_sort(&W, CC3_HET1, psqr, 24, 26, "CC3 WMbEj (Mb,Ej)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 20, 30, 20, 0, "CC3 WmBeJ (me,JB)");
    dpd_buf4_sort(&W, CC3_HET1, psqr, 27, 25, "CC3 WmBeJ (mB,eJ)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 24, 24, 24, 24, 0, "CC3 WMbeJ (Me,Jb)");
    dpd_buf4_sort(&W, CC3_HET1, psqr, 24, 25, "CC3 WMbeJ (Mb,eJ)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 27, 27, 27, 27, 0, "CC3 WmBEj (mE,jB)");
    dpd_buf4_sort(&W, CC3_HET1, psqr, 27, 26, "CC3 WmBEj (mB,Ej)");
    dpd_buf4_close(&W);

    /* also make some Wmbej (bm,ej) */
    dpd_buf4_init(&W, CC3_HET1, 0, 24, 25, 24, 25, 0, "CC3 WMbeJ (Mb,eJ)");
    dpd_buf4_sort(&W, CC3_HET1, qprs, 25, 25, "CC3 WMbeJ (bM,eJ)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC3_HET1, 0, 27, 26, 27, 26, 0, "CC3 WmBEj (mB,Ej)");
    dpd_buf4_sort(&W, CC3_HET1, qprs, 26, 26, "CC3 WmBEj (Bm,Ej)");
    dpd_buf4_close(&W);

  } /** UHF **/
  return;
}



/* Purge Wmbej matrix elements for ROHF */
void purge_HET1_Wmbej(void) {
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

  dpd_file4_init(&W, CC3_HET1, 0, 10, 10,"CC3 WMBEJ (ME,JB)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      e = W.params->roworb[h][me][1];
      esym = W.params->qsym[e];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	b = W.params->colorb[h][jb][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  dpd_file4_init(&W, CC3_HET1, 0, 10, 10,"CC3 Wmbej (me,jb)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	jsym = W.params->rsym[j];
	J = j - occ_off[jsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  dpd_file4_init(&W, CC3_HET1, 0, 10, 10,"CC3 WMbEj (ME,jb)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      e = W.params->roworb[h][me][1];
      esym = W.params->qsym[e];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	jsym = W.params->rsym[j];
	J = j - occ_off[jsym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  dpd_file4_init(&W, CC3_HET1, 0, 10, 10,"CC3 WmBeJ (me,JB)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	b = W.params->colorb[h][jb][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  dpd_file4_init(&W, CC3_HET1, 0, 10, 10,"CC3 WmBEj (mE,jB)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      e = W.params->roworb[h][me][1];
      msym = W.params->psym[m];
      esym = W.params->qsym[e];
      M = m - occ_off[msym];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	b = W.params->colorb[h][jb][1];
	jsym = W.params->rsym[j];
	bsym = W.params->ssym[b];
	J = j - occ_off[jsym];
	B = b - vir_off[bsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (E >= (virtpi[esym] - openpi[esym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);
  return;
}

}} // namespace psi::cchbar
