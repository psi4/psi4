/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void cc2_hbar_extra(void) {
  dpdfile2 t1, lt;
  dpdbuf4 A, D, E, Z, Z1;
  dpdbuf4 W1, W;
  dpdbuf4 t2, l2;

  /* LIjAb * TIjAb */
  dpd_file2_init(&lt, CC_OEI, 0, 0, 0, "Lt_IJ");

  dpd_buf4_init(&l2, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB 0 -1");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_contract442(&l2, &t2, &lt, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&t2);
  dpd_buf4_close(&l2);

  dpd_buf4_init(&l2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb 0 -1");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract442(&l2, &t2, &lt, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&t2);
  dpd_buf4_close(&l2);

  dpd_file2_close(&lt);

  /* 2 W(ME,jb) + W(Me,Jb) */
  dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_sort(&W1, CC2_HET1, rspq, 10, 10, "CC2 2 W(jb,ME) + W(Jb,Me)");
  dpd_buf4_close(&W1);

  /* CC2 WMnIj */

/*   dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>"); */
/*   dpd_buf4_copy(&A, CC2_HET1, "CC2 WMnIj (Mn,Ij)"); */
/*   dpd_buf4_close(&A); */

/*   dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA"); */

/*   /\* Wmnij <- + P(ij) t(j,e) * <mn||ie> *\/ */
/*   dpd_buf4_init(&Z, CC_TMP0, 0, 0, 0, 0, 0, 0, "CC2 ZMnIj (Mn,Ij)"); */
/*   dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>"); */
/*   dpd_contract424(&E, &t1, &Z, 3, 1, 0, 1, 0); */
/*   dpd_buf4_close(&E); */

/*   dpd_buf4_init(&W, CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj (Mn,Ij)"); */
/*   dpd_buf4_axpy(&Z, &W, 1); */
/*   dpd_buf4_close(&W); */
/*   dpd_buf4_sort_axpy(&Z, CC2_HET1, qpsr, 0, 0, "CC2 WMnIj (Mn,Ij)", 1); */
/*   dpd_buf4_close(&Z); */

/*   /\* Wmnij<- +1/2 P(ij) t(i,e) t(j,f) * <mn||ef> *\/ */
/*   dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "CC2 ZMnIf (Mn,If)"); */
/*   dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>"); */
/*   dpd_contract244(&t1, &D, &Z, 1, 2, 1, 1, 0); */
/*   dpd_buf4_close(&D); */

/*   dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 0, 0, 0, 0, "CC2 ZMnIj (Mn,Ij)"); */
/*   dpd_contract424(&Z, &t1, &Z1, 3, 1, 0, 0.5, 0); */
/*   dpd_buf4_close(&Z); */
/*   dpd_buf4_init(&W, CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj (Mn,Ij)"); */
/*   dpd_buf4_axpy(&Z1, &W, 1); */
/*   dpd_buf4_close(&W); */
/*   dpd_buf4_sort_axpy(&Z1, CC2_HET1, qpsr, 0, 0, "CC2 WMnIj (Mn,Ij)", 1); */
/*   dpd_buf4_close(&Z1); */

/*   dpd_file2_close(&t1); */
}

}} // namespace psi::ccresponse
