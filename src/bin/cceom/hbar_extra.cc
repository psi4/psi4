/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void hbar_extra(void) {
  dpdbuf4 W, W1, W2, WAmEf, WmBeJ, WmBEj, WmNIe, WMnIe;

  if (params.eom_ref == 2) {
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 20, 20, 20, 20, 0, "WMBEJ");
    dpd_buf4_sort(&W, CC_HBAR, rspq, 20, 20, "WMBEJ (JB,ME)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 30, 20, 30, 20, 0, "WmBeJ"); /* (me,JB) */
    dpd_buf4_sort(&W, CC_HBAR, rspq, 20, 30, "WmBeJ (JB,me)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 30, 30, 30, 30, 0, "Wmbej");
    dpd_buf4_sort(&W, CC_HBAR, rspq, 30, 30, "Wmbej (jb,me)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 20, 30, 20, 30, 0, "WMbEj"); /* (ME,jb) */
    dpd_buf4_sort(&W, CC_HBAR, rspq, 30, 20, "WMbEj (jb,ME)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_buf4_sort(&W, CC_HBAR, pqsr, 27, 22, "WmBiJ (mB,Ji)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 27, 22, 27, 22, 0, "WmBiJ (mB,Ji)");
    dpd_buf4_sort(&W, CC_HBAR, qprs, 26, 22, "WmBiJ (Bm,Ji)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 25, 29, 25, 29, 0, "WeIaB");
    dpd_buf4_sort(&W, CC_HBAR, qprs, 24, 29, "WeIaB (Ie,aB)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 24, 29, 24, 29, 0, "WeIaB (Ie,aB)");
    dpd_buf4_sort(&W, CC_HBAR, pqsr, 24, 28, "WeIaB (Ie,Ab)");
    dpd_buf4_close(&W);
  }

  if(params.eom_ref == 1) {

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMBEJ");
    dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "WMBEJ (JB,ME)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WmBeJ");
    dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "WmBeJ (JB,me)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "Wmbej");
    dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "Wmbej (jb,me)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbEj");
    dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "WMbEj (jb,ME)");
    dpd_buf4_close(&W);

  }

  if (params.eom_ref == 1) {  /* ROHF */

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_buf4_sort(&W, CC_HBAR, pqsr, 10, 0, "WmBiJ (mB,Ji)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WmBiJ (mB,Ji)");
    dpd_buf4_sort(&W, CC_HBAR, qprs, 11, 0, "WmBiJ (Bm,Ji)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WeIaB");
    dpd_buf4_sort(&W, CC_HBAR, qprs, 10, 5, "WeIaB (Ie,aB)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WeIaB (Ie,aB)");
    dpd_buf4_sort(&W, CC_HBAR, pqsr, 10, 5, "WeIaB (Ie,Ab)");
    dpd_buf4_close(&W);
  }

  if (params.eom_ref == 0 ) { /* RHF */
    /* 2 W(ME,jb) + W(Me,Jb) */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_copy(&W, CC_HBAR, "2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W1, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_init(&W2, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbEj");
    dpd_buf4_axpy(&W2, &W1, 2);
    dpd_buf4_close(&W2);
    dpd_buf4_sort(&W1, CC_HBAR, rspq, 10, 10, "2 W(jb,ME) + W(Jb,Me)");
    dpd_buf4_close(&W1);

    /* used in WamefSD */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf");
    dpd_buf4_scmcopy(&W, CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    dpd_buf4_sort_axpy(&W, CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    dpd_buf4_close(&W);
  }

  return;
}

}} // namespace psi::cceom
