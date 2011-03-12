/*! \file
    \ingroup ccresponse
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

void hbar_extra(void) {
  dpdfile2 lt;
  dpdbuf4 W1, W2, W;
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
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_buf4_copy(&W1, CC_HBAR, "2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_init(&W2, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_buf4_axpy(&W2, &W1, 2);
  dpd_buf4_close(&W2);
  dpd_buf4_sort(&W1, CC_HBAR, rspq, 10, 10, "2 W(jb,ME) + W(Jb,Me)");
  dpd_buf4_close(&W1);
}

}} // namespace psi::ccresponse
