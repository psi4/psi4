/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void cc2_hbar_extra(void) {
  dpdbuf4 W1, W2, W;

  if(params.ref == 0) { /** RHF **/
    
    dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Me,Jb)");
    dpd_buf4_copy(&W1, CC2_HET1, "CC2 2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_close(&W1);
    dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_init(&W2, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj (ME,jb)");
    dpd_buf4_axpy(&W2, &W1, 2);
    dpd_buf4_close(&W2);
    dpd_buf4_close(&W1);
    dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_sort(&W1, CC2_HET1, rspq, 10, 10, "CC2 2 W(jb,ME) + W(Jb,Me)");
    dpd_buf4_close(&W1);

    /*
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_buf4_scmcopy(&W, CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    dpd_buf4_sort_axpy(&W, CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
    dpd_buf4_sort(&W, CC2_HET1, rspq, 11, 5, "CC2 WAbEi (Ei,Ab)");
    dpd_buf4_close(&W);
    */
  }
}

}} // namespace psi::cceom
