/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void cc2_WabijT2(void)
{
  dpdbuf4 W;

  if(params.ref == 0) { /*** RHF ***/

    dpd_buf4_init(&W, CC2_HET1, 0, 0, 5, 0, 5, 0, "CC2 WAbIj (Ij,Ab)");
    dpd_buf4_copy(&W, CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&W);

  }

  else if(params.ref == 1) { /*** ROHF ***/

    dpd_buf4_init(&W, CC2_HET1, 0, 2, 7, 2, 7, 0, "CC2 Wabij (i>j,a>b)");
    dpd_buf4_copy(&W, CC_TAMPS, "New tIJAB");
    dpd_buf4_copy(&W, CC_TAMPS, "New tijab");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC2_HET1, 0, 0, 5, 0, 5, 0, "CC2 WAbIj (Ij,Ab)");
    dpd_buf4_copy(&W, CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&W);

  }

  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&W, CC2_HET1, 0, 2, 7, 2, 7, 0, "CC2 WABIJ (I>J,A>B)");
    dpd_buf4_copy(&W, CC_TAMPS, "New tIJAB");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC2_HET1, 0, 12, 17, 12, 17, 0, "CC2 Wabij (i>j,a>b)");
    dpd_buf4_copy(&W, CC_TAMPS, "New tijab");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC2_HET1, 0, 22, 28, 22, 28, 0, "CC2 WAbIj (Ij,Ab)");
    dpd_buf4_copy(&W, CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&W);

  }

}
}} // namespace psi::ccenergy
