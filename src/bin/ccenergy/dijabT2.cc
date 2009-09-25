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

void local_filter_T2(dpdbuf4 *T2);

void dijabT2(void)
{
  dpdbuf4 newtIJAB, newtijab, newtIjAb, tIjAb;
  dpdbuf4 dIJAB, dijab, dIjAb;

  if(params.ref == 0) { /*** RHF ***/
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_copy(&newtIjAb, CC_TAMPS, "New tIjAb Increment");
    dpd_buf4_close(&newtIjAb);

    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb Increment");
    if(params.local) {
      local_filter_T2(&newtIjAb);
    }
    else {
      dpd_buf4_init(&dIjAb, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
      dpd_buf4_dirprd(&dIjAb, &newtIjAb);
      dpd_buf4_close(&dIjAb);
    }
    dpd_buf4_close(&newtIjAb);

    /* Add the new increment to the old tIjAb to get the new tIjAb */
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_copy(&tIjAb, CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb Increment");
    dpd_buf4_axpy(&tIjAb, &newtIjAb, 1);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /*** ROHF ***/
    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&dIJAB, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
    dpd_buf4_dirprd(&dIJAB, &newtIJAB);
    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&dIJAB);

    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_buf4_init(&dijab, CC_DENOM, 0, 1, 6, 1, 6, 0, "dijab");
    dpd_buf4_dirprd(&dijab, &newtijab);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&dijab);

    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&dIjAb, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    dpd_buf4_dirprd(&dIjAb, &newtIjAb);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_close(&dIjAb);
  }
  else if(params.ref ==2) { /*** UHF ***/
    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&dIJAB, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
    dpd_buf4_dirprd(&dIJAB, &newtIJAB);
    dpd_buf4_close(&dIJAB);
    /*    dpd_buf4_print(&newtIJAB, outfile, 1); */
    dpd_buf4_close(&newtIJAB);

    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_buf4_init(&dijab, CC_DENOM, 0, 11, 16, 11, 16, 0, "dijab");
    dpd_buf4_dirprd(&dijab, &newtijab);
    dpd_buf4_close(&dijab);
    /*    dpd_buf4_print(&newtijab, outfile, 1); */
    dpd_buf4_close(&newtijab);

    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_init(&dIjAb, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
    dpd_buf4_dirprd(&dIjAb, &newtIjAb);
    dpd_buf4_close(&dIjAb);
    /*    dpd_buf4_print(&newtIjAb, outfile, 1); */
    dpd_buf4_close(&newtIjAb);
  }
}
}} // namespace psi::ccenergy
