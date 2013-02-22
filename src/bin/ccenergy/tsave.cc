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

void tsave(void)
{
  dpdfile2 t1;
  dpdbuf4 t2;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_copy(&t1, PSIF_CC_OEI, "tIA");
    dpd_file2_close(&t1);

    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_copy(&t2, PSIF_CC_TAMPS, "tIjAb");
    dpd_buf4_close(&t2);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_copy(&t1, PSIF_CC_OEI, "tIA");
    dpd_file2_close(&t1);

    dpd_file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tia");
    dpd_file2_copy(&t1, PSIF_CC_OEI, "tia");
    dpd_file2_close(&t1);

    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_copy(&t2, PSIF_CC_TAMPS, "tIJAB");
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_buf4_copy(&t2, PSIF_CC_TAMPS, "tijab");
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_copy(&t2, PSIF_CC_TAMPS, "tIjAb");
    dpd_buf4_close(&t2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_copy(&t1, PSIF_CC_OEI, "tIA");
    dpd_file2_close(&t1);

    dpd_file2_init(&t1, PSIF_CC_OEI, 0, 2, 3, "New tia");
    dpd_file2_copy(&t1, PSIF_CC_OEI, "tia");
    dpd_file2_close(&t1);

    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_copy(&t2, PSIF_CC_TAMPS, "tIJAB");
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_buf4_copy(&t2, PSIF_CC_TAMPS, "tijab");
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_copy(&t2, PSIF_CC_TAMPS, "tIjAb");
    dpd_buf4_close(&t2);

  }
}
}} // namespace psi::ccenergy
