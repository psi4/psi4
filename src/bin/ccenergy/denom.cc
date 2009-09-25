/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void local_filter_T1(dpdfile2 *T1);
void dijabT2(void);

/* apply denominators to t1 and t2 */

void denom(void)
{
  dpdfile2 newtIA, dIA, tIA, newtia, dia, tia;

  if (params.ref == 0) {
    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_copy(&newtIA, CC_OEI, "New tIA Increment");
    dpd_file2_close(&newtIA);

    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA Increment");
    if(params.local && local.filter_singles) {
      local_filter_T1(&newtIA);
    }
    else {
      dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &newtIA);
      dpd_file2_close(&dIA);
    }
    dpd_file2_close(&newtIA);

    /* Add the new increment to the old tIA to get the New tIA */
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_copy(&tIA, CC_OEI, "New tIA");
    dpd_file2_close(&tIA);
    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "New tIA Increment");
    dpd_file2_axpy(&tIA, &newtIA, 1, 0);
    dpd_file2_close(&tIA);
    dpd_file2_close(&newtIA);
  }
  else if (params.ref == 1) {
    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newtIA);
    dpd_file2_close(&dIA);
    dpd_file2_close(&newtIA); 

    dpd_file2_init(&newtia, CC_OEI, 0, 0, 1, "New tia");
    dpd_file2_init(&dia, CC_OEI, 0, 0, 1, "dia");
    dpd_file2_dirprd(&dia, &newtia);
    dpd_file2_close(&dia);
    dpd_file2_close(&newtia);
  }
  else if (params.ref == 2) {
    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newtIA);
    dpd_file2_close(&dIA);
    dpd_file2_close(&newtIA);

    dpd_file2_init(&newtia, CC_OEI, 0, 2, 3, "New tia");
    dpd_file2_init(&dia, CC_OEI, 0, 2, 3, "dia");
    dpd_file2_dirprd(&dia, &newtia);
    dpd_file2_close(&dia);
    dpd_file2_close(&newtia);
  }

  dijabT2(); 

  return;
}
}} // namespace psi::ccenergy
