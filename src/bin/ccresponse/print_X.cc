/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void print_X(char *pert, int irrep, double omega)
{
  dpdfile2 X1;
  dpdbuf4 X2;
  char lbl[32];

  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  dpd_file2_print(&X1, outfile);
  dpd_file2_close(&X1);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_print(&X2, outfile, 1);
  dpd_buf4_close(&X2);
}

}} // namespace psi::ccresponse
