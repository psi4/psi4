/*! \file
    \ingroup ccresponse
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

void cc2_sort_X(const char *pert, int irrep, double omega)
{
  dpdbuf4 X;
  char lbl[32];

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  dpd_buf4_scmcopy(&X, CC_LR, lbl, 2);
  dpd_buf4_sort_axpy(&X, CC_LR, pqsr, 0, 5, lbl, -1);
  dpd_buf4_close(&X);
}

}} // namespace psi::ccresponse
