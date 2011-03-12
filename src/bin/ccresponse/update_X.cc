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

void update_X(const char *pert, int irrep, double omega)
{
  dpdfile2 X1new, X1;
  dpdbuf4 X2new, X2;
  char lbl[32];

  sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_init(&X1new, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  dpd_file2_axpy(&X1, &X1new, 1, 0);
  dpd_file2_close(&X1);
  dpd_file2_close(&X1new);

  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_axpy(&X2, &X2new, 1);
  dpd_buf4_close(&X2);
  dpd_buf4_close(&X2new);
}

}} // namespace psi::ccresponse
