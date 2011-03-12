/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/

/*! \defgroup ccresponse ccresponse: Coupled-cluster response module */

#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double HXY(const char *pert_x, int irrep_x, double omega_x, 
	   const char *pert_y, int irrep_y, double omega_y)
{
  double polar;
  dpdfile2 X1, Y1, z;
  dpdbuf4 I;
  char lbl[32];

  sprintf(lbl, "Z_%s_IA", pert_y);
  dpd_file2_init(&z, CC_TMP0, irrep_y, 0, 1, lbl);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
  dpd_file2_init(&Y1, CC_OEI, irrep_y, 0, 1, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_dot24(&Y1, &I, &z, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_file2_close(&Y1);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  polar = 2.0 * dpd_file2_dot(&X1, &z);
  dpd_file2_close(&X1);

  dpd_file2_close(&z);

  return polar;
}

}} // namespace psi::ccresponse
