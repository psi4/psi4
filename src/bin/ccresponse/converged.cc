/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <cmath>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double converged(const char *pert, int irrep, double omega)
{
  dpdfile2 X1, X1new;
  dpdbuf4 X2, X2new;
  double rms=0.0, value;
  int row, col, h, nirreps; 
  char lbl[32];

  nirreps = moinfo.nirreps;

  sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_init(&X1new, CC_OEI, irrep, 0, 1, lbl);
  dpd_file2_mat_init(&X1new);
  dpd_file2_mat_rd(&X1new);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);

  for(h=0; h < nirreps; h++)
    for(row=0; row < X1.params->rowtot[h]; row++)
      for(col=0; col < X1.params->coltot[h^irrep]; col++) {
	value = X1new.matrix[h][row][col] - X1.matrix[h][row][col];
	rms += value * value;
      }
  dpd_file2_mat_close(&X1new);
  dpd_file2_close(&X1new);
  dpd_file2_mat_close(&X1);
  dpd_file2_close(&X1);

  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&X2new, h);
    dpd_buf4_mat_irrep_rd(&X2new, h);
    dpd_buf4_mat_irrep_init(&X2, h);
    dpd_buf4_mat_irrep_rd(&X2, h);

    for(row=0; row < X2.params->rowtot[h]; row++)
      for(col=0; col < X2.params->coltot[h^irrep]; col++) {
	value = X2new.matrix[h][row][col] - X2.matrix[h][row][col];
	rms += value * value;
      }

    dpd_buf4_mat_irrep_close(&X2new, h);
    dpd_buf4_mat_irrep_close(&X2, h);
  }
  dpd_buf4_close(&X2new);
  dpd_buf4_close(&X2);

  return sqrt(rms);
}

}} // namespace psi::ccresponse
