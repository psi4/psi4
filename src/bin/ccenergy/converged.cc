/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

int converged(void)
{
  int row,col,h,nirreps;
  double rms=0.0;
  dpdfile2 T1, T1old;
  dpdbuf4 T2, T2old;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    dpd_file2_init(&T1old, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1old);
    dpd_file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
    dpd_file2_mat_close(&T1old);
    dpd_file2_close(&T1old);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&T2old, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      dpd_buf4_mat_irrep_init(&T2old, h);
      dpd_buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&T2, h);
      dpd_buf4_mat_irrep_close(&T2old, h);
    }
    dpd_buf4_close(&T2old);
    dpd_buf4_close(&T2);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    dpd_file2_init(&T1old, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1old);
    dpd_file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
    dpd_file2_mat_close(&T1old);
    dpd_file2_close(&T1old);

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "New tia");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    dpd_file2_init(&T1old, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_mat_init(&T1old);
    dpd_file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
    dpd_file2_mat_close(&T1old);
    dpd_file2_close(&T1old);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&T2old, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      dpd_buf4_mat_irrep_init(&T2old, h);
      dpd_buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&T2, h);
      dpd_buf4_mat_irrep_close(&T2old, h);
    }
    dpd_buf4_close(&T2old);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_buf4_init(&T2old, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      dpd_buf4_mat_irrep_init(&T2old, h);
      dpd_buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&T2, h);
      dpd_buf4_mat_irrep_close(&T2old, h);
    }
    dpd_buf4_close(&T2old);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&T2old, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      dpd_buf4_mat_irrep_init(&T2old, h);
      dpd_buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&T2, h);
      dpd_buf4_mat_irrep_close(&T2old, h);
    }
    dpd_buf4_close(&T2old);
    dpd_buf4_close(&T2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    dpd_file2_init(&T1old, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1old);
    dpd_file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
    dpd_file2_mat_close(&T1old);
    dpd_file2_close(&T1old);

    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "New tia");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    dpd_file2_init(&T1old, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_mat_init(&T1old);
    dpd_file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
    dpd_file2_mat_close(&T1old);
    dpd_file2_close(&T1old);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&T2old, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      dpd_buf4_mat_irrep_init(&T2old, h);
      dpd_buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&T2, h);
      dpd_buf4_mat_irrep_close(&T2old, h);
    }
    dpd_buf4_close(&T2old);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_buf4_init(&T2old, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      dpd_buf4_mat_irrep_init(&T2old, h);
      dpd_buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&T2, h);
      dpd_buf4_mat_irrep_close(&T2old, h);
    }
    dpd_buf4_close(&T2old);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_init(&T2old, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      dpd_buf4_mat_irrep_init(&T2old, h);
      dpd_buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&T2, h);
      dpd_buf4_mat_irrep_close(&T2old, h);
    }
    dpd_buf4_close(&T2old);
    dpd_buf4_close(&T2);
  }

  rms = sqrt(rms);
  moinfo.conv = rms;

  if(rms < params.convergence) return 1;
  else return 0;
}
}} // namespace psi::ccenergy
