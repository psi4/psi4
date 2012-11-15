/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/* Computes the D2 diagnostic as defined in:
 *
 * I.M.B. Nielsen and C.L. Janssen, Chem. Phys. Lett. 310, 568 (1999).
 *
 * */

double d2diag_rhf(void)
{
  int h, nirreps, i;
  double **Co, *Eo, max;
  double **Cv, *Ev;
  dpdbuf4 Tikab, Tjkab;
  dpdbuf4 Tijac, Tijbc;
  dpdfile2 To, Tv;

  nirreps = moinfo.nirreps;
  max = 0.0;

  dpd_buf4_init(&Tikab, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&Tjkab, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_file2_init(&To, CC_TMP0, 0, 0, 0, "To");
  dpd_buf4_init(&Tijac, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&Tijbc, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_file2_init(&Tv, CC_TMP0, 0, 1, 1, "Tv");

  // Build diagnostic matrices To and Tv //
  dpd_contract442(&Tikab, &Tjkab, &To, 0, 0, 1, 0);
  dpd_contract442(&Tijac, &Tijbc, &Tv, 3, 3, 1, 0);

  dpd_buf4_close(&Tikab);
  dpd_buf4_close(&Tjkab);
  dpd_file2_close(&To);
  dpd_buf4_close(&Tijac);
  dpd_buf4_close(&Tijbc);
  dpd_file2_close(&Tv);

  dpd_file2_init(&To, CC_TMP0, 0, 0, 0, "To");
  dpd_file2_mat_init(&To);
  dpd_file2_mat_rd(&To);

  dpd_file2_init(&Tv, CC_TMP0, 0, 1, 1, "Tv");
  dpd_file2_mat_init(&Tv);
  dpd_file2_mat_rd(&Tv);

  for(h=0; h < nirreps; h++) {
    if(To.params->rowtot[h]) {
      // Diagonalize To //
      Eo = init_array(To.params->rowtot[h]);
      Co = block_matrix(To.params->rowtot[h], To.params->rowtot[h]);
      sq_rsp(To.params->rowtot[h], To.params->rowtot[h], To.matrix[h], Eo, 0, Co, 1e-12);
      // Find maximum To eigenvalue //
      for(i=0; i < To.params->rowtot[h]; i++) {
        if(Eo[i] > max) max = Eo[i];
      }
    free_block(Co);
    free(Eo);
    }

    if(Tv.params->rowtot[h]) {
      // Diagonalize Tv //
      Ev = init_array(Tv.params->rowtot[h]);
      Cv = block_matrix(Tv.params->rowtot[h], Tv.params->rowtot[h]);
      sq_rsp(Tv.params->rowtot[h], Tv.params->rowtot[h], Tv.matrix[h], Ev, 0, Cv, 1e-12);

      // Find maximum Tv eigenvalue //
      for(i=0; i < Tv.params->rowtot[h]; i++) {
        if(Ev[i] > max) max = Ev[i];
      }

      free_block(Cv);
      free(Ev);
    }
  }

  dpd_file2_mat_close(&To);
  dpd_file2_mat_close(&Tv);
  dpd_file2_close(&To);
  dpd_file2_close(&Tv);

  /*
  // Original algorithm
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_sort(&T2, CC_TMP0, qrsp, 10, 11, "tjAbI");
  dpd_buf4_close(&T2);

//  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
//  dpd_buf4_sort(&T2, CC_TMP0, rpqs, 11, 10, "tAIjb");
//  dpd_buf4_sort(&T2, CC_TMP0, pqsr, 0, 5, "tIjbA");
//  dpd_buf4_close(&T2);

  dpd_buf4_init(&Tijab, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&Tjabi, CC_TMP0, 0, 10, 11, 10, 11, 0, "tjAbI");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Tijab, h);
    dpd_buf4_mat_irrep_rd(&Tijab, h);
//    dpd_buf4_mat_irrep_shift13(&Tijab, h);

    dpd_buf4_mat_irrep_init(&Tjabi, h);
    dpd_buf4_mat_irrep_rd(&Tjabi, h);
    dpd_bufg_mat_irrep_shift31(&Tjabi, h);

    nrows = moinfo.occpi[h];
    ncols = moinfo.occpi[h] * moinfo.virtpi[h] * moinfo.virtpi[h];
    To = block_matrix(nrows, nrows);
    C_DGEMM('n', 'n', nrows, nrows, ncols, 1.0, Tijab.shift.matrix[h][h][0],
            ncols, Tjabi.shift.matrix[h][h][0], nrows, 0.0, To[0], nrows);

    Eo = init_array(nrows);
    Co = block_matrix(nrows, nrows);
    sq_rsp(nrows, nrows, To, Eo, 0, Co, 1e-12);

    for(i=0; i < nrows; i++) {
      if(Eo[i] > max) max = Eo[i];
    }

    dpd_buf4_mat_irrep_close(&Tijab, h);
    dpd_buf4_mat_irrep_close(&Tjabi, h);

    free_block(To);
    free_block(Co);
    free(Eo);
  }

  dpd_buf4_close(&Tijab);
  dpd_buf4_close(&Tjabi);
  // END original algorithm
  */

  max = sqrt(max);

  return max;

}

double d2diag(void)
{
  double norm = 0.0;

  if(params.ref == 0) { /** RHF **/
    norm = d2diag_rhf();
  }
  return norm;

  // There's no open shell definitions, but I've set it up just incase -TJM
}
}} // namespace psi::ccenergy
