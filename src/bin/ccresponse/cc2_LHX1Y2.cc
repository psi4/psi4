/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double cc2_LHX1Y2(const char *pert_x, int irrep_x, double omega_x, 
	          const char *pert_y, int irrep_y, double omega_y)
{
  dpdfile2 z, z1, X1, l1, F;
  dpdbuf4 Z, Z1, Z2, I, Y2, L2, W;
  char lbl[32];
  double polar;
  int nirreps, Gbm, Gef, Gjf, Ge, Gf, Gj, bm, ef, jf;
  int *occpi, *virtpi, **W_col_offset, **Z_col_offset, offset;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;

  sprintf(lbl, "Z_%s_MI", pert_y);
  dpd_file2_init(&z1, CC_TMP0, irrep_y, 0, 0, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract442(&I, &Y2, &z1, 0, 0, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&I);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&z1, &X1, &z, 1, 1, -1, 0);
  dpd_file2_close(&X1);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);

  sprintf(lbl, "Z_%s_AE", pert_y);
  dpd_file2_init(&z1, CC_TMP0, irrep_y, 1, 1, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract442(&Y2, &I, &z1, 3, 3, -1, 0);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&I);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&X1, &z1, &z, 0, 0, 1, 1);
  dpd_file2_close(&X1);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);


  sprintf(lbl, "Z_%s_ME", pert_x);
  dpd_file2_init(&z1, CC_TMP0, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_dot24(&X1, &I, &z1, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_file2_close(&X1);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&z1, &Y2, &z, 0, 0, 1, 1);
  dpd_buf4_close(&Y2);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  dpd_file2_init(&l1, CC_LAMPS, 0, 0, 1, "LIA 0 -1");
  polar = 2.0 * dpd_file2_dot(&z, &l1);
  dpd_file2_close(&l1);
  dpd_file2_close(&z);

  /*  fprintf(outfile, "L(1)HX1Y2 = %20.12f\n", polar); */

  return polar;
}

}} // namespace psi::ccresponse
