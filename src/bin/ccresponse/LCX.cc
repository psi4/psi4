/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double LCX(const char *pert_c, int irrep_c, 
	   const char *pert_x, int irrep_x, double omega)
{
  double polar=0.0;
  dpdfile2 X1, mu1, z1, l1, mu, lt, xc;
  dpdbuf4 X2, mu2, z2, l2, Z;
  char lbl[32];

  /*** Mu * X1 ***/

  sprintf(lbl, "%s_IA", pert_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  polar += 2.0 * dpd_file2_dot(&mu1, &X1);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  /*** L1 * MuBAR * X1 + L1 * MuBAR * X2 ***/

  dpd_file2_init(&z1, CC_TMP0, 0, 0, 1, "z_IA");

  sprintf(lbl, "%sBAR_MI", pert_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 0, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&mu1, &X1, &z1, 1, 1, -1, 0);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  sprintf(lbl, "%sBAR_AE", pert_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 1, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&X1, &mu1, &z1, 0, 0, 1, 1);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  sprintf(lbl, "%sBAR_ME", pert_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 1, lbl);
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert_x, omega);
  dpd_buf4_init(&X2, CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&mu1, &X2, &z1, 0, 0, 1, 1);
  dpd_buf4_close(&X2);
  dpd_file2_close(&mu1);

  dpd_file2_init(&l1, CC_LAMPS, 0, 0, 1, "LIA 0 -1");
  polar += 2.0 * dpd_file2_dot(&z1, &l1);
  dpd_file2_close(&l1);

  dpd_file2_close(&z1);

  /*** L2 * MuBAR * X1 + L2 * MuBAR * X2 ***/

  dpd_file2_init(&xc, CC_TMP0, 0, 0, 0, "XC_IJ");
  sprintf(lbl, "%s_IA", pert_c);
  dpd_file2_init(&mu, CC_OEI, irrep_c, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&X1, &mu, &xc, 0, 0, 2.0, 0.0);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu);
  dpd_file2_close(&xc);

  dpd_file2_init(&lt, CC_OEI, 0, 0, 0, "Lt_IJ");
  dpd_file2_init(&xc, CC_TMP0, 0, 0, 0, "XC_IJ");
  polar += dpd_file2_dot(&lt, &xc);
  dpd_file2_close(&xc);
  dpd_file2_close(&lt);

  dpd_buf4_init(&z2, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_scm(&z2, 0);
  dpd_buf4_close(&z2);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);

  dpd_buf4_init(&z2, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "%sBAR_MbIj", pert_c);
  dpd_buf4_init(&mu2, CC_LR, irrep_c, 10, 0, 10, 0, 0, lbl);
  dpd_contract244(&X1, &mu2, &Z, 0, 0, 1, 1, 0);
  dpd_buf4_close(&mu2);
  dpd_buf4_axpy(&Z, &z2, -1);
  dpd_buf4_sort(&Z, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z, &z2, -1);
  dpd_buf4_close(&Z);

  dpd_file2_close(&X1);

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega);
  dpd_buf4_init(&X2, CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");

  sprintf(lbl, "%sBAR_AE", pert_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 1, 1, lbl);
  dpd_contract424(&X2, &mu1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&mu1);
  dpd_buf4_axpy(&Z, &z2, 1);
  dpd_buf4_sort(&Z, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z, &z2, 1);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");

  sprintf(lbl, "%sBAR_MI", pert_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 0, lbl);
  dpd_contract244(&mu1, &X2, &Z, 0, 0, 0, 1, 0);
  dpd_file2_close(&mu1);
  dpd_buf4_axpy(&Z, &z2, -1);
  dpd_buf4_sort(&Z, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z, &z2, -1);
  dpd_buf4_close(&Z);

  dpd_buf4_close(&X2);

  dpd_buf4_init(&l2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  polar += dpd_buf4_dot(&l2, &z2);
  dpd_buf4_close(&l2);

  dpd_buf4_close(&z2);

  if(params.sekino) {  /* disconnected piece for Sekino-Bartlett modelIII */
    /* L2 * MUBAR * X1 */
    sprintf(lbl, "%sZ_IA", pert_c);
    dpd_file2_init(&z1, CC_TMP0, irrep_c, 0, 1, lbl);
    sprintf(lbl, "%sBAR_IA", pert_c);
    dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 1, lbl);
    dpd_buf4_init(&l2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    dpd_dot24(&mu1, &l2, &z1, 0, 0, 1, 0);
    dpd_buf4_close(&l2);
    dpd_file2_close(&mu1);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
    dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
    polar += 2.0 * dpd_file2_dot(&X1, &z1);
    dpd_file2_close(&X1);
    dpd_file2_close(&z1);
  }

  return polar;
}

}} // namespace psi::ccresponse
