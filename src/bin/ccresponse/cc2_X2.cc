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

namespace psi { namespace CCRESPONSE {

void denom2(dpdbuf4 *X2, double omega);
void local_filter_T2(dpdbuf4 *T2);

void cc2_X2_build(const char *pert, int irrep, double omega)
{
  dpdfile2 X1, z, F, t1;
  dpdbuf4 X2, X2new, Z, Z1, Z2, W, I;
  char lbl[32];

  sprintf(lbl, "%sBAR_IjAb", pert);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_copy(&X2new, CC_LR, lbl);
  dpd_buf4_close(&X2new);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  /*** D-S ***/

  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&W, CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
  dpd_contract244(&X1, &W, &Z, 0, 0, 1, 1, 0);
  dpd_buf4_close(&W);
  dpd_buf4_axpy(&Z, &X2new, -1);
  dpd_buf4_close(&X2new);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, -1);
  dpd_buf4_close(&Z);


  sprintf(lbl, "Z(Ab,Ij) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 5, 0, 5, 0, 0, lbl);
  dpd_buf4_init(&W, CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
  dpd_contract244(&X1, &W, &Z, 1, 2, 1, 1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, rspq, 0, 5, lbl, 1);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_sort(&Z, CC_TMP0, srqp, 0, 5, lbl);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_axpy(&Z, &X2new, 1);
  dpd_buf4_close(&X2new);
  dpd_buf4_close(&Z);

  dpd_file2_close(&X1);

  /*** D-D ***/

  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  dpd_buf4_axpy(&X2, &X2new, -omega);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
  dpd_contract424(&X2, &F, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_buf4_axpy(&Z, &X2new, 1);
  sprintf(lbl, "Z(jI,bA) %s", pert);
  dpd_buf4_sort(&Z, CC_TMP0, qpsr, 0, 5, lbl);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_axpy(&Z, &X2new, 1);
  dpd_buf4_close(&Z);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
  dpd_contract244(&F, &X2, &Z, 0, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_buf4_axpy(&Z, &X2new, -1);
  sprintf(lbl, "Z(jI,bA) %s", pert);
  dpd_buf4_sort(&Z, CC_TMP0, qpsr, 0, 5, lbl);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_axpy(&Z, &X2new, -1);
  dpd_buf4_close(&Z);

  dpd_buf4_close(&X2);

  /** Filter and apply denominator **/
  if(params.local) local_filter_T2(&X2new);
  else denom2(&X2new, omega);
  dpd_buf4_close(&X2new);
}

}} // namespace psi::CCRESPONSE
