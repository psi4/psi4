/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

double LRi_dot(int IRR, int R_index);
void LRi_minus(int IRR, int R_index, double overlap, double R0);

/* this function orthogonalizes the current L vector
  against the previously converged R vectors.
  This really helps for multiple excited-state L's.
  - ROHF and UHF still need to be added */

void ortho_Rs(struct L_Params *pL_params, int current_L) {
  int L_state_index, L_root, L_irr;
  int R_state_index, R_root, R_irr;
  double **O, tval, overlap;
  int L, R;

  if (params.ref != 0) return;

  L_irr  = pL_params[current_L].irrep;
  L_root = pL_params[current_L].root;
  
  for (R=1; R<params.nstates; ++R) {
    if (R == current_L) continue;
    R_irr  = pL_params[R].irrep;
    R_root = pL_params[R].root;

    if (L_irr != R_irr) continue;

    if (params.ref == 0) 
      overlap = LRi_dot(L_irr, R_root);

    if (L_root == -1)
      overlap += pL_params[R].R0;

     /* fprintf(outfile,"overlap with R[%d][%d]: %15.10lf\n", R_irr, R_root, overlap);  */
    LRi_minus(L_irr, R_root, overlap, pL_params[R].R0);
   /* overlap = LRi_dot(L_irr, R_root);
    if (L_root == -1)
      overlap += pL_params[R].R0;
    fprintf(outfile,"overlap with R[%d][%d]: %15.10lf\n", R_irr, R_root, overlap);  */
  }
  return;
}

double LRi_dot(int IRR, int R_index) {
  dpdfile2 R1, L1;
  dpdbuf4 R2, L2;
  double overlap;
  char R1A_lbl[32], lbl[32];

  sprintf(R1A_lbl, "RIA %d %d", IRR, R_index);
  dpd_file2_init(&R1, CC_RAMPS, IRR, 0, 1, R1A_lbl);
  dpd_file2_init(&L1, CC_LAMBDA, IRR, 0, 1, "New LIA");
  overlap = 2.0 * dpd_file2_dot(&L1, &R1);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  sprintf(lbl, "2RIjAb - RIjbA %d %d", IRR, R_index);
  dpd_buf4_init(&R2, CC_RAMPS, IRR, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&L2, CC_LAMBDA, IRR, 0, 5, 0, 5, 0, "New LIjAb");
  overlap += dpd_buf4_dot(&L2, &R2);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&R2);

  return overlap;
}

void LRi_minus(int IRR, int R_index, double overlap, double R0) {
  dpdfile2 R1, L1;
  dpdbuf4 R2, L2;
  char L1A_lbl[32], R1A_lbl[32], lbl[32];

  sprintf(R1A_lbl, "RIA %d %d", IRR, R_index);
  dpd_file2_init(&R1, CC_RAMPS, IRR, 0, 1, R1A_lbl);
  dpd_file2_init(&L1, CC_LAMBDA, IRR, 0, 1, "New LIA");
  dpd_file2_axpy(&R1, &L1, -overlap/(1.0 - R0*R0), 0);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  sprintf(lbl, "RIjAb %d %d", IRR, R_index);
  dpd_buf4_init(&R2, CC_RAMPS, IRR, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&L2, CC_LAMBDA, IRR, 0, 5, 0, 5, 0, "New LIjAb");
  dpd_buf4_axpy(&R2, &L2, -overlap/(1.0 - R0*R0));
  dpd_buf4_close(&L2);
  dpd_buf4_close(&R2);

  dpd_file2_init(&L1, CC_LAMBDA, IRR, 0, 1, "New LIA");
  dpd_file2_copy(&L1, CC_LAMBDA, "New Lia");
  dpd_file2_close(&L1);
  /*
  dpd_buf4_init(&L2, CC_LAMBDA, IRR, 2, 7, 0, 5, 1, "New LIjAb");
  dpd_buf4_copy(&L2, CC_LAMBDA, "New LIJAB");
  dpd_buf4_copy(&L2, CC_LAMBDA, "New Lijab");
  dpd_buf4_close(&L2);
  */
  return;
}


}} // namespace psi::cclambda
