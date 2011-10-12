/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <exception.h>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* build_Z_UHF():  Solve the orbital Z-vector equations for UHF refs:
**
**    A(AI,BJ) D(orb)(B,J) = -X(A,I)
**
** where A(AI,EM) is the orbital Hessian computed in build_A(), X(A,I)
** is the orbital rotation gradient computed in build_X(), and
** D(orb)(E,M) is the final Z-vector we want.
**
*/

void build_Z_UHF(void)
{
  dpdbuf4 A_AA, A_BB, A_AB;
  dpdfile2 X, D;
  double **A, *Z;
  int num_ai, a, i, ai, bj;
  int h, nirreps, count, dim_A, dim_B;
  int *ipiv, info;
  int *avirtpi, *aoccpi;
  int *bvirtpi, *boccpi;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi; avirtpi = moinfo.avirtpi;
  boccpi = moinfo.boccpi; bvirtpi = moinfo.bvirtpi;

  /* compute the number of ai pairs */
  num_ai = 0;
  for(h=0; h < nirreps; h++) {
    num_ai += avirtpi[h] * aoccpi[h];
    num_ai += bvirtpi[h] * boccpi[h];
  }

  /* Place all the elements of the orbital rotation gradient, X, into a
     linear array, Z */
  Z = init_array(num_ai);

  dpd_file2_init(&X, CC_OEI, 0, 1, 0, "XAI");
  dpd_file2_mat_init(&X);
  dpd_file2_mat_rd(&X);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < X.params->rowtot[h]; a++)
      for(i=0; i < X.params->coltot[h]; i++)
        Z[count++] = -X.matrix[h][a][i];
  dpd_file2_mat_close(&X);
  dpd_file2_close(&X);

  dpd_file2_init(&X, CC_OEI, 0, 3, 2, "Xai");
  dpd_file2_mat_init(&X);
  dpd_file2_mat_rd(&X);
  for(h=0; h < nirreps; h++)
    for(a=0; a < X.params->rowtot[h]; a++)
      for(i=0; i < X.params->coltot[h]; i++)
        Z[count++] = -X.matrix[h][a][i];
  dpd_file2_mat_close(&X);
  dpd_file2_close(&X);

  /* Now, build the full MO Hessian */
  dpd_buf4_init(&A_AA, CC_MISC, 0, 21, 21, 21, 21, 0, "A(AI,BJ)");
  dpd_buf4_init(&A_BB, CC_MISC, 0, 31, 31, 31, 31, 0, "A(ai,bj)");
  dpd_buf4_init(&A_AB, CC_MISC, 0, 21, 31, 21, 31, 0, "A(AI,bj)");


  dim_A = A_AA.params->rowtot[0];
  dim_B = A_BB.params->rowtot[0];

  if(num_ai != dim_A + dim_B) { /* error */
    fprintf(outfile, "Problem: num_ai(%d) != dim_A + dim_b (%d)\n", num_ai,
            dim_A + dim_B);
    throw PsiException("ccenergy: error", __FILE__, __LINE__);
  }

  A = block_matrix(num_ai, num_ai);

  dpd_buf4_mat_irrep_init(&A_AA, 0);
  dpd_buf4_mat_irrep_rd(&A_AA, 0);
  for(ai=0; ai < dim_A; ai++)
    for(bj=0; bj < dim_A; bj++)
      A[ai][bj] = A_AA.matrix[0][ai][bj];
  dpd_buf4_mat_irrep_close(&A_AA, 0);

  dpd_buf4_mat_irrep_init(&A_BB, 0);
  dpd_buf4_mat_irrep_rd(&A_BB, 0);
  for(ai=0; ai < dim_B; ai++)
    for(bj=0; bj < dim_B; bj++)
      A[ai+dim_A][bj+dim_A] = A_BB.matrix[0][ai][bj];
  dpd_buf4_mat_irrep_close(&A_BB, 0);

  dpd_buf4_mat_irrep_init(&A_AB, 0);
  dpd_buf4_mat_irrep_rd(&A_AB, 0);
  for(ai=0; ai < dim_A; ai++)
    for(bj=0; bj < dim_B; bj++)
      A[ai][bj+dim_A] = A[bj+dim_A][ai] = A_AB.matrix[0][ai][bj];
  dpd_buf4_mat_irrep_close(&A_AB, 0);

  dpd_buf4_close(&A_AA);
  dpd_buf4_close(&A_BB);
  dpd_buf4_close(&A_AB);

  /*
  ipiv = init_int_array(num_ai);
  info = C_DGESV(num_ai, 1, &(A[0][0]), num_ai, &(ipiv[0]), &(Z[0]), num_ai);
  if(info) {
    fprintf(outfile, "\nError in DGESV return in build_Z_UHF: %d.\n", info);
    exit(PSI_RETURN_FAILURE);
  }

  free(ipiv);
  free_block(A);
  */
  pople(A, Z, num_ai, 1, 1e-12, outfile, 0);

  /*
  for(ai=0; ai < num_ai; ai++) fprintf(outfile, "Z[%d] = %20.15f\n", ai, Z[ai]);
  */

  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++)
      for(i=0; i < D.params->coltot[h]; i++) {
        if(fabs(Z[count]) > 1e3) D.matrix[h][a][i] = 0.0;
        else D.matrix[h][a][i] = Z[count];
        count++;
      }
  dpd_file2_mat_wrt(&D);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 3, 2, "D(orb)(a,i)");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_mat_init(&D);
  for(h=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++)
      for(i=0; i < D.params->coltot[h]; i++) {
        if(fabs(Z[count]) > 1e3) D.matrix[h][a][i] = 0.0;
        else D.matrix[h][a][i] = Z[count];
        count++;
      }
  dpd_file2_mat_wrt(&D);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* We're done with Z */
  free(Z);
}



}} // namespace psi::ccdensity
