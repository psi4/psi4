/*! \file
    \ingroup STABLE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace stable {

void diag_A_RHF(void)
{
  int h, dim, i;
  double *eps, **v;
  dpdbuf4 A, B;

  dpd_buf4_init(&A, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_buf4_init(&B, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");
  for(h=0; h < moinfo.nirreps; h++) {

    dim = A.params->rowtot[h];
    eps = init_array(dim);
    v = block_matrix(dim, dim);
    moinfo.rank[h] = dim;

    dpd_buf4_mat_irrep_init(&A, h);
    dpd_buf4_mat_irrep_rd(&A, h);
    sq_rsp(dim, dim, A.matrix[h], eps, 1, v, 1e-12);
    dpd_buf4_mat_irrep_close(&A, h);

    for(i=0; i < MIN0(dim, 5); i++)
      moinfo.A_evals[h][i] = eps[i];

    zero_mat(v, dim, dim);
    zero_arr(eps, dim);

    dpd_buf4_mat_irrep_init(&B, h);
    dpd_buf4_mat_irrep_rd(&B, h);
    sq_rsp(dim, dim, B.matrix[h], eps, 1, v, 1e-12);
    dpd_buf4_mat_irrep_close(&B, h);

    for(i=0; i < MIN0(dim, 5); i++)
      moinfo.A_triplet_evals[h][i] = eps[i];

    free_block(v);
    free(eps);

  }
  dpd_buf4_close(&B);
  dpd_buf4_close(&A);
}

}} // namespace psi::stable
