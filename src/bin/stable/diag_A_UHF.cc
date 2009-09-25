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

extern void follow_evec_UHF(double *vf, int dim_A, int dim_B);
extern void follow_evec_UHF2(double *vf, int dim_A, int dim_B);

void diag_A_UHF(void)
{
  int h, i, dim, dim_A, dim_B;
  int nroot, root;
  int ai, bj, ck, c, k, C, K, Ksym;
  double **A, *eps, **v, *evec_to_follow;
  char lbl[32];
  dpdbuf4 A_AA, A_BB, A_AB;
  dpdfile2 B;

  dpd_buf4_init(&A_AA, PSIF_MO_HESS, 0, 21, 21, 21, 21, 0, "A(AI,BJ)");
  dpd_buf4_init(&A_BB, PSIF_MO_HESS, 0, 31, 31, 31, 31, 0, "A(ai,bj)");
  dpd_buf4_init(&A_AB, PSIF_MO_HESS, 0, 21, 31, 21, 31, 0, "A(AI,bj)");
  for(h=0; h < moinfo.nirreps; h++) {
    dim_A = A_AA.params->rowtot[h];
    dim_B = A_BB.params->rowtot[h];

    dim = dim_A + dim_B;
    moinfo.rank[h] = dim;

    A = block_matrix(dim, dim);
    eps = init_array(dim);
    v = block_matrix(dim, dim);

    dpd_buf4_mat_irrep_init(&A_AA, h);
    dpd_buf4_mat_irrep_rd(&A_AA, h);
    for(ai=0; ai < dim_A; ai++)
      for(bj=0; bj < dim_A; bj++)
	A[ai][bj] = A_AA.matrix[h][ai][bj];
    dpd_buf4_mat_irrep_close(&A_AA, h);

    dpd_buf4_mat_irrep_init(&A_BB, h);
    dpd_buf4_mat_irrep_rd(&A_BB, h);
    for(ai=0; ai < dim_B; ai++)
      for(bj=0; bj < dim_B; bj++)
	A[ai+dim_A][bj+dim_A] = A_BB.matrix[h][ai][bj];
    dpd_buf4_mat_irrep_close(&A_BB, h);

    dpd_buf4_mat_irrep_init(&A_AB, h);
    dpd_buf4_mat_irrep_rd(&A_AB, h);
    for(ai=0; ai < dim_A; ai++)
      for(bj=0; bj < dim_B; bj++)
	A[ai][bj+dim_A] = A[bj+dim_A][ai] = A_AB.matrix[h][ai][bj];
    dpd_buf4_mat_irrep_close(&A_AB, h);

    sq_rsp(dim, dim, A, eps, 1, v, 1e-12);

    for(i=0; i < MIN0(moinfo.rank[h],5); i++)
      moinfo.A_evals[h][i] = eps[i];

    if (params.num_evecs_print > 0) {
      fprintf(outfile, "First %d eigenvectors for %s block:\n",
        MIN0(moinfo.rank[h], params.num_evecs_print), moinfo.labels[h]);
      eivout(v,eps,dim,MIN0(moinfo.rank[h],params.num_evecs_print),outfile); 
    }

    /* try to follow eigenvector downhill */
    if (h==0 && params.follow_instab) {
      if (eps[0] >= 0.0) 
        fprintf(outfile, "\nNo negative eigenvalues to follow.\n");
      else {
        evec_to_follow = init_array(dim);
        for (i=0; i<dim; i++) evec_to_follow[i] = v[i][0];

        fprintf(outfile, "\nAttempting to follow eigenvector using ");
        if (params.rotation_method == 0) {
          fprintf(outfile, "orbital rotation method.\n");
          follow_evec_UHF(evec_to_follow, dim_A, dim_B);
        }  
        else {
          fprintf(outfile, "antisymmetric matrix method.\n");
          follow_evec_UHF2(evec_to_follow, dim_A, dim_B);
        } 
        free(evec_to_follow);
      }
    }

    free(eps);
    free_block(v);
    free_block(A);
  }

  dpd_buf4_close(&A_AA);
  dpd_buf4_close(&A_BB);
  dpd_buf4_close(&A_AB);

}


}} // namespace psi::stable
