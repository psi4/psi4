/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <physconst.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void diag(void)
{
  int h, i, dim, dim_A, dim_B;
  int nroot, root;
  int ai, bj, ck, c, k, C, K, Ksym;
  double **A, *eps, **v;
  char lbl[32];
  dpdbuf4 A_AA, A_BB, A_AB;
  dpdfile2 B;
  double **singlet_evals, **triplet_evals, **uhf_evals;

  if(params.ref == 0) { /**RHF **/

    singlet_evals = (double **) malloc(moinfo.nirreps * sizeof(double *));
    for(h=0; h < moinfo.nirreps; h++)
      singlet_evals[h] = init_array(params.rpi[h]);

    dpd_buf4_init(&A_AA, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");

    for(h=0; h < moinfo.nirreps; h++) {
      dim = A_AA.params->rowtot[h];
      eps = init_array(dim);
      v = block_matrix(dim, dim);

      dpd_buf4_mat_irrep_init(&A_AA, h);
      dpd_buf4_mat_irrep_rd(&A_AA, h);

      if(!strcmp(params.diag_method,"FULL")) 
	sq_rsp(dim, dim, A_AA.matrix[h], eps, 1, v, 1e-14);
      else if(!strcmp(params.diag_method,"DAVIDSON")) {
        nroot = david(A_AA.matrix[h], dim, params.rpi[h], eps, v, params.convergence, 0);

	if(nroot != params.rpi[h])
	  fprintf(outfile, "\tDavidson algorithm converged only %d roots in irrep %d.\n", 
		  nroot, h);
      }

      dpd_buf4_mat_irrep_close(&A_AA, h);

      /*
	fprintf(outfile, "RHF-CIS Singlet Eigenvectors for irrep %d:\n", h);
	print_mat(v, dim, params.rpi[h], outfile);
      */

      /* Store the eigenvectors in DPD entries and save the eigenvalues*/
      for(root=0; root < params.rpi[h]; root++) {

	sprintf(lbl, "BIA(%d)[%d] singlet", root, h);
	dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
	dpd_file2_mat_init(&B);
	for(ck=0; ck < dim; ck++) {
	  c = A_AA.params->roworb[h][ck][0];
	  k = A_AA.params->roworb[h][ck][1];

	  K = B.params->rowidx[k];
	  C = B.params->colidx[c];

	  Ksym = B.params->psym[k];

	  B.matrix[Ksym][K][C] = v[ck][root];
	}
	dpd_file2_mat_wrt(&B);
	dpd_file2_mat_close(&B);
	dpd_file2_close(&B);

	singlet_evals[h][root] = eps[root];
      }


      free(eps);
      free_block(v);
    }
    dpd_buf4_close(&A_AA);

    triplet_evals = (double **) malloc(moinfo.nirreps * sizeof(double *));
    for(h=0; h < moinfo.nirreps; h++)
      triplet_evals[h] = init_array(params.rpi[h]);

    dpd_buf4_init(&A_AA, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");

    for(h=0; h < moinfo.nirreps; h++) {
      dim = A_AA.params->rowtot[h];
      eps = init_array(dim);
      v = block_matrix(dim, dim);

      dpd_buf4_mat_irrep_init(&A_AA, h);
      dpd_buf4_mat_irrep_rd(&A_AA, h);

      if(!strcmp(params.diag_method, "FULL"))
	sq_rsp(dim, dim, A_AA.matrix[h], eps, 1, v, 1e-14);
      else if(!strcmp(params.diag_method,"DAVIDSON")) {
	nroot = david(A_AA.matrix[h], dim, params.rpi[h], eps, v, params.convergence, 0);

	if(nroot != params.rpi[h])
	  fprintf(outfile, "\tDavidson algorithm converged only %d roots in irrep %d.\n", 
		  nroot, h);
      }

      dpd_buf4_mat_irrep_close(&A_AA, h);

      /* Store the eigenvectors in DPD entries and save the eigenvalues*/
      for(root=0; root < params.rpi[h]; root++) {

	sprintf(lbl, "BIA(%d)[%d] triplet", root, h);
	dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
	dpd_file2_mat_init(&B);
	for(ck=0; ck < dim; ck++) {
	  c = A_AA.params->roworb[h][ck][0];
	  k = A_AA.params->roworb[h][ck][1];

	  K = B.params->rowidx[k];
	  C = B.params->colidx[c];

	  Ksym = B.params->psym[k];

	  B.matrix[Ksym][K][C] = v[ck][root];
	}
	dpd_file2_mat_wrt(&B);
	dpd_file2_mat_close(&B);
	dpd_file2_close(&B);

	triplet_evals[h][root] = eps[root];
      }

      free(eps);
      free_block(v);
    }
    dpd_buf4_close(&A_AA);

    moinfo.singlet_evals = singlet_evals;
    moinfo.triplet_evals = triplet_evals;

  }
  else if(params.ref == 2) { /** UHF **/

    uhf_evals = (double **) malloc(moinfo.nirreps * sizeof(double *));
    for(h=0; h < moinfo.nirreps; h++)
      uhf_evals[h] = init_array(params.rpi[h]);

    dpd_buf4_init(&A_AA, CC_MISC, 0, 21, 21, 21, 21, 0, "A(AI,BJ)");
    dpd_buf4_init(&A_BB, CC_MISC, 0, 31, 31, 31, 31, 0, "A(ai,bj)");
    dpd_buf4_init(&A_AB, CC_MISC, 0, 21, 31, 21, 31, 0, "A(AI,bj)");
    for(h=0; h < moinfo.nirreps; h++) {
      dim_A = A_AA.params->rowtot[h];
      dim_B = A_BB.params->rowtot[h];

      dim = dim_A + dim_B;

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

      if(!strcmp(params.diag_method,"FULL"))
	sq_rsp(dim, dim, A, eps, 1, v, 1e-12);
      else if(!strcmp(params.diag_method,"DAVIDSON")) {
	nroot = david(A, dim, params.rpi[h], eps, v, params.convergence, 0);

	if(nroot != params.rpi[h])
	  fprintf(outfile, "\tDavidson algorithm converged only %d roots in irrep %d.\n", 
		  nroot, h);
      }

      /*
	fprintf(outfile, "UHF-CIS Eigenvectors for irrep %d:\n", h);
	print_mat(v, dim, params.rpi[h], outfile);
      */

      /* Store the eigenvectors in DPD entries and save the eigenvalues */
      for(root=0; root < params.rpi[h]; root++) {
	sprintf(lbl, "BIA(%d)[%d]", root, h);
	dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
	dpd_file2_mat_init(&B);
	for(ck=0; ck < dim_A; ck++) {
	  c = A_AA.params->roworb[h][ck][0];
	  k = A_AA.params->roworb[h][ck][1];

	  K = B.params->rowidx[k];
	  C = B.params->colidx[c];

	  Ksym = B.params->psym[k];

	  B.matrix[Ksym][K][C] = v[ck][root];
	}
	dpd_file2_mat_wrt(&B);
	dpd_file2_mat_close(&B);
	dpd_file2_close(&B);

	sprintf(lbl, "Bia(%d)[%d]", root, h);
	dpd_file2_init(&B, CC_OEI, h, 2, 3, lbl);
	dpd_file2_mat_init(&B);
	for(ck=0; ck < dim_B; ck++) {
	  c = A_BB.params->roworb[h][ck][0];
	  k = A_BB.params->roworb[h][ck][1];

	  K = B.params->rowidx[k];
	  C = B.params->colidx[c];

	  Ksym = B.params->psym[k];

	  B.matrix[Ksym][K][C] = v[ck+dim_A][root];
	}
	dpd_file2_mat_wrt(&B);
	dpd_file2_mat_close(&B);
	dpd_file2_close(&B);

	uhf_evals[h][root] = eps[root];
      }

      free(eps);
      free_block(v);
      free_block(A);
    }
    dpd_buf4_close(&A_AA);
    dpd_buf4_close(&A_BB);
    dpd_buf4_close(&A_AB);

    moinfo.uhf_evals = uhf_evals;
  }

}

}} // namespace psi::cis
