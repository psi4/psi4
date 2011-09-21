/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libint/libint.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

/*! 
** local_init(): Set up parameters of local excitation domains.
**
** The orbital domains constructed here are based on those described
** in Broughton and Pulay, J. Comp. Chem. 14, 736-740 (1993).  The
** localization of the occupied orbitals is done elsewhere (see the
** program "local").  Pair domains are defined as the union of pairs
** of single occupied orbital domains.  "Weak pairs", which are
** defined as pair domains whose individual occupied orbital domains
** have no atoms in common, are identified (cf. int *weak_pairs).
** 
** TDC, Jan-June 2002
*/

void local_init(void)
{
  local.nso = moinfo.nso;
  local.nocc = moinfo.occpi[0]; /* active doubly occupied orbitals */
  local.nvir = moinfo.virtpi[0]; /* active virtual orbitals */

  fprintf(outfile, "\tLocalization parameters ready.\n\n");
  fflush(outfile);
}

void local_done(void)
{
  fprintf(outfile, "\tLocal parameters free.\n");
}

void local_filter_T1(dpdfile2 *T1)
{
  int i, a, ij, ii;
  int nocc, nvir;
  double *T1tilde, *T1bar;
  psio_address next;

  nocc = local.nocc;
  nvir = local.nvir;

  local.pairdom_len = init_int_array(nocc*nocc);
  local.pairdom_nrlen = init_int_array(nocc*nocc);
  local.eps_occ = init_array(nocc);
  psio_read_entry(CC_INFO, "Local Pair Domain Length", (char *) local.pairdom_len,
		  nocc*nocc*sizeof(int));
  psio_read_entry(CC_INFO, "Local Pair Domain NR Length", (char *) local.pairdom_nrlen,
		  nocc*nocc*sizeof(int));
  psio_read_entry(CC_INFO, "Local Occupied Orbital Energies", (char *) local.eps_occ,
		  nocc*sizeof(double));

  local.W = (double ***) malloc(nocc * nocc * sizeof(double **));
  local.V = (double ***) malloc(nocc * nocc * sizeof(double **));
  local.eps_vir = (double **) malloc(nocc * nocc * sizeof(double *));
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++) {
    local.eps_vir[ij] = init_array(local.pairdom_nrlen[ij]);
    psio_read(CC_INFO, "Local Virtual Orbital Energies", (char *) local.eps_vir[ij],
	      local.pairdom_nrlen[ij]*sizeof(double), next, &next);
  }
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++) {
    local.V[ij] = block_matrix(nvir,local.pairdom_len[ij]);
    psio_read(CC_INFO, "Local Residual Vector (V)", (char *) local.V[ij][0],
	      nvir*local.pairdom_len[ij]*sizeof(double), next, &next);
  }
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++) {
    local.W[ij] = block_matrix(local.pairdom_len[ij],local.pairdom_nrlen[ij]);
    psio_read(CC_INFO, "Local Transformation Matrix (W)", (char *) local.W[ij][0],
	      local.pairdom_len[ij]*local.pairdom_nrlen[ij]*sizeof(double), next, &next);
  }

  dpd_file2_mat_init(T1);
  dpd_file2_mat_rd(T1);

  for(i=0; i < nocc; i++) {
    ii = i * nocc + i;  /* diagonal element of pair matrices */

    if(!local.pairdom_len[ii]) {
      throw PsiException("loca_filter_t1: pair is zero laength which makes no sense",__FILE__,__LINE__);
    }

    T1tilde = init_array(local.pairdom_len[ii]);
    T1bar = init_array(local.pairdom_nrlen[ii]);

    /* Transform the virtuals to the redundant projected virtual basis */
    C_DGEMV('t', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii], 
	    &(T1->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

    /* Transform the virtuals to the non-redundant virtual basis */
    C_DGEMV('t', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]), local.pairdom_nrlen[ii], 
	    &(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

    /* Apply the denominators */
    for(a=0; a < local.pairdom_nrlen[ii]; a++)
      T1bar[a] /= (local.eps_occ[i] - local.eps_vir[ii][a]);

    /* Transform the new T1's to the redundant projected virtual basis */
    C_DGEMV('n', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]), local.pairdom_nrlen[ii],
	    &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);


    /* Transform the new T1's to the MO basis */
    C_DGEMV('n', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii], 
	    &(T1tilde[0]), 1, 0.0, &(T1->matrix[0][i][0]), 1);

    free(T1bar);
    free(T1tilde);

  }

  dpd_file2_mat_wrt(T1);
  dpd_file2_mat_close(T1);

  for(i=0; i < nocc*nocc; i++) {
    free_block(local.W[i]);
    free_block(local.V[i]);
    free(local.eps_vir[i]);
  }
  free(local.W);
  free(local.V);
  free(local.eps_vir);

  free(local.eps_occ);
  free(local.pairdom_len);
  free(local.pairdom_nrlen);
}

void local_filter_T2(dpdbuf4 *T2)
{
  int ij, i, j, a, b;
  int nso, nocc, nvir;
  double **X1, **X2, **T2tilde, **T2bar;
  psio_address next;

  nso = local.nso;
  nocc = local.nocc;
  nvir = local.nvir;

  local.pairdom_len = init_int_array(nocc*nocc);
  local.pairdom_nrlen = init_int_array(nocc*nocc);
  local.weak_pairs = init_int_array(nocc*nocc);
  local.eps_occ = init_array(nocc);
  psio_read_entry(CC_INFO, "Local Pair Domain Length", (char *) local.pairdom_len,
		  nocc*nocc*sizeof(int));
  psio_read_entry(CC_INFO, "Local Pair Domain NR Length", (char *) local.pairdom_nrlen,
		  nocc*nocc*sizeof(int));
  psio_read_entry(CC_INFO, "Local Occupied Orbital Energies", (char *) local.eps_occ,
		  nocc*sizeof(double));
  psio_read_entry(CC_INFO, "Local Weak Pairs", (char *) local.weak_pairs,
		  nocc*nocc*sizeof(int));
  local.W = (double ***) malloc(nocc * nocc * sizeof(double **));
  local.V = (double ***) malloc(nocc * nocc * sizeof(double **));
  local.eps_vir = (double **) malloc(nocc * nocc * sizeof(double *));
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++) {
    local.eps_vir[ij] = init_array(local.pairdom_nrlen[ij]);
    psio_read(CC_INFO, "Local Virtual Orbital Energies", (char *) local.eps_vir[ij],
	      local.pairdom_nrlen[ij]*sizeof(double), next, &next);
  }
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++) {
    local.V[ij] = block_matrix(nvir,local.pairdom_len[ij]);
    psio_read(CC_INFO, "Local Residual Vector (V)", (char *) local.V[ij][0],
	      nvir*local.pairdom_len[ij]*sizeof(double), next, &next);
  }
  next = PSIO_ZERO;
  for(ij=0; ij < nocc*nocc; ij++) {
    local.W[ij] = block_matrix(local.pairdom_len[ij],local.pairdom_nrlen[ij]);
    psio_read(CC_INFO, "Local Transformation Matrix (W)", (char *) local.W[ij][0],
	      local.pairdom_len[ij]*local.pairdom_nrlen[ij]*sizeof(double), next, &next);
  }

  /* Grab the MO-basis T2's */
  dpd_buf4_mat_irrep_init(T2, 0);
  dpd_buf4_mat_irrep_rd(T2, 0);

  X1 = block_matrix(nso,nvir);
  X2 = block_matrix(nvir,nso);
  T2tilde = block_matrix(nso,nso);
  T2bar = block_matrix(nvir, nvir);

  for(i=0,ij=0; i < nocc; i++) {
    for(j=0; j < nocc; j++,ij++) {

      if(!local.weak_pairs[ij]) {

	/* Transform the virtuals to the redundant projected virtual basis */
	C_DGEMM('t', 'n', local.pairdom_len[ij], nvir, nvir, 1.0, &(local.V[ij][0][0]), local.pairdom_len[ij],
		&(T2->matrix[0][ij][0]), nvir, 0.0, &(X1[0][0]), nvir);
	C_DGEMM('n', 'n', local.pairdom_len[ij], local.pairdom_len[ij], nvir, 1.0, &(X1[0][0]), nvir,
		&(local.V[ij][0][0]), local.pairdom_len[ij], 0.0, &(T2tilde[0][0]), nso);

	/* Transform the virtuals to the non-redundant virtual basis */
	C_DGEMM('t', 'n', local.pairdom_nrlen[ij], local.pairdom_len[ij], local.pairdom_len[ij], 1.0, 
		&(local.W[ij][0][0]), local.pairdom_nrlen[ij], &(T2tilde[0][0]), nso, 0.0, &(X2[0][0]), nso);
	C_DGEMM('n', 'n', local.pairdom_nrlen[ij], local.pairdom_nrlen[ij], local.pairdom_len[ij], 1.0, 
		&(X2[0][0]), nso, &(local.W[ij][0][0]), local.pairdom_nrlen[ij], 0.0, &(T2bar[0][0]), nvir);

	/* Divide the new amplitudes by the denominators */
	for(a=0; a < local.pairdom_nrlen[ij]; a++) {
	  for(b=0; b < local.pairdom_nrlen[ij]; b++) {
	    T2bar[a][b] /= (local.eps_occ[i] + local.eps_occ[j]
			    - local.eps_vir[ij][a] - local.eps_vir[ij][b]);
	  }
	}

	/* Transform the new T2's to the redundant virtual basis */
	C_DGEMM('n', 'n', local.pairdom_len[ij], local.pairdom_nrlen[ij], local.pairdom_nrlen[ij], 1.0, 
		&(local.W[ij][0][0]), local.pairdom_nrlen[ij], &(T2bar[0][0]), nvir, 0.0, &(X1[0][0]), nvir);
	C_DGEMM('n','t', local.pairdom_len[ij], local.pairdom_len[ij], local.pairdom_nrlen[ij], 1.0, 
		&(X1[0][0]), nvir, &(local.W[ij][0][0]), local.pairdom_nrlen[ij], 0.0, &(T2tilde[0][0]), nso);

	/* Transform the new T2's to the MO basis */
	C_DGEMM('n', 'n', nvir, local.pairdom_len[ij], local.pairdom_len[ij], 1.0, 
		&(local.V[ij][0][0]), local.pairdom_len[ij], &(T2tilde[0][0]), nso, 0.0, &(X2[0][0]), nso);
	C_DGEMM('n', 't', nvir, nvir, local.pairdom_len[ij], 1.0, &(X2[0][0]), nso,
		&(local.V[ij][0][0]), local.pairdom_len[ij], 0.0, &(T2->matrix[0][ij][0]), nvir);
      }
      else  /* This must be a neglected weak pair; force it to zero */
	memset((void *) T2->matrix[0][ij], 0, nvir*nvir*sizeof(double));

    }
  }
  free_block(X1);
  free_block(X2);
  free_block(T2tilde);
  free_block(T2bar);

  /* Write the updated MO-basis T2's to disk */
  dpd_buf4_mat_irrep_wrt(T2, 0);
  dpd_buf4_mat_irrep_close(T2, 0);

  for(i=0; i < nocc*nocc; i++) {
    free_block(local.W[i]);
    free_block(local.V[i]);
    free(local.eps_vir[i]);
  }
  free(local.W);
  free(local.V);
  free(local.eps_vir);

  free(local.eps_occ);
  free(local.pairdom_len);
  free(local.pairdom_nrlen);
  free(local.weak_pairs);
}

}} // namespace psi::ccresponse
