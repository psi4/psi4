/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* for ROHF and UHF */
void precondition(dpdfile2 *RIA, dpdfile2 *Ria,
  dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb,
  double eval)
{
  dpdfile2 DIA, Dia;
  dpdbuf4 DIJAB, Dijab, DIjAb;
  int h, nirreps, i, j, a, b, ij, ab, C_irr;
  double tval;

  C_irr = RIA->my_irrep;
  nirreps = RIA->params->nirreps;

  global_dpd_->file2_mat_init(RIA);
  global_dpd_->file2_mat_rd(RIA);
  global_dpd_->file2_init(&DIA, PSIF_EOM_D, C_irr, 0, 1, "DIA");
  global_dpd_->file2_mat_init(&DIA);
  global_dpd_->file2_mat_rd(&DIA);
  for(h=0; h < nirreps; h++)
     for(i=0; i < RIA->params->rowtot[h]; i++)
        for(a=0; a < RIA->params->coltot[h^C_irr]; a++) {
           tval = eval - DIA.matrix[h][i][a];
           if (fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
        }
  global_dpd_->file2_mat_wrt(RIA);
  global_dpd_->file2_mat_close(RIA);
  global_dpd_->file2_mat_close(&DIA);
  global_dpd_->file2_close(&DIA);

  global_dpd_->file2_mat_init(Ria);
  global_dpd_->file2_mat_rd(Ria);
  if (params.eom_ref == 1)
    global_dpd_->file2_init(&Dia, PSIF_EOM_D, C_irr, 0, 1, "Dia");
  else if (params.eom_ref == 2)
    global_dpd_->file2_init(&Dia, PSIF_EOM_D, C_irr, 2, 3, "Dia");
  global_dpd_->file2_mat_init(&Dia);
  global_dpd_->file2_mat_rd(&Dia);
  for(h=0; h < nirreps; h++)
     for(i=0; i < Ria->params->rowtot[h]; i++)
        for(a=0; a < Ria->params->coltot[h^C_irr]; a++) {
           tval = eval - Dia.matrix[h][i][a];
           if (fabs(tval) > 0.0001) Ria->matrix[h][i][a] /= tval;
        }
  global_dpd_->file2_mat_wrt(Ria);
  global_dpd_->file2_mat_close(Ria);
  global_dpd_->file2_mat_close(&Dia);
  global_dpd_->file2_close(&Dia);


  global_dpd_->buf4_init(&DIJAB, PSIF_EOM_D, C_irr, 2, 7, 2, 7, 0, "DIJAB");
  for(h=0; h < RIJAB->params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(RIJAB, h);
    global_dpd_->buf4_mat_irrep_init(&DIJAB, h);
    global_dpd_->buf4_mat_irrep_rd(RIJAB, h);
    global_dpd_->buf4_mat_irrep_rd(&DIJAB, h);
    for(ij=0; ij < RIJAB->params->rowtot[h]; ij++)
       for(ab=0; ab < RIJAB->params->coltot[h^C_irr]; ab++) {
           tval = eval - DIJAB.matrix[h][ij][ab];
           if (fabs(tval) > 0.0001) RIJAB->matrix[h][ij][ab] /= tval;
      }
    global_dpd_->buf4_mat_irrep_wrt(RIJAB, h);
    global_dpd_->buf4_mat_irrep_close(RIJAB, h);
    global_dpd_->buf4_mat_irrep_close(&DIJAB, h);
  }
  global_dpd_->buf4_close(&DIJAB);


  if (params.eom_ref == 1)
    global_dpd_->buf4_init(&Dijab, PSIF_EOM_D, C_irr, 2, 7, 2, 7, 0, "Dijab");
  else if (params.eom_ref == 2)
    global_dpd_->buf4_init(&Dijab, PSIF_EOM_D, C_irr, 12, 17, 12, 17, 0, "Dijab");
  for(h=0; h < Rijab->params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(Rijab, h);
    global_dpd_->buf4_mat_irrep_init(&Dijab, h);
    global_dpd_->buf4_mat_irrep_rd(Rijab, h);
    global_dpd_->buf4_mat_irrep_rd(&Dijab, h);
    for(ij=0; ij < Rijab->params->rowtot[h]; ij++)
       for(ab=0; ab < Rijab->params->coltot[h^C_irr]; ab++) {
           tval = eval - Dijab.matrix[h][ij][ab];
           if (fabs(tval) > 0.0001) Rijab->matrix[h][ij][ab] /= tval;
      }
    global_dpd_->buf4_mat_irrep_wrt(Rijab, h);
    global_dpd_->buf4_mat_irrep_close(Rijab, h);
    global_dpd_->buf4_mat_irrep_close(&Dijab, h);
  }
  global_dpd_->buf4_close(&Dijab);


  if (params.eom_ref == 1)
    global_dpd_->buf4_init(&DIjAb, PSIF_EOM_D, C_irr, 0, 5, 0, 5, 0, "DIjAb");
  else if (params.eom_ref == 2)
    global_dpd_->buf4_init(&DIjAb, PSIF_EOM_D, C_irr, 22, 28, 22, 28, 0, "DIjAb");
  for(h=0; h < RIjAb->params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(RIjAb, h);
    global_dpd_->buf4_mat_irrep_init(&DIjAb, h);
    global_dpd_->buf4_mat_irrep_rd(RIjAb, h);
    global_dpd_->buf4_mat_irrep_rd(&DIjAb, h);
    for(ij=0; ij < RIjAb->params->rowtot[h]; ij++)
       for(ab=0; ab < RIjAb->params->coltot[h^C_irr]; ab++) {
           tval = eval - DIjAb.matrix[h][ij][ab];
           if (fabs(tval) > 0.0001) RIjAb->matrix[h][ij][ab] /= tval;
      }
    global_dpd_->buf4_mat_irrep_wrt(RIjAb, h);
    global_dpd_->buf4_mat_irrep_close(RIjAb, h);
    global_dpd_->buf4_mat_irrep_close(&DIjAb, h);
  }
  global_dpd_->buf4_close(&DIjAb);

  return;
}

void precondition_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, double eval)
{
  dpdfile2 DIA;
  dpdbuf4 DIjAb;
  int h, nirreps, i, j, a, b, ij, ab, ii, C_irr;
  double tval;

  /* Local correlation decs */
  int nso, nocc, nvir;
  double *T1tilde, *T1bar, **T2tilde, **T2bar, **X1, **X2;
  psio_address next;

  C_irr = RIA->my_irrep;
  nirreps = RIA->params->nirreps;

  if(params.local) {

    nso = local.nso;
    nocc = local.nocc;
    nvir = local.nvir;

    local.pairdom_len = init_int_array(nocc*nocc);
    local.pairdom_nrlen = init_int_array(nocc*nocc);
    local.eps_occ = init_array(nocc);
    local.weak_pairs = init_int_array(nocc*nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length", (char *) local.pairdom_len,
		    nocc*nocc*sizeof(int));
    psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length (Non-redundant basis)", (char *) local.pairdom_nrlen,
		    nocc*nocc*sizeof(int));
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) local.eps_occ,
		    nocc*sizeof(double));
    psio_read_entry(PSIF_CC_INFO, "Local Weak Pairs", (char *) local.weak_pairs,
		    nocc*nocc*sizeof(int));

    local.W = (double ***) malloc(nocc * nocc * sizeof(double **));
    local.V = (double ***) malloc(nocc * nocc * sizeof(double **));
    local.eps_vir = (double **) malloc(nocc * nocc * sizeof(double *));
    next = PSIO_ZERO;
    for(ij=0; ij < nocc*nocc; ij++) {
      local.eps_vir[ij] = init_array(local.pairdom_nrlen[ij]);
      psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *) local.eps_vir[ij],
		local.pairdom_nrlen[ij]*sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(ij=0; ij < nocc*nocc; ij++) {
      local.V[ij] = block_matrix(nvir,local.pairdom_len[ij]);
      psio_read(PSIF_CC_INFO, "Local Residual Vector (V)", (char *) local.V[ij][0],
		nvir*local.pairdom_len[ij]*sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(ij=0; ij < nocc*nocc; ij++) {
      local.W[ij] = block_matrix(local.pairdom_len[ij],local.pairdom_nrlen[ij]);
      psio_read(PSIF_CC_INFO, "Local Transformation Matrix (W)", (char *) local.W[ij][0],
		local.pairdom_len[ij]*local.pairdom_nrlen[ij]*sizeof(double), next, &next);
    }

    if(local.filter_singles) {

      global_dpd_->file2_mat_init(RIA);
      global_dpd_->file2_mat_rd(RIA);

      for(i=0; i < nocc; i++) {
	ii = i * nocc +i;

	if(!local.pairdom_len[ii]) {
	  outfile->Printf( "\n\tlocal_filter_T1: Pair ii = [%d] is zero-length, which makes no sense.\n",ii);
	  exit(2);
	}

	T1tilde = init_array(local.pairdom_len[ii]);
	T1bar = init_array(local.pairdom_nrlen[ii]);

	/* Transform the virtuals to the redundant projected virtual basis */
	C_DGEMV('t', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii],
		&(RIA->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

	/* Transform the virtuals to the non-redundant virtual basis */
	C_DGEMV('t', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]), local.pairdom_nrlen[ii],
		&(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

	for(a=0; a < local.pairdom_nrlen[ii]; a++) {
	  tval = eval + local.eps_occ[i] - local.eps_vir[ii][a];
	  if(fabs(tval) > 0.0001) T1bar[a] /= tval;
	}

	/* Transform the new T1's to the redundant projected virtual basis */
	C_DGEMV('n', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]), local.pairdom_nrlen[ii],
		&(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);

	/* Transform the new T1's to the MO basis */
	C_DGEMV('n', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii],
		&(T1tilde[0]), 1, 0.0, &(RIA->matrix[0][i][0]), 1);

	free(T1bar);
	free(T1tilde);
      }

      global_dpd_->file2_mat_wrt(RIA);
      global_dpd_->file2_mat_close(RIA);
    }
  }
  else {
    global_dpd_->file2_mat_init(RIA);
    global_dpd_->file2_mat_rd(RIA);
    global_dpd_->file2_init(&DIA, PSIF_EOM_D, C_irr, 0, 1, "DIA");
    global_dpd_->file2_mat_init(&DIA);
    global_dpd_->file2_mat_rd(&DIA);
    for(h=0; h < nirreps; h++)
      for(i=0; i < RIA->params->rowtot[h]; i++)
        for(a=0; a < RIA->params->coltot[h^C_irr]; a++) {
	  tval = eval - DIA.matrix[h][i][a];
	  if (fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
        }
    global_dpd_->file2_mat_wrt(RIA);
    global_dpd_->file2_mat_close(RIA);
    global_dpd_->file2_mat_close(&DIA);
    global_dpd_->file2_close(&DIA);

  }

  if(params.local) {

    global_dpd_->buf4_mat_irrep_init(RIjAb, 0);
    global_dpd_->buf4_mat_irrep_rd(RIjAb, 0);

    X1 = block_matrix(nso,nvir);
    X2 = block_matrix(nvir,nso);

    T2tilde = block_matrix(nso,nso);
    T2bar = block_matrix(nvir, nvir);
    for(i=0,ij=0; i < nocc; i++) {
      for(j=0; j < nocc; j++,ij++) {

	if(!local.weak_pairs[ij]) {

	  /* Transform the virtuals to the redundant projected virtual basis */
	  C_DGEMM('t', 'n', local.pairdom_len[ij], nvir, nvir, 1.0, &(local.V[ij][0][0]), local.pairdom_len[ij],
		  &(RIjAb->matrix[0][ij][0]), nvir, 0.0, &(X1[0][0]), nvir);
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
	      tval = eval + local.eps_occ[i]+local.eps_occ[j]-local.eps_vir[ij][a]-local.eps_vir[ij][b];
	      if(fabs(tval) > 0.0001) T2bar[a][b] /= tval;
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
		  &(local.V[ij][0][0]), local.pairdom_len[ij], 0.0, &(RIjAb->matrix[0][ij][0]), nvir);
	}
	else  /* This must be a neglected weak pair; force it to zero */
	  memset((void *) RIjAb->matrix[0][ij], 0, nvir*nvir*sizeof(double));
      }
    }
    free_block(T2tilde);
    free_block(T2bar);

    free_block(X1);
    free_block(X2);

    global_dpd_->buf4_mat_irrep_wrt(RIjAb, 0);
    global_dpd_->buf4_mat_irrep_close(RIjAb, 0);

    /* Free Local Memory */
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
  else {

    global_dpd_->buf4_init(&DIjAb, PSIF_EOM_D, C_irr, 0, 5, 0, 5, 0, "DIjAb");
    for(h=0; h < RIjAb->params->nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(RIjAb, h);
      global_dpd_->buf4_mat_irrep_init(&DIjAb, h);
      global_dpd_->buf4_mat_irrep_rd(RIjAb, h);
      global_dpd_->buf4_mat_irrep_rd(&DIjAb, h);
      for(ij=0; ij < RIjAb->params->rowtot[h]; ij++)
	for(ab=0; ab < RIjAb->params->coltot[h^C_irr]; ab++) {
	  tval = eval - DIjAb.matrix[h][ij][ab];
	  if (fabs(tval) > 0.0001) RIjAb->matrix[h][ij][ab] /= tval;
	}
      global_dpd_->buf4_mat_irrep_wrt(RIjAb, h);
      global_dpd_->buf4_mat_irrep_close(RIjAb, h);
      global_dpd_->buf4_mat_irrep_close(&DIjAb, h);
    }
    global_dpd_->buf4_close(&DIjAb);
  }

  return;
}

}} // namespace psi::cceom
