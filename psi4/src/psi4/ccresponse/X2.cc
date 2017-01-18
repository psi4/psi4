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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void denom2(dpdbuf4 *X2, double omega);
void local_filter_T2(dpdbuf4 *T2);

void X2_build(const char *pert, int irrep, double omega)
{
  dpdfile2 X1, z, F, t1;
  dpdbuf4 X2, X2new, Z, Z1, Z2, W, T2, I;
  char lbl[32];
  int Gej, Gab, Gij, Ge, Gj, Gi, nrows, length, E, e, II;
  int Gbm, Gfe, bm, b, m, Gb, Gm, Gf, B, M, fe, f, ef, ncols;
  double *X;
  dpdbuf4 S, A, B_s;
  int ij, Gc, C, c, cc;
  int rows_per_bucket, nbuckets, row_start, rows_left, nlinks;
  psio_address next;
  double **X_diag, **B_diag;

  sprintf(lbl, "%sBAR_IjAb", pert);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_copy(&X2new, PSIF_CC_LR, lbl);
  global_dpd_->buf4_close(&X2new);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  /*** D-S ***/

  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  global_dpd_->contract244(&X1, &W, &Z, 0, 0, 1, 1, 0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_axpy(&Z, &X2new, -1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  /*  dpd_contract244(&X1, &W, &Z, 1, 0, 0, 1, 0); */
  /* ooc code below added 7/28/05, -TDC */
  global_dpd_->file2_mat_init(&X1);
  global_dpd_->file2_mat_rd(&X1);
  for(Gej=0; Gej < moinfo.nirreps; Gej++) {
    Gab = Gej; /* W is totally symmetric */
    Gij = Gab ^ irrep;
    global_dpd_->buf4_mat_irrep_init(&Z, Gij);
    global_dpd_->buf4_mat_irrep_shift13(&Z, Gij);
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gj = Ge ^ Gej;
      Gi = Gj ^ Gij;
      nrows = moinfo.occpi[Gj];
      length = nrows * W.params->coltot[Gab];
      global_dpd_->buf4_mat_irrep_init_block(&W, Gej, nrows);
      for(E=0; E < moinfo.virtpi[Ge]; E++) {
	e = moinfo.vir_off[Ge] + E;
	global_dpd_->buf4_mat_irrep_rd_block(&W, Gej, W.row_offset[Gej][e], nrows);
	for(II=0; II < moinfo.occpi[Gi]; II++) {
	  if(length)
	    C_DAXPY(length, X1.matrix[Gi][II][E], W.matrix[Gej][0], 1, Z.shift.matrix[Gij][Gi][II], 1);
	}
      }
      global_dpd_->buf4_mat_irrep_close_block(&W, Gej, nrows);
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, Gij);
    global_dpd_->buf4_mat_irrep_close(&Z, Gij);
  }
  global_dpd_->file2_mat_close(&X1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_axpy(&Z, &X2new, 1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, 1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z);

  sprintf(lbl, "z(N,I) %s", pert);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep, 0, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
  global_dpd_->dot23(&X1, &W, &z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&W);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->contract244(&z, &T2, &Z, 0, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->file2_close(&z);
  global_dpd_->buf4_axpy(&Z, &X2new, -1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z);

  sprintf(lbl, "z(A,E) %s", pert);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep, 1, 1, lbl);
  /*   dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); */
  /*  dpd_dot24(&X1, &W, &z, 0, 0, 1, 0); */
  /* ooc code below added 7/28/05, -TDC */
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->file2_scm(&z, 0);
  global_dpd_->file2_mat_init(&z);
  global_dpd_->file2_mat_init(&X1);
  global_dpd_->file2_mat_rd(&X1);
  for(Gbm=0; Gbm < moinfo.nirreps; Gbm++) {
    Gfe = Gbm; /* W is totally symmetric */
    global_dpd_->buf4_mat_irrep_row_init(&W, Gbm);
    X = init_array(W.params->coltot[Gfe]);
    for(bm=0; bm < W.params->rowtot[Gbm]; bm++) {
      global_dpd_->buf4_mat_irrep_row_rd(&W, Gbm, bm);
      b = W.params->roworb[Gbm][bm][0];
      m = W.params->roworb[Gbm][bm][1];
      Gb = W.params->psym[b];
      Gm = Gbm ^ Gb;
      Ge = Gm ^ irrep;
      Gf = Gfe ^ Ge;
      B = b - moinfo.vir_off[Gb];
      M = m - moinfo.occ_off[Gm];
      zero_arr(X, W.params->coltot[Gfe]);
      for(fe=0; fe < W.params->coltot[Gfe]; fe++) {
	f = W.params->colorb[Gfe][fe][0];
	e = W.params->colorb[Gfe][fe][1];
	ef = W.params->colidx[e][f];
	X[fe] = 2.0 * W.matrix[Gbm][0][fe] - W.matrix[Gbm][0][ef];
      }
      nrows = moinfo.virtpi[Gf];
      ncols = moinfo.virtpi[Ge];
      if(nrows & ncols)
	C_DGEMV('n',nrows,ncols,1,&X[W.col_offset[Gfe][Gf]],ncols,
		X1.matrix[Gm][M],1,1,z.matrix[Gb][B],1);
    }
    free(X);
    global_dpd_->buf4_mat_irrep_row_close(&W, Gbm);
  }
  global_dpd_->file2_mat_close(&X1);
  global_dpd_->file2_mat_wrt(&z);
  global_dpd_->file2_mat_close(&z);
  /* end ooc additions, 7/28/05, -TDC */
  global_dpd_->buf4_close(&W);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->contract424(&T2, &z, &Z, 3, 1, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->file2_close(&z);
  global_dpd_->buf4_axpy(&Z, &X2new, 1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, 1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z);

  global_dpd_->file2_close(&X1);

  /*** D-D ***/

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  global_dpd_->buf4_axpy(&X2, &X2new, -omega);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "FAE");
  global_dpd_->contract424(&X2, &F, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&F);
  global_dpd_->buf4_axpy(&Z, &X2new, 1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, 1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "FMI");
  global_dpd_->contract244(&F, &X2, &Z, 0, 0, 0, 1, 0);
  global_dpd_->file2_close(&F);
  global_dpd_->buf4_axpy(&Z, &X2new, -1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
  global_dpd_->contract444(&W, &X2, &X2new, 1, 1, 1, 1);
  global_dpd_->buf4_close(&W);

  if(params.abcd == "OLD") {
    sprintf(lbl, "Z(Ab,Ij) %s", pert);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 5, 0, 5, 0, 0, lbl);
    global_dpd_->buf4_init(&I, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract444(&I, &X2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
    global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
    global_dpd_->buf4_close(&Z);
  }
  else if(params.abcd == "NEW") {
    timer_on("ABCD:new");

    global_dpd_->buf4_close(&X2);

    timer_on("ABCD:S");
    sprintf(lbl, "X_%s_(+)(ij,ab) (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 3, 8, 3, 8, 0, lbl);
    global_dpd_->buf4_init(&I, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
    sprintf(lbl, "S_%s_(ab,ij)", pert);
    global_dpd_->buf4_init(&S, PSIF_CC_TMP0, irrep, 8, 3, 8, 3, 0, lbl);
    global_dpd_->contract444(&I, &X2, &S, 0, 0, 0.5, 0);
    global_dpd_->buf4_close(&S);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&X2);
    timer_off("ABCD:S");

    /* X_diag(ij,c)  = 2 * X(ij,cc)*/
    /* NB: Gcc = 0 and B is totally symmetry, so Gab = 0 */
    /* But Gij = irrep ^ Gab = irrep */
    sprintf(lbl, "X_%s_(+)(ij,ab) (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 3, 8, 3, 8, 0, lbl);
    global_dpd_->buf4_mat_irrep_init(&X2, irrep);
    global_dpd_->buf4_mat_irrep_rd(&X2, irrep);
    X_diag = global_dpd_->dpd_block_matrix(X2.params->rowtot[irrep], moinfo.nvirt);
    for(ij=0; ij < X2.params->rowtot[irrep]; ij++)
      for(Gc=0; Gc < moinfo.nirreps; Gc++)
	for(C=0; C < moinfo.virtpi[Gc]; C++) {
	  c = C + moinfo.vir_off[Gc];
	  cc = X2.params->colidx[c][c];
	  X_diag[ij][c] = X2.matrix[irrep][ij][cc];
	}
    global_dpd_->buf4_mat_irrep_close(&X2, irrep);

    global_dpd_->buf4_init(&B_s, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
    sprintf(lbl, "S_%s_(ab,ij)", pert);
    global_dpd_->buf4_init(&S, PSIF_CC_TMP0, irrep, 8, 3, 8, 3, 0, lbl);
    global_dpd_->buf4_mat_irrep_init(&S, 0);
    global_dpd_->buf4_mat_irrep_rd(&S, 0);

    rows_per_bucket = dpd_memfree()/(B_s.params->coltot[0] + moinfo.nvirt);
    if(rows_per_bucket > B_s.params->rowtot[0]) rows_per_bucket = B_s.params->rowtot[0];
    nbuckets = (int) ceil((double) B_s.params->rowtot[0]/(double) rows_per_bucket);
    rows_left = B_s.params->rowtot[0] % rows_per_bucket;

    B_diag = global_dpd_->dpd_block_matrix(rows_per_bucket, moinfo.nvirt);
    next = PSIO_ZERO;
    ncols = X2.params->rowtot[irrep];
    nlinks = moinfo.nvirt;
    for(m=0; m < (rows_left ? nbuckets-1:nbuckets); m++) {
      row_start = m * rows_per_bucket;
      nrows = rows_per_bucket;
      if(nrows && ncols && nlinks) {
	psio_read(PSIF_CC_BINTS,"B(+) <ab|cc>",(char *) B_diag[0],nrows*nlinks*sizeof(double),next, &next);
	C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks,
		X_diag[0], nlinks, 1, S.matrix[0][row_start], ncols);
      }

    }
    if(rows_left) {
      row_start = m * rows_per_bucket;
      nrows = rows_left;
      if(nrows && ncols && nlinks) {
	psio_read(PSIF_CC_BINTS,"B(+) <ab|cc>",(char *) B_diag[0],nrows*nlinks*sizeof(double),next, &next);
	C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks,
		X_diag[0], nlinks, 1, S.matrix[0][row_start], ncols);
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&S, 0);
    global_dpd_->buf4_mat_irrep_close(&S, 0);
    global_dpd_->buf4_close(&S);
    global_dpd_->buf4_close(&B_s);
    global_dpd_->free_dpd_block(B_diag, rows_per_bucket, moinfo.nvirt);
    global_dpd_->free_dpd_block(X_diag, X2.params->rowtot[irrep], moinfo.nvirt);
    global_dpd_->buf4_close(&X2);

    timer_on("ABCD:A");
    sprintf(lbl, "X_%s_(-)(ij,ab) (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 4, 9, 4, 9, 0, lbl);
    global_dpd_->buf4_init(&I, PSIF_CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
    sprintf(lbl, "A_%s_(ab,ij)", pert);
    global_dpd_->buf4_init(&A, PSIF_CC_TMP0, irrep, 9, 4, 9, 4, 0, lbl);
    global_dpd_->contract444(&I, &X2, &A, 0, 0, 0.5, 0);
    global_dpd_->buf4_close(&A);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&X2);
    timer_off("ABCD:A");

    timer_on("ABCD:axpy");
    global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
    sprintf(lbl, "S_%s_(ab,ij)", pert);
    global_dpd_->buf4_init(&S, PSIF_CC_TMP0, irrep, 5, 0, 8, 3, 0, lbl);
    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_sort_axpy(&S, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
    global_dpd_->buf4_close(&S);
    sprintf(lbl, "A_%s_(ab,ij)", pert);
    global_dpd_->buf4_init(&A, PSIF_CC_TMP0, irrep, 5, 0, 9, 4, 0, lbl);
    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_sort_axpy(&A, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
    global_dpd_->buf4_close(&A);
    global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */

    timer_off("ABCD:axpy");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    timer_off("ABCD:new");
  }

  sprintf(lbl, "Z(Mb,Ij) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 10, 0, 10, 0, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->contract444(&I, &X2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&I);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&t1, &Z, &Z1, 0, 0, 1, 1, 0);
  global_dpd_->file2_close(&t1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_axpy(&Z1, &X2new, -1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z1);

  sprintf(lbl, "Z(Ij,Mn) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 0, 0, 0, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract444(&X2, &I, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  global_dpd_->contract444(&Z, &T2, &X2new, 0, 1, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_close(&X2);

  sprintf(lbl, "Z(Ib,jA) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 10, 10, 10, 10, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
  global_dpd_->contract444(&X2, &W, &Z, 0, 1, 1, 0);
  global_dpd_->buf4_close(&X2);
  global_dpd_->buf4_close(&W);
  sprintf(lbl, "X(IA,jb) III %s", pert);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 10, 10, lbl);
  global_dpd_->buf4_close(&Z);
  sprintf(lbl, "X(IA,jb) I %s", pert);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%s_(2IAjb-IbjA) (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
  global_dpd_->contract444(&X2, &W, &Z1, 0, 1, 0.5, 0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&X2);
  sprintf(lbl, "Z(Ib,jA) %s", pert);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep, 10, 10, 10, 10, 0, lbl);
  global_dpd_->buf4_axpy(&Z2, &Z1, 0.5);
  global_dpd_->buf4_close(&Z2);
  sprintf(lbl, "X(IA,jb) III %s", pert);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep, 10, 10, 10, 10, 0, lbl);
  global_dpd_->buf4_axpy(&Z2, &Z1, 1);
  global_dpd_->buf4_close(&Z2);
  sprintf(lbl, "X(Ij,Ab) I+III %s", pert);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 0, 5, lbl);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_axpy(&Z1, &X2new, 1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_LR, qpsr, 0, 5, lbl, 1); /* II+IV */
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z1);

  sprintf(lbl, "z(F,A) %s", pert);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep, 1, 1, lbl);
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract442(&I, &X2, &z, 2, 2, 1, 0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&X2);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->contract424(&T2, &z, &Z, 3, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->file2_close(&z);
  global_dpd_->buf4_axpy(&Z, &X2new, -1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z);

  sprintf(lbl, "z(N,I) %s", pert);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep, 0, 0, lbl);
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract442(&I, &X2, &z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&X2);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->contract244(&z, &T2, &Z, 0, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->file2_close(&z);
  global_dpd_->buf4_axpy(&Z, &X2new, -1);
  global_dpd_->buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  global_dpd_->buf4_close(&Z);

  if(params.local) local_filter_T2(&X2new);
  else denom2(&X2new, omega);
  global_dpd_->buf4_close(&X2new);
}

}} // namespace psi::ccresponse
