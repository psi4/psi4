/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
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
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_contract244(&X1, &W, &Z, 0, 0, 1, 1, 0);
  dpd_buf4_close(&W);
  dpd_buf4_axpy(&Z, &X2new, -1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, -1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  /*  dpd_contract244(&X1, &W, &Z, 1, 0, 0, 1, 0); */
  /* ooc code below added 7/28/05, -TDC */
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);
  for(Gej=0; Gej < moinfo.nirreps; Gej++) {
    Gab = Gej; /* W is totally symmetric */
    Gij = Gab ^ irrep;
    dpd_buf4_mat_irrep_init(&Z, Gij);
    dpd_buf4_mat_irrep_shift13(&Z, Gij);
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gj = Ge ^ Gej;
      Gi = Gj ^ Gij;
      nrows = moinfo.occpi[Gj];
      length = nrows * W.params->coltot[Gab];
      dpd_buf4_mat_irrep_init_block(&W, Gej, nrows);
      for(E=0; E < moinfo.virtpi[Ge]; E++) {
	e = moinfo.vir_off[Ge] + E;
	dpd_buf4_mat_irrep_rd_block(&W, Gej, W.row_offset[Gej][e], nrows);
	for(II=0; II < moinfo.occpi[Gi]; II++) {
	  if(length)
	    C_DAXPY(length, X1.matrix[Gi][II][E], W.matrix[Gej][0], 1, Z.shift.matrix[Gij][Gi][II], 1);
	}
      }
      dpd_buf4_mat_irrep_close_block(&W, Gej, nrows);
    }
    dpd_buf4_mat_irrep_wrt(&Z, Gij);
    dpd_buf4_mat_irrep_close(&Z, Gij);
  }
  dpd_file2_mat_close(&X1);
  dpd_buf4_close(&W);
  dpd_buf4_axpy(&Z, &X2new, 1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, 1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z);

  sprintf(lbl, "z(N,I) %s", pert);
  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
  dpd_dot23(&X1, &W, &z, 0, 0, 1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&z, &T2, &Z, 0, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&z);
  dpd_buf4_axpy(&Z, &X2new, -1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, -1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z);

  sprintf(lbl, "z(A,E) %s", pert);
  dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, lbl);
  /*   dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); */
  /*  dpd_dot24(&X1, &W, &z, 0, 0, 1, 0); */
  /* ooc code below added 7/28/05, -TDC */
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_file2_scm(&z, 0);
  dpd_file2_mat_init(&z);
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);
  for(Gbm=0; Gbm < moinfo.nirreps; Gbm++) {
    Gfe = Gbm; /* W is totally symmetric */
    dpd_buf4_mat_irrep_row_init(&W, Gbm);
    X = init_array(W.params->coltot[Gfe]);
    for(bm=0; bm < W.params->rowtot[Gbm]; bm++) {
      dpd_buf4_mat_irrep_row_rd(&W, Gbm, bm);
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
    dpd_buf4_mat_irrep_row_close(&W, Gbm);
  }
  dpd_file2_mat_close(&X1);
  dpd_file2_mat_wrt(&z);
  dpd_file2_mat_close(&z);
  /* end ooc additions, 7/28/05, -TDC */
  dpd_buf4_close(&W);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&T2, &z, &Z, 3, 1, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&z);
  dpd_buf4_axpy(&Z, &X2new, 1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, 1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z);

  dpd_file2_close(&X1);

  /*** D-D ***/

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  dpd_buf4_axpy(&X2, &X2new, -omega);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "FAE");
  dpd_contract424(&X2, &F, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_buf4_axpy(&Z, &X2new, 1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, 1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "FMI");
  dpd_contract244(&F, &X2, &Z, 0, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_buf4_axpy(&Z, &X2new, -1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, -1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z);

  dpd_buf4_init(&W, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
  dpd_contract444(&W, &X2, &X2new, 1, 1, 1, 1);
  dpd_buf4_close(&W);

  if(params.abcd == "OLD") {
    sprintf(lbl, "Z(Ab,Ij) %s", pert);
    dpd_buf4_init(&Z, CC_TMP0, irrep, 5, 0, 5, 0, 0, lbl);
    dpd_buf4_init(&I, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract444(&I, &X2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    dpd_buf4_sort_axpy(&Z, CC_LR, rspq, 0, 5, lbl, 1);
    dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
    dpd_buf4_close(&Z);
  }
  else if(params.abcd == "NEW") {
    timer_on("ABCD:new");

    dpd_buf4_close(&X2);

    timer_on("ABCD:S");
    sprintf(lbl, "X_%s_(+)(ij,ab) (%5.3f)", pert, omega);
    dpd_buf4_init(&X2, CC_LR, irrep, 3, 8, 3, 8, 0, lbl);
    dpd_buf4_init(&I, CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
    sprintf(lbl, "S_%s_(ab,ij)", pert);
    dpd_buf4_init(&S, CC_TMP0, irrep, 8, 3, 8, 3, 0, lbl);
    dpd_contract444(&I, &X2, &S, 0, 0, 0.5, 0);
    dpd_buf4_close(&S);
    dpd_buf4_close(&I);
    dpd_buf4_close(&X2);
    timer_off("ABCD:S");

    /* X_diag(ij,c)  = 2 * X(ij,cc)*/
    /* NB: Gcc = 0 and B is totally symmetry, so Gab = 0 */
    /* But Gij = irrep ^ Gab = irrep */
    sprintf(lbl, "X_%s_(+)(ij,ab) (%5.3f)", pert, omega);
    dpd_buf4_init(&X2, CC_LR, irrep, 3, 8, 3, 8, 0, lbl);
    dpd_buf4_mat_irrep_init(&X2, irrep);
    dpd_buf4_mat_irrep_rd(&X2, irrep);
    X_diag = dpd_block_matrix(X2.params->rowtot[irrep], moinfo.nvirt);
    for(ij=0; ij < X2.params->rowtot[irrep]; ij++)
      for(Gc=0; Gc < moinfo.nirreps; Gc++)
	for(C=0; C < moinfo.virtpi[Gc]; C++) {
	  c = C + moinfo.vir_off[Gc];
	  cc = X2.params->colidx[c][c];
	  X_diag[ij][c] = X2.matrix[irrep][ij][cc];
	}
    dpd_buf4_mat_irrep_close(&X2, irrep);

    dpd_buf4_init(&B_s, CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
    sprintf(lbl, "S_%s_(ab,ij)", pert);
    dpd_buf4_init(&S, CC_TMP0, irrep, 8, 3, 8, 3, 0, lbl);
    dpd_buf4_mat_irrep_init(&S, 0);
    dpd_buf4_mat_irrep_rd(&S, 0);

    rows_per_bucket = dpd_memfree()/(B_s.params->coltot[0] + moinfo.nvirt);
    if(rows_per_bucket > B_s.params->rowtot[0]) rows_per_bucket = B_s.params->rowtot[0];
    nbuckets = (int) ceil((double) B_s.params->rowtot[0]/(double) rows_per_bucket);
    rows_left = B_s.params->rowtot[0] % rows_per_bucket;

    B_diag = dpd_block_matrix(rows_per_bucket, moinfo.nvirt);
    next = PSIO_ZERO;
    ncols = X2.params->rowtot[irrep];
    nlinks = moinfo.nvirt;
    for(m=0; m < (rows_left ? nbuckets-1:nbuckets); m++) {
      row_start = m * rows_per_bucket;
      nrows = rows_per_bucket;
      if(nrows && ncols && nlinks) {
	psio_read(CC_BINTS,"B(+) <ab|cc>",(char *) B_diag[0],nrows*nlinks*sizeof(double),next, &next);
	C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks,
		X_diag[0], nlinks, 1, S.matrix[0][row_start], ncols);
      }

    }
    if(rows_left) {
      row_start = m * rows_per_bucket;
      nrows = rows_left;
      if(nrows && ncols && nlinks) {
	psio_read(CC_BINTS,"B(+) <ab|cc>",(char *) B_diag[0],nrows*nlinks*sizeof(double),next, &next);
	C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks,
		X_diag[0], nlinks, 1, S.matrix[0][row_start], ncols);
      }
    }
    dpd_buf4_mat_irrep_wrt(&S, 0);
    dpd_buf4_mat_irrep_close(&S, 0);
    dpd_buf4_close(&S);
    dpd_buf4_close(&B_s);
    dpd_free_block(B_diag, rows_per_bucket, moinfo.nvirt);
    dpd_free_block(X_diag, X2.params->rowtot[irrep], moinfo.nvirt);
    dpd_buf4_close(&X2);

    timer_on("ABCD:A");
    sprintf(lbl, "X_%s_(-)(ij,ab) (%5.3f)", pert, omega);
    dpd_buf4_init(&X2, CC_LR, irrep, 4, 9, 4, 9, 0, lbl);
    dpd_buf4_init(&I, CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
    sprintf(lbl, "A_%s_(ab,ij)", pert);
    dpd_buf4_init(&A, CC_TMP0, irrep, 9, 4, 9, 4, 0, lbl);
    dpd_contract444(&I, &X2, &A, 0, 0, 0.5, 0);
    dpd_buf4_close(&A);
    dpd_buf4_close(&I);
    dpd_buf4_close(&X2);
    timer_off("ABCD:A");

    timer_on("ABCD:axpy");
    dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
    sprintf(lbl, "S_%s_(ab,ij)", pert);
    dpd_buf4_init(&S, CC_TMP0, irrep, 5, 0, 8, 3, 0, lbl);
    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    dpd_buf4_sort_axpy(&S, CC_LR, rspq, 0, 5, lbl, 1);
    dpd_buf4_close(&S);
    sprintf(lbl, "A_%s_(ab,ij)", pert);
    dpd_buf4_init(&A, CC_TMP0, irrep, 5, 0, 9, 4, 0, lbl);
    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    dpd_buf4_sort_axpy(&A, CC_LR, rspq, 0, 5, lbl, 1);
    dpd_buf4_close(&A);
    dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */

    timer_off("ABCD:axpy");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    timer_off("ABCD:new");
  }

  sprintf(lbl, "Z(Mb,Ij) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 10, 0, 10, 0, 0, lbl);
  dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_contract444(&I, &X2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z1, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&t1, &Z, &Z1, 0, 0, 1, 1, 0);
  dpd_file2_close(&t1);
  dpd_buf4_close(&Z);
  dpd_buf4_axpy(&Z1, &X2new, -1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z1, CC_LR, qpsr, 0, 5, lbl, -1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z1);

  sprintf(lbl, "Z(Ij,Mn) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 0, 0, 0, 0, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract444(&X2, &I, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_contract444(&Z, &T2, &X2new, 0, 1, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&Z);

  dpd_buf4_close(&X2);

  sprintf(lbl, "Z(Ib,jA) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 10, 10, 10, 10, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
  dpd_contract444(&X2, &W, &Z, 0, 1, 1, 0);
  dpd_buf4_close(&X2);
  dpd_buf4_close(&W);
  sprintf(lbl, "X(IA,jb) III %s", pert);
  dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 10, lbl);
  dpd_buf4_close(&Z);
  sprintf(lbl, "X(IA,jb) I %s", pert);
  dpd_buf4_init(&Z1, CC_TMP0, irrep, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%s_(2IAjb-IbjA) (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
  dpd_contract444(&X2, &W, &Z1, 0, 1, 0.5, 0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&X2);
  sprintf(lbl, "Z(Ib,jA) %s", pert);
  dpd_buf4_init(&Z2, CC_TMP0, irrep, 10, 10, 10, 10, 0, lbl);
  dpd_buf4_axpy(&Z2, &Z1, 0.5);
  dpd_buf4_close(&Z2);
  sprintf(lbl, "X(IA,jb) III %s", pert);
  dpd_buf4_init(&Z2, CC_TMP0, irrep, 10, 10, 10, 10, 0, lbl);
  dpd_buf4_axpy(&Z2, &Z1, 1);
  dpd_buf4_close(&Z2);
  sprintf(lbl, "X(Ij,Ab) I+III %s", pert);
  dpd_buf4_sort(&Z1, CC_TMP0, prqs, 0, 5, lbl);
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_axpy(&Z1, &X2new, 1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z1, CC_LR, qpsr, 0, 5, lbl, 1); /* II+IV */
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z1);

  sprintf(lbl, "z(F,A) %s", pert);
  dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, lbl);
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&I, &X2, &z, 2, 2, 1, 0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&X2);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&T2, &z, &Z, 3, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&z);
  dpd_buf4_axpy(&Z, &X2new, -1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, -1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z);

  sprintf(lbl, "z(N,I) %s", pert);
  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, lbl);
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&I, &X2, &z, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&X2);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  dpd_buf4_init(&Z, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&z, &T2, &Z, 0, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&z);
  dpd_buf4_axpy(&Z, &X2new, -1);
  dpd_buf4_close(&X2new);  /* Need to close X2new to avoid collisions */
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_sort_axpy(&Z, CC_LR, qpsr, 0, 5, lbl, -1);
  dpd_buf4_init(&X2new, CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
  dpd_buf4_close(&Z);

  if(params.local) local_filter_T2(&X2new);
  else denom2(&X2new, omega);
  dpd_buf4_close(&X2new);
}

}} // namespace psi::ccresponse
