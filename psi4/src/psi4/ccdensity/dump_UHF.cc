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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* DUMP_UHF(): Mulliken-order the UHF-CC two-electron density and dump
** it to a file for subsequent backtransformation.  Basically all we
** have to do is swap indices two and three, e.g.
**
** G'(pr,qs) = G(pq,rs)
**
** In order for the Mulliken-ordered density to be valid for the
** backtransformation algorithm used in TRANSQT, the final density
** must have eight-fold permutational symmetry like the original
** integrals.
*/

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define IOFF_MAX 32641

void dump_UHF(struct iwlbuf *AA, struct iwlbuf *BB, struct iwlbuf *AB, struct RHO_Params rho_params)
{
  int nirreps, nmo, h, row, col;
  int *qt_aocc, *qt_avir;
  int *qt_bocc, *qt_bvir;
  int *ioff;
  int p, q, r, s, P, Q, R, S, pr, qs;
  double value;
  dpdbuf4 G;

  /* initialize the ioff array */
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(p=1; p < IOFF_MAX; p++) ioff[p] = ioff[p-1] + p;

  qt_aocc = moinfo.qt_aocc;  qt_avir = moinfo.qt_avir;
  qt_bocc = moinfo.qt_bocc;  qt_bvir = moinfo.qt_bvir;
  nirreps = moinfo.nirreps;
  nmo = moinfo.nmo;

  psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
	/* psio_write_entry(PSIF_MO_OPDM, "MO-basis Alpha OPDM", (char *) moinfo.opdm_a[0],
		   sizeof(double)*nmo*nmo);
  psio_write_entry(PSIF_MO_OPDM, "MO-basis Beta OPDM", (char *) moinfo.opdm_b[0],
		   sizeof(double)*nmo*nmo); */
  psio_write_entry(PSIF_MO_OPDM, rho_params.opdm_a_lbl, (char *) moinfo.opdm_a[0],
		   sizeof(double)*nmo*nmo);
  psio_write_entry(PSIF_MO_OPDM, rho_params.opdm_b_lbl, (char *) moinfo.opdm_b[0],
		   sizeof(double)*nmo*nmo);
  psio_close(PSIF_MO_OPDM, 1);

if (!params.onepdm) {
  psio_open(PSIF_MO_LAG, PSIO_OPEN_OLD);
  psio_write_entry(PSIF_MO_LAG, "MO-basis Alpha Lagrangian", (char *) moinfo.I_a[0],
		   sizeof(double)*nmo*nmo);
  psio_write_entry(PSIF_MO_LAG, "MO-basis Beta Lagrangian", (char *) moinfo.I_b[0],
		   sizeof(double)*nmo*nmo);
  psio_close(PSIF_MO_LAG, 1);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 0, 0, "G(IJ,KL)");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "G(IJ,KL)");
  global_dpd_->buf4_dump(&G, AA, qt_aocc, qt_aocc, qt_aocc, qt_aocc, 1, 0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 12, 12, 0, "Gijkl");
  global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 10, 10, "G(ij,kl)");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "G(ij,kl)");
  global_dpd_->buf4_dump(&G, BB, qt_bocc, qt_bocc, qt_bocc, qt_bocc, 1, 0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_aocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_bocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_aocc[r];
	s = G.params->colorb[h][col][1]; S = qt_bocc[s];

	value = 2.0 * G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_aocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_aocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_aocc[r];
	s = G.params->colorb[h][col][1]; S = qt_avir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(AA, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_bocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_bocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_bocc[r];
	s = G.params->colorb[h][col][1]; S = qt_bvir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(BB, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_aocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_bocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_aocc[r];
	s = G.params->colorb[h][col][1]; S = qt_bvir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, P, R, Q, S, value, 0, "NULL", 0);
	iwl_buf_wrt_val(AB, P, R, S, Q, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_bocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_aocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_bocc[r];
	s = G.params->colorb[h][col][1]; S = qt_avir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, Q, S, P, R, value, 0, "NULL", 0);
	iwl_buf_wrt_val(AB, S, Q, P, R, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 20, 20, "G(IA,JB)");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "G(IA,JB)");
  global_dpd_->buf4_symm(&G);
  global_dpd_->buf4_dump(&G, AA, qt_aocc, qt_avir, qt_aocc, qt_avir, 1, 0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 15, 10, 15, 0, "G(ij,ab)");
  global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 30, 30, "G(ia,jb)");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "G(ia,jb)");
  global_dpd_->buf4_symm(&G);
  global_dpd_->buf4_dump(&G, BB, qt_bocc, qt_bvir, qt_bocc, qt_bvir, 1, 0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_aocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_bocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_avir[r];
	s = G.params->colorb[h][col][1]; S = qt_bvir[s];

	value = 2.0 * G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_aocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_avir[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_aocc[r];
	s = G.params->colorb[h][col][1]; S = qt_avir[s];

	value = 0.5 * G.matrix[h][row][col];

	iwl_buf_wrt_val(AA, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_bocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_bvir[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_bocc[r];
	s = G.params->colorb[h][col][1]; S = qt_bvir[s];

	value = 0.5 * G.matrix[h][row][col];

	iwl_buf_wrt_val(BB, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_aocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_bvir[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_aocc[r];
	s = G.params->colorb[h][col][1]; S = qt_bvir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_bocc[p];
      q = G.params->roworb[h][row][1]; Q = qt_avir[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_bocc[r];
	s = G.params->colorb[h][col][1]; S = qt_avir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, Q, S, P, R, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 5, 21, 7, 0, "GCIAB");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_avir[p];
      q = G.params->roworb[h][row][1]; Q = qt_aocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_avir[r];
	s = G.params->colorb[h][col][1]; S = qt_avir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(AA, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 15, 31, 17, 0, "Gciab");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_bvir[p];
      q = G.params->roworb[h][row][1]; Q = qt_bocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_bvir[r];
	s = G.params->colorb[h][col][1]; S = qt_bvir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(BB, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_avir[p];
      q = G.params->roworb[h][row][1]; Q = qt_bocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_avir[r];
	s = G.params->colorb[h][col][1]; S = qt_bvir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, P, R, Q, S, value, 0, "NULL", 0);
	iwl_buf_wrt_val(AB, P, R, S, Q, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_bvir[p];
      q = G.params->roworb[h][row][1]; Q = qt_aocc[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_bvir[r];
	s = G.params->colorb[h][col][1]; S = qt_avir[s];

	value = G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, Q, S, P, R, value, 0, "NULL", 0);
	iwl_buf_wrt_val(AB, S, Q, P, R, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 7, 7, 0, "GABCD");
  global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 5, 5, "G(AB,CD)");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 5, 5, 5, 5, 0, "G(AB,CD)");
  global_dpd_->buf4_symm(&G);
  global_dpd_->buf4_dump(&G, AA, qt_avir, qt_avir, qt_avir, qt_avir, 1, 0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 15, 15, 17, 17, 0, "Gabcd");
  global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 15, 15, "G(ab,cd)");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 15, 15, 15, 15, 0, "G(ab,cd)");
  global_dpd_->buf4_symm(&G);
  global_dpd_->buf4_dump(&G, BB, qt_bvir, qt_bvir, qt_bvir, qt_bvir, 1, 0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
  for(h=0; h < G.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0]; P = qt_avir[p];
      q = G.params->roworb[h][row][1]; Q = qt_bvir[q];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0]; R = qt_avir[r];
	s = G.params->colorb[h][col][1]; S = qt_bvir[s];

	value = 2.0 * G.matrix[h][row][col];

	iwl_buf_wrt_val(AB, P, R, Q, S, value, 0, "NULL", 0);
      }
    }
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);
}

  free(ioff);
}

}} // namespace psi::ccdensity
