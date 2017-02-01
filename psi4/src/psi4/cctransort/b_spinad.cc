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

#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.hpp"
namespace psi { namespace cctransort {

void b_spinad(std::shared_ptr<PSIO> psio)
{
  dpdbuf4 B, Bs, Ba;

  // This should probably replace the (VV|VV) sort in sort_tei so that we limit disk usage
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  global_dpd_->buf4_init(&Bs, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
  global_dpd_->buf4_scm(&Bs, 0.0);
  global_dpd_->buf4_init(&Ba, PSIF_CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
  global_dpd_->buf4_scm(&Ba, 0.0);
  for(int h=0; h < B.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_row_init(&B, h);
    global_dpd_->buf4_mat_irrep_row_init(&Bs, h);
    global_dpd_->buf4_mat_irrep_row_init(&Ba, h);
    for(int ab=0; ab < Bs.params->rowtot[h]; ab++) {
      int a = Bs.params->roworb[h][ab][0];
      int b = Bs.params->roworb[h][ab][1];
      int AB = B.params->rowidx[a][b];
      global_dpd_->buf4_mat_irrep_row_rd(&B, h, AB);
      for(int cd=0; cd < Bs.params->coltot[h]; cd++ ) {
        int c = Bs.params->colorb[h][cd][0];
        int d = Bs.params->colorb[h][cd][1];
        int CD = B.params->colidx[c][d];
        int DC = B.params->colidx[d][c];
        Bs.matrix[h][0][cd] = B.matrix[h][0][CD] + B.matrix[h][0][DC];
        Ba.matrix[h][0][cd] = B.matrix[h][0][CD] - B.matrix[h][0][DC];
      }
      global_dpd_->buf4_mat_irrep_row_wrt(&Bs, h, ab);
      global_dpd_->buf4_mat_irrep_row_wrt(&Ba, h, ab);
    }
    global_dpd_->buf4_mat_irrep_row_close(&Ba, h);
    global_dpd_->buf4_mat_irrep_row_close(&Bs, h);
    global_dpd_->buf4_mat_irrep_row_close(&B, h);
  }
  global_dpd_->buf4_close(&Ba);
  global_dpd_->buf4_close(&Bs);
  global_dpd_->buf4_close(&B);

  /* Generate <ab|cc> components of B(+) */

  global_dpd_->buf4_init(&Bs, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
  int nvirt = 0;
  for(int h=0; h < Bs.params->nirreps; h++) nvirt += Bs.params->ppi[h];

  int rows_per_bucket = dpd_memfree()/(Bs.params->coltot[0] + nvirt);
  if(rows_per_bucket > Bs.params->rowtot[0]) rows_per_bucket = Bs.params->rowtot[0];
  int nbuckets = (int) ceil((double) Bs.params->rowtot[0]/(double) rows_per_bucket);
  int rows_left = Bs.params->rowtot[0] % rows_per_bucket;

  global_dpd_->buf4_mat_irrep_init_block(&Bs, 0, rows_per_bucket);
  double **B_diag = global_dpd_->dpd_block_matrix(rows_per_bucket, nvirt);

  psio_address next = PSIO_ZERO;
  int m;
  for(m=0; m < (rows_left ? nbuckets-1:nbuckets); m++) {
    int row_start = m * rows_per_bucket;
    global_dpd_->buf4_mat_irrep_rd_block(&Bs, 0, row_start, rows_per_bucket);
    for(int ab=0; ab < rows_per_bucket; ab++)
      for(int Gc=0; Gc < Bs.params->nirreps; Gc++)
        for(int C=0; C < Bs.params->rpi[Gc]; C++) {
          int c = Bs.params->roff[Gc] + C;
          int cc = Bs.params->colidx[c][c];
          B_diag[ab][c] = Bs.matrix[0][ab][cc];
        }
    psio->write(PSIF_CC_BINTS, "B(+) <ab|cc>", (char *) B_diag[0], rows_per_bucket*nvirt*sizeof(double), next, &next);
  }
  if(rows_left) {
    int row_start = m * rows_per_bucket;
    global_dpd_->buf4_mat_irrep_rd_block(&Bs, 0, row_start, rows_left);
    for(int ab=0; ab < rows_left; ab++)
      for(int Gc=0; Gc < Bs.params->nirreps; Gc++)
        for(int C=0; C < Bs.params->rpi[Gc]; C++) {
          int c = Bs.params->roff[Gc] + C;
          int cc = Bs.params->colidx[c][c];
          B_diag[ab][c] = Bs.matrix[0][ab][cc];
        }
    psio->write(PSIF_CC_BINTS, "B(+) <ab|cc>", (char *) B_diag[0], rows_left*nvirt*sizeof(double), next, &next);
  }
  global_dpd_->free_dpd_block(B_diag, rows_per_bucket, nvirt);
  global_dpd_->buf4_mat_irrep_close_block(&Bs, 0, rows_per_bucket);
  global_dpd_->buf4_close(&Bs);
}

}} // namespace psi::cctranssort
