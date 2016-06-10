/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi {
  namespace transqt2 {

void transtwo_rhf(void)
{
  int nirreps, nso, nmo;
  double ***C, **scratch, **TMP;
  int h, p, q, r, s, pq, rs, Gs, Gr, PQ, RS;
  int nrows, ncols, nlinks;
  unsigned long int memfree, rows_per_bucket, rows_left, this_bucket_rows;
  int nbuckets, n;
  dpdbuf4 J, K;
  struct iwlbuf MBuff;
  int stat;
  int *C_offset;

  nirreps = moinfo.nirreps;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  C = moinfo.C;
  C_offset = moinfo.C_offset;

  TMP = block_matrix(nso,nso);

  if(params.print_lvl) {
    outfile->Printf( "\tStarting first half-transformation.\n");
    
  }

  psio_open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
  psio_open(PSIF_HALFT0, PSIO_OPEN_NEW);

  timer_on("RHF:1sthalf");
  global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (pq,rs)");
  global_dpd_->buf4_init(&K, PSIF_HALFT0, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  for(h=0; h < nirreps; h++) {

    memfree = (unsigned long int) (dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
    if(J.params->coltot[h]) {
      rows_per_bucket = memfree/(2 * J.params->coltot[h]);
      if(rows_per_bucket > J.params->rowtot[h]) rows_per_bucket = (unsigned long int) J.params->rowtot[h];
      nbuckets = (int) ceil(((double) J.params->rowtot[h])/((double) rows_per_bucket));
      rows_left = (unsigned long int) (J.params->rowtot[h] % rows_per_bucket);
    }
    else {
      nbuckets = 0;
      rows_per_bucket = 0;
      rows_left = 0;
    }

    if(params.print_lvl > 1) {
      outfile->Printf( "\th = %d; memfree         = %lu\n", h, memfree);
      outfile->Printf( "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      outfile->Printf( "\th = %d; rows_left       = %lu\n", h, rows_left);
      outfile->Printf( "\th = %d; nbuckets        = %d\n", h, nbuckets);
      
    }

    global_dpd_->buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    global_dpd_->buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < nbuckets; n++) {
      if(nbuckets == 1) this_bucket_rows = rows_per_bucket;
      else this_bucket_rows = (n < nbuckets-1) ? rows_per_bucket : rows_left;
      global_dpd_->buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, this_bucket_rows);
      for(pq=0; pq < this_bucket_rows; pq++) {
        for(Gr=0; Gr < nirreps; Gr++) {
          Gs = h^Gr;
          nrows = moinfo.sopi[Gr];
          ncols = moinfo.actpi[Gs];
          nlinks = moinfo.sopi[Gs];
          rs = J.col_offset[h][Gr];
          if(nrows && ncols && nlinks)
            C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][pq][rs],nlinks,
                    &(C[Gs][0][C_offset[Gs]]),moinfo.mopi[Gs],0.0,TMP[0],nso);

          nrows = moinfo.actpi[Gr];
          ncols = moinfo.actpi[Gs];
          nlinks = moinfo.sopi[Gr];
          rs = K.col_offset[h][Gr];
          if(nrows && ncols && nlinks)
            C_DGEMM('t','n',nrows,ncols,nlinks,1.0,&(C[Gr][0][C_offset[Gr]]),moinfo.mopi[Gr],TMP[0],nso,
                    0.0,&K.matrix[h][pq][rs],ncols);
        } /* Gr */
      } /* pq */
      global_dpd_->buf4_mat_irrep_wrt_block(&K, h, n*rows_per_bucket, this_bucket_rows);
    }

    global_dpd_->buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    global_dpd_->buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  global_dpd_->buf4_close(&K);
  global_dpd_->buf4_close(&J);

  psio_close(PSIF_SO_PRESORT, 0);

  timer_off("RHF:1sthalf");
  timer_on("RHF:midsort");
  if(params.print_lvl) {
    outfile->Printf( "\tSorting half-transformed integrals.\n");
    
  }

  psio_open(PSIF_HALFT1, PSIO_OPEN_NEW);

  global_dpd_->buf4_init(&K, PSIF_HALFT0, 0, 3, 8, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  global_dpd_->buf4_sort(&K, PSIF_HALFT1, rspq, 8, 3, "Half-Transformed Ints (ij,pq)");
  global_dpd_->buf4_close(&K);

  psio_close(PSIF_HALFT0, 0);
  timer_off("RHF:midsort");
  timer_on("RHF:2ndhalf");

  if(params.print_lvl) {
    outfile->Printf( "\tStarting second half-transformation.\n");
    
  }
  iwl_buf_init(&MBuff, PSIF_MO_TEI, params.tolerance, 0, 0);

  global_dpd_->buf4_init(&J, PSIF_HALFT1, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  global_dpd_->buf4_init(&K, PSIF_CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)");
  for(h=0; h < nirreps; h++) {

    memfree = (unsigned long int) (dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
    if(J.params->coltot[h]) {
      rows_per_bucket = memfree/(2 * J.params->coltot[h]);
      if(rows_per_bucket > J.params->rowtot[h]) rows_per_bucket = (unsigned long int) J.params->rowtot[h];
      nbuckets = (int) ceil(((double) J.params->rowtot[h])/((double) rows_per_bucket));
      rows_left = (unsigned long int) (J.params->rowtot[h] % rows_per_bucket);
    }
    else {
      nbuckets = 0;
      rows_per_bucket = 0;
      rows_left = 0;
    }

    if(params.print_lvl > 1) {
      outfile->Printf( "\th = %d; memfree         = %lu\n", h, memfree);
      outfile->Printf( "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      outfile->Printf( "\th = %d; rows_left       = %lu\n", h, rows_left);
      outfile->Printf( "\th = %d; nbuckets        = %d\n", h, nbuckets);
      
    }

    global_dpd_->buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    global_dpd_->buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < nbuckets; n++) {
      if(nbuckets == 1) this_bucket_rows = rows_per_bucket;
      else this_bucket_rows = (n < nbuckets-1) ? rows_per_bucket : rows_left;
      global_dpd_->buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, this_bucket_rows);
      for(pq=0; pq < this_bucket_rows; pq++) {
        for(Gr=0; Gr < nirreps; Gr++) {
          Gs = h^Gr;
          nrows = moinfo.sopi[Gr];
          ncols = moinfo.actpi[Gs];
          nlinks = moinfo.sopi[Gs];
          rs = J.col_offset[h][Gr];
          if(nrows && ncols && nlinks)
            C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][pq][rs],nlinks,
                    &(C[Gs][0][C_offset[Gs]]),moinfo.mopi[Gs],0.0,TMP[0],nso);

          nrows = moinfo.actpi[Gr];
          ncols = moinfo.actpi[Gs];
          nlinks = moinfo.sopi[Gr];
          rs = K.col_offset[h][Gr];
          if(nrows && ncols && nlinks)
            C_DGEMM('t','n',nrows,ncols,nlinks,1.0,&(C[Gr][0][C_offset[Gr]]),moinfo.mopi[Gr],TMP[0],nso,
                    0.0,&K.matrix[h][pq][rs],ncols);
        } /* Gr */

        p = moinfo.pitz2corr_two[K.params->roworb[h][pq+n*rows_per_bucket][0]];
        q = moinfo.pitz2corr_two[K.params->roworb[h][pq+n*rows_per_bucket][1]];
        PQ = INDEX(p,q);
        for(rs=0; rs < K.params->coltot[h]; rs++) {
          r = moinfo.pitz2corr_two[K.params->colorb[h][rs][0]];
          s = moinfo.pitz2corr_two[K.params->colorb[h][rs][1]];
          RS = INDEX(r,s);
          if(r >= s && RS <= PQ)
            iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs], params.print_tei, "outfile", 0);
        } /* rs */
      } /* pq */
    }
    global_dpd_->buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    global_dpd_->buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  global_dpd_->buf4_close(&K);
  global_dpd_->buf4_close(&J);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  timer_off("RHF:2ndhalf");

  free_block(TMP);

  psio_close(PSIF_HALFT1, 0);

  if(params.print_lvl) {
    outfile->Printf( "\tTwo-electron integral transformation complete.\n");
    
  }
}

  } // namespace transqt2
} // namespace psi