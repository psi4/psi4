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
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi {
extern FILE* outfile;
  namespace transqt2 {

void transtwo_uhf(void)
{
  int nirreps, nso, nmo;
  double ***C, ***C_a, ***C_b, **scratch, **TMP;
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
  C_a = moinfo.C_a;
  C_b = moinfo.C_b;
  C_offset = moinfo.C_offset;

  TMP = block_matrix(nso,nso);

  /*** AA/AB two-electron integral transformation ***/

  if(params.print_lvl) {
    fprintf(outfile, "\tStarting AA/AB first half-transformation.\n");
    fflush(outfile);
  }

  psio_open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
  psio_open(PSIF_HALFT0, PSIO_OPEN_NEW);

  C = C_a; /* Use alpha MOs for this half-transformation */

  dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (pq,rs)");
  dpd_buf4_init(&K, PSIF_HALFT0, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
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
      fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memfree);
      fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rows_left);
      fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nbuckets);
      fflush(outfile);
    }

    dpd_buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < nbuckets; n++) {
      if(nbuckets == 1) this_bucket_rows = rows_per_bucket;
      else this_bucket_rows = (n < nbuckets-1) ? rows_per_bucket : rows_left;
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, this_bucket_rows);
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
      dpd_buf4_mat_irrep_wrt_block(&K, h, n*rows_per_bucket, this_bucket_rows);
    }
    dpd_buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  psio_close(PSIF_SO_PRESORT, 1); /* must keep the presort file for the upcoming BB transformation */

  if(params.print_lvl) {
    fprintf(outfile, "\tSorting AA/AB half-transformed integrals.\n");
    fflush(outfile);
  }

  psio_open(PSIF_HALFT1, PSIO_OPEN_NEW);

  dpd_buf4_init(&K, PSIF_HALFT0, 0, 3, 8, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  dpd_buf4_sort(&K, PSIF_HALFT1, rspq, 8, 3, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_close(&K);

  psio_close(PSIF_HALFT0, 0);

  if(params.print_lvl) {
    fprintf(outfile, "\tStarting AA second half-transformation.\n");
    fflush(outfile);
  }
  iwl_buf_init(&MBuff, PSIF_MO_AA_TEI, params.tolerance, 0, 0);

  C = C_a; /* Usa alpha MOs for this half-transformation */

  dpd_buf4_init(&J, PSIF_HALFT1, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_init(&K, CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)");
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
      fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memfree);
      fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rows_left);
      fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nbuckets);
      fflush(outfile);
    }

    dpd_buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < nbuckets; n++) {
      if(nbuckets == 1) this_bucket_rows = rows_per_bucket;
      else this_bucket_rows = (n < nbuckets-1) ? rows_per_bucket : rows_left;
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, this_bucket_rows);
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

        p = moinfo.pitz2corr_two_A[K.params->roworb[h][pq+n*rows_per_bucket][0]];
        q = moinfo.pitz2corr_two_A[K.params->roworb[h][pq+n*rows_per_bucket][1]];
        PQ = INDEX(p,q);
        for(rs=0; rs < K.params->coltot[h]; rs++) {
          r = moinfo.pitz2corr_two_A[K.params->colorb[h][rs][0]];
          s = moinfo.pitz2corr_two_A[K.params->colorb[h][rs][1]];
          RS = INDEX(r,s);
          if(r >= s && RS <= PQ)
            iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs], params.print_tei, outfile, 0);
        } /* rs */
      } /* pq */
    }
    dpd_buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  if(params.print_lvl) {
    fprintf(outfile, "\tStarting AB second half-transformation.\n");
    fflush(outfile);
  }
  iwl_buf_init(&MBuff, PSIF_MO_AB_TEI, params.tolerance, 0, 0);

  C = C_b; /* Usa beta MOs for this half-transformation */

  dpd_buf4_init(&J, PSIF_HALFT1, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_init(&K, CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)");
  for(h=0; h < nirreps; h++) {

    if (J.params->coltot[h]) {
      memfree = (unsigned long int) (dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
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
      fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memfree);
      fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rows_left);
      fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nbuckets);
      fflush(outfile);
    }

    dpd_buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < nbuckets; n++) {
      if(nbuckets == 1) this_bucket_rows = rows_per_bucket;
      else this_bucket_rows = (n < nbuckets-1) ? rows_per_bucket : rows_left;
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, this_bucket_rows);
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

        p = moinfo.pitz2corr_two_A[K.params->roworb[h][pq+n*rows_per_bucket][0]];
        q = moinfo.pitz2corr_two_A[K.params->roworb[h][pq+n*rows_per_bucket][1]];
        PQ = INDEX(p,q);
        for(rs=0; rs < K.params->coltot[h]; rs++) {
          r = moinfo.pitz2corr_two_B[K.params->colorb[h][rs][0]];
          s = moinfo.pitz2corr_two_B[K.params->colorb[h][rs][1]];
          RS = INDEX(r,s);
          if(r >= s)
            iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs], params.print_tei, outfile, 0);
        } /* rs */
      } /* pq */
    }
    dpd_buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  psio_close(PSIF_HALFT1, 0);

  /*** AA/AB two-electron integral transformation complete ***/

  /*** BB two-electron integral transformation ***/

  if(params.print_lvl) {
    fprintf(outfile, "\tStarting BB first half-transformation.\n");
    fflush(outfile);
  }

  psio_open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
  psio_open(PSIF_HALFT0, PSIO_OPEN_NEW);

  C = C_b; /* Use beta MOs for this half-transformation */

  dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (pq,rs)");
  dpd_buf4_init(&K, PSIF_HALFT0, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  for(h=0; h < nirreps; h++) {

    if (J.params->coltot[h]) {
      memfree = (unsigned long int) (dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
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
      fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memfree);
      fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rows_left);
      fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nbuckets);
      fflush(outfile);
    }

    dpd_buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < nbuckets; n++) {
      if(nbuckets == 1) this_bucket_rows = rows_per_bucket;
      else this_bucket_rows = (n < nbuckets-1) ? rows_per_bucket : rows_left;
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, this_bucket_rows);
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
      dpd_buf4_mat_irrep_wrt_block(&K, h, n*rows_per_bucket, this_bucket_rows);
    }
    dpd_buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  psio_close(PSIF_SO_PRESORT, 0);

  if(params.print_lvl) {
    fprintf(outfile, "\tSorting BB half-transformed integrals.\n");
    fflush(outfile);
  }

  psio_open(PSIF_HALFT1, PSIO_OPEN_NEW);

  dpd_buf4_init(&K, PSIF_HALFT0, 0, 3, 8, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  dpd_buf4_sort(&K, PSIF_HALFT1, rspq, 8, 3, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_close(&K);

  psio_close(PSIF_HALFT0, 0);

  if(params.print_lvl) {
    fprintf(outfile, "\tStarting BB second half-transformation.\n");
    fflush(outfile);
  }
  iwl_buf_init(&MBuff, PSIF_MO_BB_TEI, params.tolerance, 0, 0);

  C = C_b; /* Usa beta MOs for this half-transformation */

  dpd_buf4_init(&J, PSIF_HALFT1, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_init(&K, CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)");
  for(h=0; h < nirreps; h++) {

    if (J.params->coltot[h]) {
      memfree = (unsigned long int) (dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
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
      fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memfree);
      fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rows_left);
      fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nbuckets);
      fflush(outfile);
    }

    dpd_buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < nbuckets; n++) {
      if(nbuckets == 1) this_bucket_rows = rows_per_bucket;
      else this_bucket_rows = (n < nbuckets-1) ? rows_per_bucket : rows_left;
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, this_bucket_rows);
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

        p = moinfo.pitz2corr_two_B[K.params->roworb[h][pq+n*rows_per_bucket][0]];
        q = moinfo.pitz2corr_two_B[K.params->roworb[h][pq+n*rows_per_bucket][1]];
        PQ = INDEX(p,q);
        for(rs=0; rs < K.params->coltot[h]; rs++) {
          r = moinfo.pitz2corr_two_B[K.params->colorb[h][rs][0]];
          s = moinfo.pitz2corr_two_B[K.params->colorb[h][rs][1]];
          RS = INDEX(r,s);
          if(r >= s && RS <= PQ)
            iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs], params.print_tei, outfile, 0);
        } /* rs */
      } /* pq */
    }
    dpd_buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  psio_close(PSIF_HALFT1, 0);

  /*** BB two-electron integral transformation complete ***/

  free_block(TMP);

  if(params.print_lvl) {
    fprintf(outfile, "\tTwo-electron integral transformation complete.\n");
    fflush(outfile);
  }
}

  } // namespace transqt2
} // namespace psi
