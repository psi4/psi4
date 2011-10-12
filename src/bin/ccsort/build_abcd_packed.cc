/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <ccfiles.h>
#include <psifiles.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"
#include <psi4-dec.h>

namespace psi { namespace ccsort {

void idx_permute_multipass(dpdfile4 *File, int this_bucket,
                           int **bucket_map, unsigned long int **bucket_offset,
                           int p, int q, int r, int s,
                           int perm_pr, int perm_qs, int perm_prqs,
                           double value, FILE *outfile);

int build_abcd_packed(int inputfile, double tolerance, int keep)
{
  struct iwlbuf InBuf;
  int lastbuf;
  long int memoryb, memoryd, core_left, row_length;
  unsigned int h, nirreps, n, row, nump, numq, nbuckets;
  int **bucket_map;
  unsigned long **bucket_offset, **bucket_rowdim, **bucket_size;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s;
  int ab, cd, c, d, CD, DC;
  double value, abcd, abdc;
  struct iwlbuf *SortBuf;
  psio_address next;
  dpdfile4 B;
  dpdbuf4 B_s, B_a;
  double **B_diag;
  int rows_per_bucket, rows_left, row_start, nvirt;
  int m, Gc, C, cc;

  dpd_file4_init_nocache(&B, CC_BINTS, 0, 8, 5, "B(+) <ab|cd>");  /* junk target */
  dpd_buf4_init(&B_s, CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
  dpd_buf4_scm(&B_s, 0.0);
  dpd_buf4_init(&B_a, CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
  dpd_buf4_scm(&B_a, 0.0);

  nirreps = B.params->nirreps;

  memoryb = Process::environment.get_memory();
  memoryd = memoryb/sizeof(double);

  /* It's annoying that I have to compute this here */
  for(h=0,nump=0,numq=0; h < B.params->nirreps; h++) {
    nump += B.params->ppi[h];
    numq += B.params->qpi[h];
  }
  bucket_map = init_int_matrix(nump,numq);
  for(p=0; p < nump; p++)
    for(q=0; q < numq; q++)
      bucket_map[p][q] = -1;  /* initialize with junk */

  /* Room for one bucket to begin with */
  bucket_offset = (unsigned long int **) malloc(sizeof(unsigned long int *));
  bucket_offset[0] =
             (unsigned long int *) malloc(sizeof(unsigned long int)*nirreps);
  for(int hh=0; hh < nirreps; hh++) bucket_offset[0][hh] = 0;

  bucket_rowdim = (unsigned long int **) malloc(sizeof(unsigned long int *));
  bucket_rowdim[0] =
             (unsigned long int *) malloc(sizeof(unsigned long int)*nirreps);
  for(int hh=0; hh < nirreps; hh++) bucket_rowdim[0][hh] = 0;

  bucket_size = (unsigned long int **) malloc(sizeof(unsigned long int *));
  bucket_size[0] =
             (unsigned long int *) malloc(sizeof(unsigned long int)*nirreps);
  for(int hh=0; hh < nirreps; hh++) bucket_size[0][hh] = 0;

  /* Figure out how many passes we need and where each p,q goes */
  for(h=0,core_left=memoryd,nbuckets=1; h < nirreps; h++) {

    row_length = (long int) B.params->coltot[h^(B.my_irrep)];

    for(row=0; row < B.params->rowtot[h]; row++) {

      if((core_left - row_length) >= 0) {
        core_left -= row_length;
        bucket_rowdim[nbuckets-1][h]++;
        bucket_size[nbuckets-1][h] += row_length;
      }
      else {
        nbuckets++;
        core_left = memoryd - row_length;

        /* Make room for another bucket */
        bucket_offset = (unsigned long int **) realloc((void *) bucket_offset,
                         nbuckets * sizeof(unsigned long int *));
        bucket_offset[nbuckets-1] =
             (unsigned long int *) malloc(sizeof(unsigned long int)*nirreps);
        for(int hh=0; hh < nirreps; hh++) bucket_offset[nbuckets-1][hh] = 0;
        bucket_offset[nbuckets-1][h] = row;

        bucket_rowdim = (unsigned long int **) realloc((void *) bucket_rowdim,
                         nbuckets * sizeof(unsigned long int *));
        bucket_rowdim[nbuckets-1] =
             (unsigned long int *) malloc(sizeof(unsigned long int)*nirreps);
        for(int hh=0; hh < nirreps; hh++) bucket_rowdim[nbuckets-1][hh] = 0;
        bucket_rowdim[nbuckets-1][h] = 1;

        bucket_size = (unsigned long int **) realloc((void *) bucket_size,
                         nbuckets * sizeof(unsigned long int *));
        bucket_size[nbuckets-1] =
             (unsigned long int *) malloc(sizeof(unsigned long int)*nirreps);
        for(int hh=0; hh < nirreps; hh++) bucket_size[nbuckets-1][hh] = 0;
        bucket_size[nbuckets-1][h] = row_length;
      }

      p = B.params->roworb[h][row][0];
      q = B.params->roworb[h][row][1];
      bucket_map[p][q] = nbuckets - 1;
    }
  }

  fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", B.label, nbuckets);
  fflush(outfile);

  next = PSIO_ZERO;
  for(n=0; n < nbuckets; n++) { /* nbuckets = number of passes */

    /* Prepare target matrix */
    for(h=0; h < nirreps; h++) {
      B.matrix[h] = block_matrix(bucket_rowdim[n][h], B.params->coltot[h]);
    }

    iwl_buf_init(&InBuf, inputfile, tolerance, 1, 1);

    lblptr = InBuf.labels;
    valptr = InBuf.values;
    lastbuf = InBuf.lastbuf;

    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      idx_permute_multipass(&B,n,bucket_map,bucket_offset,
                            p,q,r,s,1,1,1,value,outfile);

    } /* end loop through current buffer */

    /* Now run through the rest of the buffers in the file */
    while (!lastbuf) {
      iwl_buf_fetch(&InBuf);
      lastbuf = InBuf.lastbuf;

      for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
        p = (int) lblptr[idx++];
        q = (int) lblptr[idx++];
        r = (int) lblptr[idx++];
        s = (int) lblptr[idx++];

        value = (double) valptr[InBuf.idx];

        idx_permute_multipass(&B,n,bucket_map,bucket_offset,
                              p,q,r,s,1,1,1,value,outfile);

      } /* end loop through current buffer */
    } /* end loop over reading buffers */

    iwl_buf_close(&InBuf, 1); /* close buffer for next pass */

    for(h=0; h < nirreps;h++) {
      if(bucket_size[n][h]) {

/* 	psio_write(B.filenum, B.label, (char *) B.matrix[h][0], */
/* 		   bucket_size[n][h]*((long int) sizeof(double)), next, &next); */

        dpd_buf4_mat_irrep_row_init(&B_s, h);
        dpd_buf4_mat_irrep_row_init(&B_a, h);
        for(ab=0; ab < bucket_rowdim[n][h]; ab++) {
          for(cd=0; cd < B_s.params->coltot[h]; cd++) {
            c = B_s.params->colorb[h][cd][0];
            d = B_s.params->colorb[h][cd][1];
            CD = B.params->colidx[c][d];
            DC = B.params->colidx[d][c];
            abcd = B.matrix[h][ab][CD];
            abdc = B.matrix[h][ab][DC];
            B_s.matrix[h][0][cd] = abcd + abdc;
            B_a.matrix[h][0][cd] = abcd - abdc;
          }
          dpd_buf4_mat_irrep_row_wrt(&B_s, h, ab+bucket_offset[n][h]);
          dpd_buf4_mat_irrep_row_wrt(&B_a, h, ab+bucket_offset[n][h]);
        }
        dpd_buf4_mat_irrep_row_close(&B_s, h);
        dpd_buf4_mat_irrep_row_close(&B_a, h);
      }
      free_block(B.matrix[h]);
    }

  } /* end loop over buckets/passes */

  /* Get rid of the input integral file */
  psio_open(inputfile, PSIO_OPEN_OLD);
  psio_close(inputfile, keep);

  free_int_matrix(bucket_map);

  for(n=0; n < nbuckets; n++) {
    free(bucket_offset[n]);
    free(bucket_rowdim[n]);
    free(bucket_size[n]);
  }
  free(bucket_offset);
  free(bucket_rowdim);
  free(bucket_size);

  dpd_file4_close(&B);
  dpd_buf4_close(&B_s);
  dpd_buf4_close(&B_a);

  /* Generate <ab|cc> components of B(+) */
  for(h=0,nvirt=0; h < moinfo.nirreps; h++) nvirt += moinfo.virtpi[h];
  dpd_buf4_init(&B_s, CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");

  rows_per_bucket = dpd_memfree()/(B_s.params->coltot[0] + nvirt);
  if(rows_per_bucket > B_s.params->rowtot[0]) rows_per_bucket = B_s.params->rowtot[0];
  nbuckets = (int) ceil((double) B_s.params->rowtot[0]/(double) rows_per_bucket);
  rows_left = B_s.params->rowtot[0] % rows_per_bucket;

  dpd_buf4_mat_irrep_init_block(&B_s, 0, rows_per_bucket);
  B_diag = dpd_block_matrix(rows_per_bucket, nvirt);

  next = PSIO_ZERO;
  for(m=0; m < (rows_left ? nbuckets-1:nbuckets); m++) {
    row_start = m * rows_per_bucket;
    dpd_buf4_mat_irrep_rd_block(&B_s, 0, row_start, rows_per_bucket);
    for(ab=0; ab < rows_per_bucket; ab++)
      for(Gc=0; Gc < moinfo.nirreps; Gc++)
        for(C=0; C < moinfo.virtpi[Gc]; C++) {
          c = moinfo.vir_off[Gc] + C;
          cc = B_s.params->colidx[c][c];
          B_diag[ab][c] = B_s.matrix[0][ab][cc];
        }
    psio_write(CC_BINTS, "B(+) <ab|cc>", (char *) B_diag[0], rows_per_bucket*nvirt*sizeof(double), next, &next);
  }
  if(rows_left) {
    row_start = m * rows_per_bucket;
    dpd_buf4_mat_irrep_rd_block(&B_s, 0, row_start, rows_left);
    for(ab=0; ab < rows_left; ab++)
      for(Gc=0; Gc < moinfo.nirreps; Gc++)
        for(C=0; C < moinfo.virtpi[Gc]; C++) {
          c = moinfo.vir_off[Gc] + C;
          cc = B_s.params->colidx[c][c];
          B_diag[ab][c] = B_s.matrix[0][ab][cc];
        }
    psio_write(CC_BINTS, "B(+) <ab|cc>", (char *) B_diag[0], rows_left*nvirt*sizeof(double), next, &next);
  }
  dpd_free_block(B_diag, rows_per_bucket, nvirt);
  dpd_buf4_mat_irrep_close_block(&B_s, 0, rows_per_bucket);
  dpd_buf4_close(&B_s);

  return 0;
}

}} // namespace psi::ccsort
