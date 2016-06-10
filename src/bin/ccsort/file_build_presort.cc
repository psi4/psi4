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
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#include <psi4-dec.h>
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void idx_permute_presort(dpdfile4 *File, int this_bucket,
                           int **bucket_map, unsigned long int **bucket_offset,
                           int p, int q, int r, int s,
                           double value, std::string OutFileRMR);

int file_build_presort(dpdfile4 *File, int inputfile, double tolerance, int keep)
{
  struct iwlbuf InBuf;
  int lastbuf;
  long int memoryb, memoryd, core_left, row_length;
  unsigned int h, nirreps, n, row, nump, numq, nbuckets;
  int **bucket_map;
  unsigned long int **bucket_offset, **bucket_rowdim, **bucket_size;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s;
  double value;
  struct iwlbuf *SortBuf;
  psio_address next;

  nirreps = File->params->nirreps;

  memoryb = Process::environment.get_memory();
  memoryd = memoryb/sizeof(double);

  /* It's annoying that I have to compute this here */
  for(h=0,nump=0,numq=0; h < File->params->nirreps; h++) {
    nump += File->params->ppi[h];
    numq += File->params->qpi[h];
  }
  bucket_map = init_int_matrix(nump,numq);

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

    row_length = (long int) File->params->coltot[h^(File->my_irrep)];

    for(row=0; row < File->params->rowtot[h]; row++) {

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

      p = File->params->roworb[h][row][0];
      q = File->params->roworb[h][row][1];
      bucket_map[p][q] = nbuckets - 1;
    }
  }

  outfile->Printf( "\tSorting File: %s nbuckets = %d\n", File->label, nbuckets);
  

  next = PSIO_ZERO;
  for(n=0; n < nbuckets; n++) { /* nbuckets = number of passes */

    /* Prepare target matrix */
    for(h=0; h < nirreps; h++) {
      File->matrix[h] = block_matrix(bucket_rowdim[n][h], File->params->coltot[h]);
    }

    iwl_buf_init(&InBuf, inputfile, tolerance, 1, 1);

    lblptr = InBuf.labels;
    valptr = InBuf.values;
    lastbuf = InBuf.lastbuf;

    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

/*       outfile->Printf( "%d %d %d %d %20.12f\n", p, q, r, s, value); */

      idx_permute_presort(File,n,bucket_map,bucket_offset,p,q,r,s,value,"outfile");

    } /* end loop through current buffer */

    /* Now run through the rest of the buffers in the file */
    while (!lastbuf) {
      iwl_buf_fetch(&InBuf);
      lastbuf = InBuf.lastbuf;

      for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
        p = abs((int) lblptr[idx++]);
        q = (int) lblptr[idx++];
        r = (int) lblptr[idx++];
        s = (int) lblptr[idx++];

        value = (double) valptr[InBuf.idx];

/* 	outfile->Printf( "%d %d %d %d %20.12f\n", p, q, r, s, value); */

        idx_permute_presort(File,n,bucket_map,bucket_offset,p,q,r,s,value,"outfile");

      } /* end loop through current buffer */
    } /* end loop over reading buffers */

    iwl_buf_close(&InBuf, 1); /* close buffer for next pass */

    for(h=0; h < nirreps;h++) {
      if(bucket_size[n][h])
        psio_write(File->filenum, File->label, (char *) File->matrix[h][0],
                   bucket_size[n][h]*((long int) sizeof(double)), next, &next);
      free_block(File->matrix[h]);
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

  return 0;
}

}} // namespace psi::ccsort