/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
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
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psi4-dec.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void idx_permute(dpdfile4 *File, struct iwlbuf *OutBuf,
                 int **bucket_map, int p, int q, int r, int s,
                 int perm_pr, int perm_qs, int perm_prqs,
                 double value, std::string out);

int file_build(dpdfile4 *File, int inputfile, double tolerance,
               int perm_pr, int perm_qs, int perm_prqs, int keep)
{
  struct iwlbuf InBuf;
  int lastbuf;
  long int memoryb, memoryd;
  int h, nirreps, n, row, col, nump, numq, row_length, core_left, nbuckets;
  int **bucket_map, **bucket_offset, **bucket_rowdim, **bucket_size, offset;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s;
  double value;
  struct iwlbuf *SortBuf, *OutBuf;
  psio_address next;

  nirreps = File->params->nirreps;

  memoryb = Process::environment.get_memory();
  //fndcor(&memoryb, infile, outfile);
  memoryd = memoryb/sizeof(double);

  /* It's annoying that I have to compute this here */
  for(h=0,nump=0,numq=0; h < File->params->nirreps; h++) {
      nump += File->params->ppi[h];
      numq += File->params->qpi[h];
    }
  bucket_map = init_int_matrix(nump,numq);

  /* Room for one bucket to begin with */
  bucket_offset = (int **) malloc(sizeof(int *));
  bucket_offset[0] = init_int_array(nirreps);
  bucket_rowdim = (int **) malloc(sizeof(int *));
  bucket_rowdim[0] = init_int_array(nirreps);
  bucket_size = (int **) malloc(sizeof(int *));
  bucket_size[0] = init_int_array(nirreps);

  /* Figure out how many buckets we need and where each p,q goes */
  for(h=0,core_left=memoryd,nbuckets=1; h < nirreps; h++) {

      row_length = File->params->coltot[h^(File->my_irrep)];

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
	      int **p;

	      p = static_cast<int **>(realloc(static_cast<void *>(bucket_offset),
                                              nbuckets * sizeof(int *)));
	      if(p == nullptr) {
		throw PsiException("file_build: allocation error", __FILE__, __LINE__);
	      } else {
		bucket_offset = p;
	      }
	      bucket_offset[nbuckets-1] = init_int_array(nirreps);
              bucket_offset[nbuckets-1][h] = row;

	      p = static_cast<int **>(realloc(static_cast<void *>(bucket_rowdim),
                                              nbuckets * sizeof(int *)));
	      if(p == nullptr) {
		throw PsiException("file_build: allocation error", __FILE__, __LINE__);
	      } else {
		bucket_rowdim = p;
	      }
              bucket_rowdim[nbuckets-1] = init_int_array(nirreps);
              bucket_rowdim[nbuckets-1][h] = 1;

	      p = static_cast<int **>(realloc(static_cast<void *>(bucket_size),
                                              nbuckets * sizeof(int *)));
	      if(p == nullptr) {
		throw PsiException("file_build: allocation error", __FILE__, __LINE__);
	      } else {
		bucket_size = p;
	      }
              bucket_size[nbuckets-1] = init_int_array(nirreps);
              bucket_size[nbuckets-1][h] = row_length;
            }

          p = File->params->roworb[h][row][0];
          q = File->params->roworb[h][row][1];
          bucket_map[p][q] = nbuckets - 1;
        }
    }

  outfile->Printf( "\tSorting File: %s nbuckets = %d\n",
          File->label, nbuckets);

  /* Set up IWL buffers for sorting */
  SortBuf = (struct iwlbuf *) malloc(nbuckets * sizeof(struct iwlbuf));
  for(n=0; n < nbuckets; n++)
      iwl_buf_init(&SortBuf[n], PSIF_CC_MAX+1+n, tolerance, 0, 0);

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

      idx_permute(File,SortBuf,bucket_map,p,q,r,s,
                  perm_pr,perm_qs,perm_prqs,value,"outfile");

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

      idx_permute(File,SortBuf,bucket_map,p,q,r,s,
                  perm_pr,perm_qs,perm_prqs,value,"outfile");

      } /* end loop through current buffer */
    } /* end loop over reading buffers */

  iwl_buf_close(&InBuf, keep);

  for(n=0; n < nbuckets; n++) {
      iwl_buf_flush(&SortBuf[n], 1);
      iwl_buf_close(&SortBuf[n], 1);
    }

  free_int_matrix(bucket_map);

  /* Now sort each buffer and send it to the final target */
  next = PSIO_ZERO;
  for(n=0; n < nbuckets; n++) {
      iwl_buf_init(&SortBuf[n], PSIF_CC_MAX+1+n, tolerance, 1, 0);
      lastbuf = 0;

      for(h=0; h < nirreps; h++) {
          File->matrix[h] = block_matrix(bucket_rowdim[n][h],
                                         File->params->coltot[h]);
        }

      while(!lastbuf) {
          iwl_buf_fetch(&SortBuf[n]);
          lastbuf = SortBuf[n].lastbuf;

          for(idx=4*SortBuf[n].idx; SortBuf[n].idx < SortBuf[n].inbuf;
              SortBuf[n].idx++) {
              p = (int) SortBuf[n].labels[idx++];
              q = (int) SortBuf[n].labels[idx++];
              r = (int) SortBuf[n].labels[idx++];
              s = (int) SortBuf[n].labels[idx++];

              value = (double) SortBuf[n].values[SortBuf[n].idx];

              h = File->params->psym[p] ^ File->params->qsym[q];

              row = File->params->rowidx[p][q];
              col = File->params->colidx[r][s];

              offset = bucket_offset[n][h];

              File->matrix[h][row-offset][col] = value;
            }
        }


      for(h=0; h < nirreps;h++) {
          if(bucket_size[n][h])
             psio_write(File->filenum, File->label, (char *) File->matrix[h][0],
                        bucket_size[n][h]*sizeof(double), next, &next);
          free_block(File->matrix[h]);
        }

      iwl_buf_close(&SortBuf[n], 0);
    }

  for(n=0; n < nbuckets; n++) {
      free(bucket_offset[n]);
      free(bucket_rowdim[n]);
      free(bucket_size[n]);
    }
  free(bucket_offset);
  free(bucket_rowdim);
  free(bucket_size);

  free(SortBuf);

  return 0;
}

}} // namespace psi::ccdensity
