/*! \file
    \ingroup TRANSQT2
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
#define EXTERN
#include "globals.h"

namespace psi {
extern FILE* outfile;
  namespace transqt2 {

void frozen_core(int,int,int,int,double,double *,double *,double *,
                 double *,int);

void idx_permute_presort(dpdfile4 *,int,int **,int **,int,int,int,int,
                           double,FILE *);

int file_build_presort(dpdfile4 *File, int inputfile, double tolerance,
                       long int memoryb, int keep, int fzc, double *D_a,
                       double *D_b, double *fock_a, double *fock_b, int ref)
{
  struct iwlbuf InBuf;
  int lastbuf;
  long int memoryd, core_left, row_length;
  int h, nirreps, n, row, nump, numq, nbuckets;
  int **bucket_map, **bucket_offset, **bucket_rowdim;
  long int **bucket_size;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s, pq, rs;
  double value;
  struct iwlbuf *SortBuf;
  psio_address next;

  nirreps = File->params->nirreps;

  memoryd = memoryb/sizeof(double);

  /* It's annoying that I have to compute this here */
  for(h=0,nump=0,numq=0; h < nirreps; h++) {
    nump += File->params->ppi[h];
    numq += File->params->qpi[h];
  }
  bucket_map = init_int_matrix(nump,numq);

  /* Room for one bucket to begin with */
  bucket_offset = (int **) malloc(sizeof(int *));
  bucket_offset[0] = init_int_array(nirreps);
  bucket_rowdim = (int **) malloc(sizeof(int *));
  bucket_rowdim[0] = init_int_array(nirreps);
  bucket_size = (long int **) malloc(sizeof(long int *));
  bucket_size[0] = init_long_int_array(nirreps);

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
        bucket_offset = (int **) realloc((void *) bucket_offset,
                                         nbuckets * sizeof(int *));
        bucket_offset[nbuckets-1] = init_int_array(nirreps);
        bucket_offset[nbuckets-1][h] = row;

        bucket_rowdim = (int **) realloc((void *) bucket_rowdim,
                                         nbuckets * sizeof(int *));
        bucket_rowdim[nbuckets-1] = init_int_array(nirreps);
        bucket_rowdim[nbuckets-1][h] = 1;

        bucket_size = (long int **) realloc((void *) bucket_size,
                                            nbuckets * sizeof(long int *));
        bucket_size[nbuckets-1] = init_long_int_array(nirreps);
        bucket_size[nbuckets-1][h] = row_length;
      }

      p = File->params->roworb[h][row][0];
      q = File->params->roworb[h][row][1];
      bucket_map[p][q] = nbuckets - 1;
    }
  }

  if(params.print_lvl) {
    fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", File->label, nbuckets);
    fflush(outfile);
  }

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

      idx_permute_presort(File,n,bucket_map,bucket_offset,p,q,r,s,value,outfile);

      if(fzc && !n) /* build frozen-core operator only on first pass*/
        frozen_core(p,q,r,s,value,D_a,D_b,fock_a,fock_b,ref);

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

        idx_permute_presort(File,n,bucket_map,bucket_offset,p,q,r,s,value,outfile);

        if(fzc && !n) /* build frozen-core operator only on first pass */
          frozen_core(p,q,r,s,value,D_a,D_b,fock_a,fock_b,ref);

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

void frozen_core(int p, int q, int r, int s, double value, double *D_a, double *D_b,
                 double *fock_a, double *fock_b, int ref)
{
  int a, b, c, d, ab, cd, ad, bc;
  int dum,found=0;
  int al[8], bl[8], cl[8], dl[8];

  if(ref==0 || ref==1) { /* RHF/ROHF */
    a = al[0] = p;
    b = bl[0] = q;
    c = cl[0] = r;
    d = dl[0] = s;
    ab = INDEX(a,b);
    cd = INDEX(c,d);
    bc = INDEX(b,c);
    ad = INDEX(a,d);
    fock_a[cd] += 2.0 * D_a[ab] * value;
    if(b>=c) fock_a[bc] -= D_a[ad] * value;

    a = al[1] = q;
    b = bl[1] = p;
    c = cl[1] = r;
    d = dl[1] = s;
    if(!(a==al[0] && b==bl[0] && c==cl[0] && d==dl[0])) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) fock_a[cd] += 2.0 * D_a[ab] * value;
      if(b>=c) fock_a[bc] -= D_a[ad] * value;
    }

    a = al[2] = p;
    b = bl[2] = q;
    c = cl[2] = s;
    d = dl[2] = r;
    for(dum=0,found=0; dum < 2 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) fock_a[cd] += 2.0 * D_a[ab] * value;
      if(b>=c) fock_a[bc] -= D_a[ad] * value;
    }

    a = al[3] = q;
    b = bl[3] = p;
    c = cl[3] = s;
    d = dl[3] = r;
    for(dum=0,found=0; dum < 3 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) fock_a[cd] += 2.0 * D_a[ab] * value;
      if(b>=c) fock_a[bc] -= D_a[ad] * value;
    }

    a = al[4] = r;
    b = bl[4] = s;
    c = cl[4] = p;
    d = dl[4] = q;
    for(dum=0,found=0; dum < 4 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) fock_a[cd] += 2.0 * D_a[ab] * value;
      if(b>=c) fock_a[bc] -= D_a[ad] * value;
    }

    a = al[5] = r;
    b = bl[5] = s;
    c = cl[5] = q;
    d = dl[5] = p;
    for(dum=0,found=0; dum < 5 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) fock_a[cd] += 2.0 * D_a[ab] * value;
      if(b>=c) fock_a[bc] -= D_a[ad] * value;
    }

    a = al[6] = s;
    b = bl[6] = r;
    c = cl[6] = p;
    d = dl[6] = q;
    for(dum=0,found=0; dum < 6 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) fock_a[cd] += 2.0 * D_a[ab] * value;
      if(b>=c) fock_a[bc] -= D_a[ad] * value;
    }

    a = al[7] = s;
    b = bl[7] = r;
    c = cl[7] = q;
    d = dl[7] = p;
    for(dum=0,found=0; dum < 7 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) fock_a[cd] += 2.0 * D_a[ab] * value;
      if(b>=c) fock_a[bc] -= D_a[ad] * value;
    }
  }
  else { /* UHF */
    a = al[0] = p;
    b = bl[0] = q;
    c = cl[0] = r;
    d = dl[0] = s;
    ab = INDEX(a,b);
    cd = INDEX(c,d);
    bc = INDEX(b,c);
    ad = INDEX(a,d);
    fock_a[cd] += (D_a[ab] + D_b[ab]) * value;
    fock_b[cd] += (D_a[ab] + D_b[ab]) * value;
    if(b>=c) {
      fock_a[bc] -= D_a[ad] * value;
      fock_b[bc] -= D_b[ad] * value;
    }

    a = al[1] = q;
    b = bl[1] = p;
    c = cl[1] = r;
    d = dl[1] = s;
    if(!(a==al[0] && b==bl[0] && c==cl[0] && d==dl[0])) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) {
        fock_a[cd] += (D_a[ab] + D_b[ab]) * value;
        fock_b[cd] += (D_a[ab] + D_b[ab]) * value;
      }
      if(b>=c) {
        fock_a[bc] -= D_a[ad] * value;
        fock_b[bc] -= D_b[ad] * value;
      }
    }

    a = al[2] = p;
    b = bl[2] = q;
    c = cl[2] = s;
    d = dl[2] = r;
    for(dum=0,found=0; dum < 2 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) {
        fock_a[cd] += (D_a[ab] + D_b[ab]) * value;
        fock_b[cd] += (D_a[ab] + D_b[ab]) * value;
      }
      if(b>=c) {
        fock_a[bc] -= D_a[ad] * value;
        fock_b[bc] -= D_b[ad] * value;
      }
    }

    a = al[3] = q;
    b = bl[3] = p;
    c = cl[3] = s;
    d = dl[3] = r;
    for(dum=0,found=0; dum < 3 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) {
        fock_a[cd] += (D_a[ab] + D_b[ab]) * value;
        fock_b[cd] += (D_a[ab] + D_b[ab]) * value;
      }
      if(b>=c) {
        fock_a[bc] -= D_a[ad] * value;
        fock_b[bc] -= D_b[ad] * value;
      }
    }

    a = al[4] = r;
    b = bl[4] = s;
    c = cl[4] = p;
    d = dl[4] = q;
    for(dum=0,found=0; dum < 4 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) {
        fock_a[cd] += (D_a[ab] + D_b[ab]) * value;
        fock_b[cd] += (D_a[ab] + D_b[ab]) * value;
      }
      if(b>=c) {
        fock_a[bc] -= D_a[ad] * value;
        fock_b[bc] -= D_b[ad] * value;
      }
    }

    a = al[5] = r;
    b = bl[5] = s;
    c = cl[5] = q;
    d = dl[5] = p;
    for(dum=0,found=0; dum < 5 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) {
        fock_a[cd] += (D_a[ab] + D_b[ab]) * value;
        fock_b[cd] += (D_a[ab] + D_b[ab]) * value;
      }
      if(b>=c) {
        fock_a[bc] -= D_a[ad] * value;
        fock_b[bc] -= D_b[ad] * value;
      }
    }

    a = al[6] = s;
    b = bl[6] = r;
    c = cl[6] = p;
    d = dl[6] = q;
    for(dum=0,found=0; dum < 6 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) {
        fock_a[cd] += (D_a[ab] + D_b[ab]) * value;
        fock_b[cd] += (D_a[ab] + D_b[ab]) * value;
      }
      if(b>=c) {
        fock_a[bc] -= D_a[ad] * value;
        fock_b[bc] -= D_b[ad] * value;
      }
    }

    a = al[7] = s;
    b = bl[7] = r;
    c = cl[7] = q;
    d = dl[7] = p;
    for(dum=0,found=0; dum < 7 && !found; dum++)
      if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
    if(!found) {
      ab = INDEX(a,b);
      cd = INDEX(c,d);
      bc = INDEX(b,c);
      ad = INDEX(a,d);
      if(c>=d) {
        fock_a[cd] += (D_a[ab] + D_b[ab]) * value;
        fock_b[cd] += (D_a[ab] + D_b[ab]) * value;
      }
      if(b>=c) {
        fock_a[bc] -= D_a[ad] * value;
        fock_b[bc] -= D_b[ad] * value;
      }
    }
  }
}

  } // namespace transqt2
} // namespace psi
