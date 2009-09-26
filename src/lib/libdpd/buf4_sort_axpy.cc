/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {

/*
** dpd_buf4_sort_axpy(): A general DPD buffer sorting function that also adds 
** the result to a target dpdbuf4 that already exists.  Like buf4_sort(), this will
** (eventually) handle all 24 possible permutations of four-index
** buffers.  This code assumes that all symmetry-blocks of both source and
** target buffers can be stored in core.
**
** Arguments:
**   dpdbuf4 *InBuf: A pointer to the alread-initialized input
**     buffer.
**   int outfilenum: The PSI unit number for the target data.
**   enum indices index: The desired sorting pattern (see dpd.h).
**   int pqnum: The index combination for the bra indices for the new
**     dpd file4.  N.B. this is NOT error checked for consistency
**     with the input buffer.
**   int rsnum: The index combination for the ket indices for the new
**     dpd file4.  N.B. this is NOT error checked for consistency
**     with the input buffer.
**   char *label: A string labelling for this buffer.
**   double alpha: A value by which the InBuf data will be multiplied before 
**     it is added to the OutBuf.
**
** Note that the notation for some of these multi-index sorts may be
** confusing.  When the caller requests ordering "psqr" for example,
** the desired arrangement returned is indeed Out[ps][qr] =
** In[pq][rs].  However, the index notation *within the code below*
** will appear as Out[pq][rs] = In[pr][sq]. To compute this easily,
** take the desired ordering psqr, translate its indices to p->p,
** s->q, q->r, and r->s, and finally write the source indices as
** prsq.  (This "problem" arises because I want to use the same index
** notation Out[pq][rs] for the target for every rearrangement.  For
** all pair permutations, the notation is straightforward. )
**
** -Daniel, June 2002
**
** NB: Because the output dpdbuf4 is initialized within this code, one SHOULD NOT
** use this function to sort a buffer into itself.
**
** Added fully out-of-core (multipass) sorting algorithms to a number 
** of important cases: pqsr, prqs, prsq, rpqs, srqp, and rspq.  These 
** use the maximum amount of core, but require multiple passes through
** the input integral list.  These have been tested only with RHF orbitals.
** I need to add comparable code to buf4_sort.c.
**
** -TDC, March 2004
**
**
** Added OOC algorithm for rpsq sort.
** -TDC, August 2005
*/

int dpd_buf4_sort_axpy(dpdbuf4 *InBuf, int outfilenum, enum indices index,
		       int pqnum, int rsnum, const char *label, double alpha)
{
  int h,nirreps, row, col, my_irrep, r_irrep;
  int p, q, r, s, P, Q, R, S, pq, rs, sr, pr, qs, qp, rq, qr, ps, sp, rp, sq;
  int Gp, Gq, Gr, Gs, Gpq, Grs, Gsr, Gpr, Gqs, Grq, Gqr, Gps, Gsp, Grp, Gsq;
  dpdbuf4 OutBuf;
  long int rowtot, coltot, core_total, maxrows;
  int incore;
  int Grow, Gcol;
  int out_rows_per_bucket, out_nbuckets, out_rows_left, out_row_start, n;
  int in_rows_per_bucket, in_nbuckets, in_rows_left, in_row_start, m;
  int rows_per_bucket, nbuckets, rows_left, row_start;

  nirreps = InBuf->params->nirreps;
  my_irrep = InBuf->file.my_irrep;

#ifdef DPD_TIMER
  timer_on("buf4_sort");
#endif

  dpd_buf4_init(&OutBuf, outfilenum, my_irrep, pqnum, rsnum,
		pqnum, rsnum, 0, label);

  /* select in-core vs. out-of-core algorithms */
  incore = 1;
  core_total = 0;
  for(h=0; h < nirreps; h++) {
    coltot = InBuf->params->coltot[h^my_irrep];
    if(coltot) {
      maxrows = DPD_BIGNUM/coltot;
      if(maxrows < 1) {
	fprintf(stderr, "\nLIBDPD Error: too many rows in buf4_sort_axpy.\n");
	dpd_error("buf4_sort_axpy", stderr);
      }
    }
    else maxrows = DPD_BIGNUM;
    rowtot = InBuf->params->rowtot[h];
    for(; rowtot > maxrows; rowtot -= maxrows) {
      if(core_total > (core_total + 2*maxrows*coltot)) incore = 0;
      else core_total += 2*maxrows*coltot;
    }
    if(core_total > (core_total + 2*rowtot*coltot)) incore = 0;
    core_total += 2*rowtot*coltot;
  }
  if(core_total > dpd_memfree()) incore = 0;

  /* Init input and output buffers and read in all blocks of both */
#ifdef DPD_TIMER
  if(index == rspq) timer_on("axpy:alloc");
#endif
  if(incore) {
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&OutBuf, h);
      dpd_buf4_mat_irrep_rd(&OutBuf, h);

      dpd_buf4_mat_irrep_init(InBuf, h);
      dpd_buf4_mat_irrep_rd(InBuf, h);
    }
  }
#ifdef DPD_TIMER
  if(index == rspq) timer_off("axpy:alloc");
#endif

  switch(index) {
  case pqrs:
    fprintf(stderr, "\nDPD sort error: invalid index ordering.\n");
    dpd_error("buf_sort", stderr);
    break;

  case pqsr:

#ifdef DPD_TIMER
    timer_on("pqsr");
#endif

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;

	/* p->p; q->q; s->r; r->s = pqsr */
      
	for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	  p = OutBuf.params->roworb[h][pq][0];
	  q = OutBuf.params->roworb[h][pq][1];

	  row = InBuf->params->rowidx[p][q];
	      
	  for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	    r = OutBuf.params->colorb[r_irrep][rs][0];
	    s = OutBuf.params->colorb[r_irrep][rs][1];

	    sr = InBuf->params->colidx[s][r];
	      
	    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[h][row][sr];
	  }
	}
      }
    }
    else {

      for(Gpq=0; Gpq < nirreps; Gpq++) {
	Grs = Gpq^my_irrep;

	/* determine how many rows of OutBuf/InBuf we can store in half of the core */
	rows_per_bucket = dpd_memfree()/(2 * OutBuf.params->coltot[Grs]);
	if(rows_per_bucket > OutBuf.params->rowtot[Gpq])
	  rows_per_bucket = OutBuf.params->rowtot[Gpq];
	nbuckets = (int) ceil((double) OutBuf.params->rowtot[Gpq]/(double) rows_per_bucket);
	rows_left = OutBuf.params->rowtot[Gpq] % rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(&OutBuf, Gpq, rows_per_bucket);
	dpd_buf4_mat_irrep_init_block(InBuf, Gpq, rows_per_bucket);

	for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {

	  row_start = n*rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, row_start, rows_per_bucket);
	  dpd_buf4_mat_irrep_rd_block(InBuf, Gpq, row_start, rows_per_bucket);

	  for(pq=0; pq < rows_per_bucket; pq++) {
	    for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
	      r = OutBuf.params->colorb[Grs][rs][0];
	      s = OutBuf.params->colorb[Grs][rs][1];
	      sr = InBuf->params->colidx[s][r];
	      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Gpq][pq][sr];
	    }
	  }


	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, row_start, rows_per_bucket);
	}
	if(rows_left) {
	  row_start = n*rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, row_start, rows_left);
	  dpd_buf4_mat_irrep_rd_block(InBuf, Gpq, row_start, rows_left);

	  for(pq=0; pq < rows_per_bucket; pq++) {
	    for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
	      r = OutBuf.params->colorb[Grs][rs][0];
	      s = OutBuf.params->colorb[Grs][rs][1];
	      sr = InBuf->params->colidx[s][r];
	      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Gpq][pq][sr];
	    }
	  }
	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, row_start, rows_left);
	}

	dpd_buf4_mat_irrep_close_block(&OutBuf, Gpq, rows_per_bucket);
	dpd_buf4_mat_irrep_close_block(InBuf, Gpq, rows_per_bucket);
      }
    }

#ifdef DPD_TIMER
    timer_off("pqsr");
#endif
    break;

  case prqs:

#ifdef DPD_TIMER
    timer_on("prqs");
#endif

    /* p->p; r->q; q->r; s->s = prqs */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;

	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    /* Irreps on the source */
	    Gpr = Gp^Gr;  Gqs = Gq^Gs;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  pr = InBuf->params->rowidx[P][R];
			  
		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    qs = InBuf->params->colidx[Q][S];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gpr][pr][qs];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {

      for(Gpq=0; Gpq < nirreps; Gpq++) {
	Grs = Gpq^my_irrep;

	/* determine how many rows of OutBuf we can store in half of the core */
	out_rows_per_bucket = dpd_memfree()/(2 * OutBuf.params->coltot[Grs]);
	if(out_rows_per_bucket > OutBuf.params->rowtot[Gpq])
	  out_rows_per_bucket = OutBuf.params->rowtot[Gpq];
	out_nbuckets = (int) ceil((double) OutBuf.params->rowtot[Gpq]/(double) out_rows_per_bucket);
	out_rows_left = OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

	for(n=0; n < (out_rows_left ? out_nbuckets-1 : out_nbuckets); n++) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	  for(Grow=0; Grow < nirreps; Grow++) {
	    Gcol = Grow^my_irrep;

	    /* determine how many rows of InBuf we can store in the other half of the core */
	    in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gcol]);
	    if(in_rows_per_bucket > InBuf->params->rowtot[Grow])
	      in_rows_per_bucket = InBuf->params->rowtot[Grow];
	    in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grow]/(double) in_rows_per_bucket);
	    in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

	    /* allocate space for the bucket of rows */
	    dpd_buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

	    for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

	      for(pq=0; pq < out_rows_per_bucket; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gpr = Gp^Gr;

		  if(Gpr == Grow) {
		    pr = InBuf->params->rowidx[p][r] - in_row_start;
		    /* check if the current value is in the current in_bucket or not */
		    if(pr >= 0 && pr < in_rows_per_bucket) {
		      qs = InBuf->params->colidx[q][s];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][pr][qs];
		    }
		  }
		}
	      }

	    }
	    if(in_rows_left) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

	      for(pq=0; pq < out_rows_per_bucket; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gpr = Gp^Gr;

		  if(Gpr == Grow) {
		    pr = InBuf->params->rowidx[p][r] - in_row_start;
		    /* check if the current value is in core or not */
		    if(pr >= 0 && pr < in_rows_left) {
		      qs = InBuf->params->colidx[q][s];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][pr][qs];
		    }
		  }
		}
	      }
	    }

	    dpd_buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	}
	if(out_rows_left) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_left);

	  for(Grow=0; Grow < nirreps; Grow++) {
	    Gcol = Grow^my_irrep;

	    /* determine how many rows of InBuf we can store in the other half of the core */
	    in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gcol]);
	    if(in_rows_per_bucket > InBuf->params->rowtot[Grow])
	      in_rows_per_bucket = InBuf->params->rowtot[Grow];
	    in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grow]/(double) in_rows_per_bucket);
	    in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

	    /* allocate space for the bucket of rows */
	    dpd_buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

	    for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

	      for(pq=0; pq < out_rows_left; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gpr = Gp^Gr;

		  if(Gpr == Grow) {
		    pr = InBuf->params->rowidx[p][r] - in_row_start;
		    /* check if the current value is in the current in_bucket or not */
		    if(pr >= 0 && pr < in_rows_per_bucket) {
		      qs = InBuf->params->colidx[q][s];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][pr][qs];
		    }
		  }
		}
	      }

	    }
	    if(in_rows_left) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

	      for(pq=0; pq < out_rows_left; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gpr = Gp^Gr;

		  if(Gpr == Grow) {
		    pr = InBuf->params->rowidx[p][r] - in_row_start;
		    /* check if the current value is in core or not */
		    if(pr >= 0 && pr < in_rows_left) {
		      qs = InBuf->params->colidx[q][s];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][pr][qs];
		    }
		  }
		}
	      }
	    }

	    dpd_buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_left);
	}

	dpd_buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);
      }

    }
	  

#ifdef DPD_TIMER
    timer_off("prqs");
#endif
    break;

  case prsq:

#ifdef DPD_TIMER
    timer_on("prsq");
#endif

    /* p->p; r->q; s->r; q->s = psqr */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;

	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Gps = Gp^Gs;  Gqr = Gq^Gr;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  qr = InBuf->params->colidx[Q][R];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    ps = InBuf->params->rowidx[P][S];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gps][ps][qr];

		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {

      for(Gpq=0; Gpq < nirreps; Gpq++) {
	Grs = Gpq^my_irrep;

	/* determine how many rows of OutBuf we can store in half of the core */
	out_rows_per_bucket = dpd_memfree()/(2 * OutBuf.params->coltot[Grs]);
	if(out_rows_per_bucket > OutBuf.params->rowtot[Gpq])
	  out_rows_per_bucket = OutBuf.params->rowtot[Gpq];
	out_nbuckets = (int) ceil((double) OutBuf.params->rowtot[Gpq]/(double) out_rows_per_bucket);
	out_rows_left = OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

	for(n=0; n < (out_rows_left ? out_nbuckets-1 : out_nbuckets); n++) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	  for(Grow=0; Grow < nirreps; Grow++) {
	    Gcol = Grow^my_irrep;

	    /* determine how many rows of InBuf we can store in the other half of the core */
	    in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gcol]);
	    if(in_rows_per_bucket > InBuf->params->rowtot[Grow])
	      in_rows_per_bucket = InBuf->params->rowtot[Grow];
	    in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grow]/(double) in_rows_per_bucket);
	    in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

	    /* allocate space for the bucket of rows */
	    dpd_buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

	    for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

	      for(pq=0; pq < out_rows_per_bucket; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gps = Gp^Gs;

		  if(Gps == Grow) {
		    ps = InBuf->params->rowidx[p][s] - in_row_start;
		    /* check if the current value is in the current in_bucket or not */
		    if(ps >= 0 && ps < in_rows_per_bucket) {
		      qr = InBuf->params->colidx[q][r];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][ps][qr];
		    }
		  }
		}
	      }

	    }
	    if(in_rows_left) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

	      for(pq=0; pq < out_rows_per_bucket; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gps = Gp^Gs;

		  if(Gps == Grow) {
		    ps = InBuf->params->rowidx[p][s] - in_row_start;
		    /* check if the current value is in core or not */
		    if(ps >= 0 && ps < in_rows_left) {
		      qr = InBuf->params->colidx[q][r];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][ps][qr];
		    }
		  }
		}
	      }
	    }

	    dpd_buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	}
	if(out_rows_left) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_left);

	  for(Grow=0; Grow < nirreps; Grow++) {
	    Gcol = Grow^my_irrep;

	    /* determine how many rows of InBuf we can store in the other half of the core */
	    in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gcol]);
	    if(in_rows_per_bucket > InBuf->params->rowtot[Grow])
	      in_rows_per_bucket = InBuf->params->rowtot[Grow];
	    in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grow]/(double) in_rows_per_bucket);
	    in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

	    /* allocate space for the bucket of rows */
	    dpd_buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

	    for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

	      for(pq=0; pq < out_rows_left; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gps = Gp^Gs;

		  if(Gps == Grow) {
		    ps = InBuf->params->rowidx[p][s] - in_row_start;
		    /* check if the current value is in the current in_bucket or not */
		    if(ps >= 0 && ps < in_rows_per_bucket) {
		      qr = InBuf->params->colidx[q][r];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][ps][qr];
		    }
		  }
		}
	      }

	    }
	    if(in_rows_left) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

	      for(pq=0; pq < out_rows_left; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gps = Gp^Gs;

		  if(Gps == Grow) {
		    ps = InBuf->params->rowidx[p][s] - in_row_start;
		    /* check if the current value is in core or not */
		    if(ps >= 0 && ps < in_rows_left) {
		      qr = InBuf->params->colidx[q][r];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][ps][qr];
		    }
		  }
		}
	      }
	    }

	    dpd_buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_left);
	}

	dpd_buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);
      }

    }

#ifdef DPD_TIMER
    timer_off("prsq");
#endif
    break;

  case psqr:

#ifdef DPD_TIMER
    timer_on("psqr");
#endif

    /* p->p; s->q; q->r; r->s = prsq */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	     
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Gpr = Gp^Gr;  Gsq = Gs^Gq;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  pr = InBuf->params->rowidx[P][R];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    sq = InBuf->params->colidx[S][Q];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gpr][pr][sq];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for psqr sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }

#ifdef DPD_TIMER
    timer_off("psqr");
#endif
    break;

  case psrq:

#ifdef DPD_TIMER
    timer_on("psrq");
#endif

    /* p->p; s->q; r->r; q->s = psrq */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Gps = Gp^Gs;  Grq = Gr^Gq;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  rq = InBuf->params->colidx[R][Q];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    ps = InBuf->params->rowidx[P][S];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gps][ps][rq];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for psrq sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }

#ifdef DPD_TIMER
    timer_off("psrq");
#endif
    break;

  case qprs:

#ifdef DPD_TIMER
    timer_on("qprs");
#endif

    /* q->p; p->q; r->r; s->s = qprs */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;

	for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	  p = OutBuf.params->roworb[h][pq][0];
	  q = OutBuf.params->roworb[h][pq][1];
	  qp = InBuf->params->rowidx[q][p];

	  for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	    r = OutBuf.params->colorb[r_irrep][rs][0];
	    s = OutBuf.params->colorb[r_irrep][rs][1];

	    col = InBuf->params->colidx[r][s];
		  
	    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[h][qp][col];
	  }
	}

      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for qprs sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }

#ifdef DPD_TIMER
    timer_off("qprs");
#endif
    break;

  case qpsr:

#ifdef DPD_TIMER
    timer_on("qpsr");
#endif

    /* q->p; p->q; s->r; r->s = qpsr */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;

	for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	  p = OutBuf.params->roworb[h][pq][0];
	  q = OutBuf.params->roworb[h][pq][1];
	  qp = InBuf->params->rowidx[q][p];

	  for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	    r = OutBuf.params->colorb[r_irrep][rs][0];
	    s = OutBuf.params->colorb[r_irrep][rs][1];
	    sr = InBuf->params->colidx[s][r];
		  
	    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[h][qp][sr];
	  }
	}

      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for qpsr sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }

#ifdef DPD_TIMER
    timer_off("qpsr");
#endif
    break;

  case qrps:
#ifdef DPD_TIMER
    timer_on("qrps");
#endif

    /* q->p; r->q; p->r; s->s = rpqs */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Grp = Gr^Gp; Gqs = Gq^Gs;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  rp = InBuf->params->rowidx[R][P];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    qs = InBuf->params->colidx[Q][S];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Grp][rp][qs];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }    
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for qrps sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }

#ifdef DPD_TIMER
    timer_off("qrps");
#endif
    break;

  case qrsp:

#ifdef DPD_TIMER
    timer_on("qrsp");
#endif

    /* q->p; r->q; s->r; p->s = spqr */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Gsp = Gs^Gp; Gqr = Gq^Gr;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  qr = InBuf->params->colidx[Q][R];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    sp = InBuf->params->rowidx[S][P];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gsp][sp][qr];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for qrsp sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }

#ifdef DPD_TIMER
    timer_off("qrsp");
#endif
    break;

  case qspr:
    fprintf(stderr,"\nDPD sort error: qspr index ordering not yet coded.\n");
    dpd_error("buf_sort", stderr);
    break;

  case qsrp:
    fprintf(stderr,"\nDPD sort error: qsrp index ordering not yet coded.\n");
    dpd_error("buf_sort", stderr);
    break;

  case rqps:

#ifdef DPD_TIMER
    timer_on("rqps");
#endif

    /* r->p; q->q; p->r; s->s = rqps */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Grq = Gr^Gq; Gps = Gp^Gs;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  rq = InBuf->params->rowidx[R][Q];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    ps = InBuf->params->colidx[P][S];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Grq][rq][ps];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for rqps sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }

#ifdef DPD_TIMER
    timer_off("rqps");
#endif
    break;

  case rqsp:

#ifdef DPD_TIMER
    timer_on("rqsp");
#endif

    /* r->p; q->q; s->r; p->s = sqpr */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Gsq = Gs^Gq;  Gpr = Gp^Gr;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  pr = InBuf->params->colidx[P][R];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    sq = InBuf->params->rowidx[S][Q];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gsq][sq][pr];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for rqsp sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }


#ifdef DPD_TIMER
    timer_off("rqsp");
#endif
    break;

  case rpqs:

#ifdef DPD_TIMER
    timer_on("rpqs");
#endif

    /* r->p; p->q; q->r; s->s = qrps */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Gqr = Gq^Gr;  Gps = Gp^Gs;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  qr = InBuf->params->rowidx[Q][R];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    ps = InBuf->params->colidx[P][S];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gqr][qr][ps];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {

      for(Gpq=0; Gpq < nirreps; Gpq++) {
	Grs = Gpq^my_irrep;

	/* determine how many rows of OutBuf we can store in half of the core */
	out_rows_per_bucket = dpd_memfree()/(2 * OutBuf.params->coltot[Grs]);
	if(out_rows_per_bucket > OutBuf.params->rowtot[Gpq])
	  out_rows_per_bucket = OutBuf.params->rowtot[Gpq];
	out_nbuckets = (int) ceil((double) OutBuf.params->rowtot[Gpq]/(double) out_rows_per_bucket);
	out_rows_left = OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

	for(n=0; n < (out_rows_left ? out_nbuckets-1 : out_nbuckets); n++) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	  for(Grow=0; Grow < nirreps; Grow++) {
	    Gcol = Grow^my_irrep;

	    /* determine how many rows of InBuf we can store in the other half of the core */
	    in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gcol]);
	    if(in_rows_per_bucket > InBuf->params->rowtot[Grow])
	      in_rows_per_bucket = InBuf->params->rowtot[Grow];
	    in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grow]/(double) in_rows_per_bucket);
	    in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

	    /* allocate space for the bucket of rows */
	    dpd_buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

	    for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

	      for(pq=0; pq < out_rows_per_bucket; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gqr = Gq^Gr;

		  if(Gqr == Grow) {
		    qr = InBuf->params->rowidx[q][r] - in_row_start;
		    /* check if the current value is in the current in_bucket or not */
		    if(qr >= 0 && qr < in_rows_per_bucket) {
		      ps = InBuf->params->colidx[p][s];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][qr][ps];
		    }
		  }
		}
	      }

	    }
	    if(in_rows_left) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

	      for(pq=0; pq < out_rows_per_bucket; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gqr = Gq^Gr;

		  if(Gqr == Grow) {
		    qr = InBuf->params->rowidx[q][r] - in_row_start;
		    /* check if the current value is in core or not */
		    if(qr >= 0 && qr < in_rows_left) {
		      ps = InBuf->params->colidx[p][s];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][qr][ps];
		    }
		  }
		}
	      }
	    }

	    dpd_buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	}
	if(out_rows_left) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_left);

	  for(Grow=0; Grow < nirreps; Grow++) {
	    Gcol = Grow^my_irrep;

	    /* determine how many rows of InBuf we can store in the other half of the core */
	    in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gcol]);
	    if(in_rows_per_bucket > InBuf->params->rowtot[Grow])
	      in_rows_per_bucket = InBuf->params->rowtot[Grow];
	    in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grow]/(double) in_rows_per_bucket);
	    in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

	    /* allocate space for the bucket of rows */
	    dpd_buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

	    for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

	      for(pq=0; pq < out_rows_left; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gqr = Gq^Gr;

		  if(Gqr == Grow) {
		    qr = InBuf->params->rowidx[q][r] - in_row_start;
		    /* check if the current value is in the current in_bucket or not */
		    if(qr >= 0 && qr < in_rows_per_bucket) {
		      ps = InBuf->params->colidx[p][s];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][qr][ps];
		    }
		  }
		}
	      }

	    }
	    if(in_rows_left) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

	      for(pq=0; pq < out_rows_left; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gqr = Gq^Gr;

		  if(Gqr == Grow) {
		    qr = InBuf->params->rowidx[q][r] - in_row_start;
		    /* check if the current value is in core or not */
		    if(qr >= 0 && qr < in_rows_left) {
		      ps = InBuf->params->colidx[p][s];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][qr][ps];
		    }
		  }
		}
	      }
	    }

	    dpd_buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_left);
	}

	dpd_buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);
      }
    }

#ifdef DPD_TIMER
    timer_off("rpqs");
#endif
    break;

  case rpsq:
	  
#ifdef DPD_TIMER
    timer_on("rpsq");
#endif

    /* r->p; p->q; s->r; q->s = qspr */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;

	    Gqs = Gq^Gs;  Gpr = Gp^Gr;

	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];

		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  pr = InBuf->params->colidx[P][R];

		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    qs = InBuf->params->rowidx[Q][S];

		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gqs][qs][pr];
			      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      for(Gpq=0; Gpq < nirreps; Gpq++) {
	Grs = Gpq^my_irrep;

	/* determine how many rows of OutBuf we can store in half of the core */
	out_rows_per_bucket = dpd_memfree()/(2 * OutBuf.params->coltot[Grs]);
	if(out_rows_per_bucket > OutBuf.params->rowtot[Gpq])
	  out_rows_per_bucket = OutBuf.params->rowtot[Gpq];
	out_nbuckets = (int) ceil((double) OutBuf.params->rowtot[Gpq]/(double) out_rows_per_bucket);
	out_rows_left = OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

	for(n=0; n < (out_rows_left ? out_nbuckets-1 : out_nbuckets); n++) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	  for(Grow=0; Grow < nirreps; Grow++) {
	    Gcol = Grow^my_irrep;

	    /* determine how many rows of InBuf we can store in the other half of the core */
	    in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gcol]);
	    if(in_rows_per_bucket > InBuf->params->rowtot[Grow])
	      in_rows_per_bucket = InBuf->params->rowtot[Grow];
	    in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grow]/(double) in_rows_per_bucket);
	    in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

	    /* allocate space for the bucket of rows */
	    dpd_buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

	    for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

	      for(pq=0; pq < out_rows_per_bucket; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gqs = Gq^Gs;

		  if(Gqs == Grow) {
		    qs = InBuf->params->rowidx[q][s] - in_row_start;
		    /* check if the current value is in the current in_bucket or not */
		    if(qs >= 0 && qs < in_rows_per_bucket) {
		      pr = InBuf->params->colidx[p][r];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][qs][pr];
		    }
		  }
		}
	      }

	    }
	    if(in_rows_left) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

	      for(pq=0; pq < out_rows_per_bucket; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gqs = Gq^Gs;

		  if(Gqs == Grow) {
		    qs = InBuf->params->rowidx[q][s] - in_row_start;
		    /* check if the current value is in core or not */
		    if(qs >= 0 && qs < in_rows_left) {
		      pr = InBuf->params->colidx[p][r];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][qs][pr];
		    }
		  }
		}
	      }
	    }

	    dpd_buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	}
	if(out_rows_left) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_left);

	  for(Grow=0; Grow < nirreps; Grow++) {
	    Gcol = Grow^my_irrep;

	    /* determine how many rows of InBuf we can store in the other half of the core */
	    in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gcol]);
	    if(in_rows_per_bucket > InBuf->params->rowtot[Grow])
	      in_rows_per_bucket = InBuf->params->rowtot[Grow];
	    in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grow]/(double) in_rows_per_bucket);
	    in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

	    /* allocate space for the bucket of rows */
	    dpd_buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

	    for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

	      for(pq=0; pq < out_rows_left; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gqs = Gq^Gs;

		  if(Gqs == Grow) {
		    qs = InBuf->params->rowidx[q][s] - in_row_start;
		    /* check if the current value is in the current in_bucket or not */
		    if(qs >= 0 && qs < in_rows_per_bucket) {
		      pr = InBuf->params->colidx[p][r];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][qs][pr];
		    }
		  }
		}
	      }

	    }
	    if(in_rows_left) {

	      in_row_start = m*in_rows_per_bucket;
	      dpd_buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

	      for(pq=0; pq < out_rows_left; pq++) {
		p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
		q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
		Gp = OutBuf.params->psym[p];
		Gq = Gpq^Gp;
		for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		  r = OutBuf.params->colorb[Grs][rs][0];
		  s = OutBuf.params->colorb[Grs][rs][1];
		  Gr = OutBuf.params->rsym[r];
		  Gs = Grs^Gr;

		  Gqs = Gq^Gs;

		  if(Gqs == Grow) {
		    qs = InBuf->params->rowidx[q][s] - in_row_start;
		    /* check if the current value is in core or not */
		    if(qs >= 0 && qs < in_rows_left) {
		      pr = InBuf->params->colidx[p][r];
		      OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grow][qs][pr];
		    }
		  }
		}
	      }
	    }

	    dpd_buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_left);
	}

	dpd_buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);
      }

    }

#ifdef DPD_TIMER
    timer_off("rpsq");
#endif
    break;

  case rsqp:

#ifdef DPD_TIMER
    timer_on("rsqp");
#endif

    /* r->p; s->q; q->r; p->s = srpq */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	  p = OutBuf.params->roworb[h][pq][0];
	  q = OutBuf.params->roworb[h][pq][1];

	  col = InBuf->params->colidx[p][q];
	  
	  for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	    r = OutBuf.params->colorb[r_irrep][rs][0];
	    s = OutBuf.params->colorb[r_irrep][rs][1];

	    row = InBuf->params->rowidx[s][r];
		  
	    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[h][row][col];

	  }
	}
      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for rsqp sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }

#ifdef DPD_TIMER
    timer_off("rsqp");
#endif
    break;

  case rspq:

#ifdef DPD_TIMER
    timer_on("rspq");
#endif

    /* r->p; s->q; p->r; q->s = rspq */

    /* loop over row irreps of OutBuf */
    if(incore) {
#ifdef DPD_TIMER
      timer_on("axpy:sort");
#endif
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	  p = OutBuf.params->roworb[h][pq][0];
	  q = OutBuf.params->roworb[h][pq][1];

	  col = InBuf->params->colidx[p][q];
	  
	  for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	    r = OutBuf.params->colorb[r_irrep][rs][0];
	    s = OutBuf.params->colorb[r_irrep][rs][1];

	    row = InBuf->params->rowidx[r][s];
  
	    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[r_irrep][row][col];
	  }
	}
      }
#ifdef DPD_TIMER
      timer_off("axpy:sort");
#endif
    }
    else {

      for(Gpq=0; Gpq < nirreps; Gpq++) {
	Grs = Gpq^my_irrep;

	/* determine how many rows of OutBuf we can store in half of the core */
	out_rows_per_bucket = dpd_memfree()/(2 * OutBuf.params->coltot[Grs]);
	if(out_rows_per_bucket > OutBuf.params->rowtot[Gpq])
	  out_rows_per_bucket = OutBuf.params->rowtot[Gpq];
	out_nbuckets = (int) ceil((double) OutBuf.params->rowtot[Gpq]/(double) out_rows_per_bucket);
	out_rows_left = OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

	/* determine how many rows of InBuf we can store in the other half of the core */
	in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gpq]);
	if(in_rows_per_bucket > InBuf->params->rowtot[Grs])
	  in_rows_per_bucket = InBuf->params->rowtot[Grs];
	in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grs]/(double) in_rows_per_bucket);
	in_rows_left = InBuf->params->rowtot[Grs] % in_rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(InBuf, Grs, in_rows_per_bucket);

	for(n=0; n < (out_rows_left ? out_nbuckets-1 : out_nbuckets); n++) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	  for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	    in_row_start = m*in_rows_per_bucket;
	    dpd_buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start, in_rows_per_bucket);

	    for(pq=0; pq < out_rows_per_bucket; pq++) {
	      for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		if(rs >= in_row_start && rs < in_rows_per_bucket+in_row_start)
		  OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grs][rs-in_row_start][pq+out_row_start];
	      }
	    }

	  }
	  if(in_rows_left) {

	    in_row_start = m*in_rows_per_bucket;
	    dpd_buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start, in_rows_left);

	    for(pq=0; pq < out_rows_per_bucket; pq++) {
	      for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		if(rs >= in_row_start && rs < in_rows_left+in_row_start)
		  OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grs][rs-in_row_start][pq+out_row_start];
	      }
	    }
	  }

	  /*	  mat_print(OutBuf.matrix[Gpq], out_rows_per_bucket, OutBuf.params->coltot[Grs], stderr); */
	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	}
	if(out_rows_left) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_left);

	  for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	    in_row_start = m*in_rows_per_bucket;
	    dpd_buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start, in_rows_per_bucket);

	    for(pq=0; pq < out_rows_left; pq++) {
	      for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		if(rs >= in_row_start && rs < in_rows_per_bucket+in_row_start)
		  OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grs][rs-in_row_start][pq+out_row_start];
	      }
	    }

	  }
	  if(in_rows_left) {

	    in_row_start = m*in_rows_per_bucket;
	    dpd_buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start, in_rows_left);

	    for(pq=0; pq < out_rows_left; pq++) {
	      for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		if(rs >= in_row_start && rs < in_rows_left+in_row_start)
		  OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grs][rs-in_row_start][pq+out_row_start];
	      }
	    }
	  }

	  /*	  mat_print(OutBuf.matrix[Gpq], out_rows_left, OutBuf.params->coltot[Grs], stderr); */
	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_left);
	}

	dpd_buf4_mat_irrep_close_block(InBuf, Grs, in_rows_per_bucket);

	dpd_buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);

      } /* Gpq */
    }

#ifdef DPD_TIMER
    timer_off("rspq");
#endif
    break;

  case sqrp:

#ifdef DPD_TIMER
    timer_on("sqrp");
#endif

    /* s->p; q->q; r->r; p->s = sqrp */

    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	    Gs = Gr^r_irrep;
		  
	    Gsq = Gs^Gq;  Grp = Gr^Gp;
		  
	    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	      P = OutBuf.params->poff[Gp] + p;
	      for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
		Q = OutBuf.params->qoff[Gq] + q;
		pq = OutBuf.params->rowidx[P][Q];
			  
		for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		  R = OutBuf.params->roff[Gr] + r;
		  rp = InBuf->params->colidx[R][P];
			      
		  for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		    S = OutBuf.params->soff[Gs] + s;
		    rs = OutBuf.params->colidx[R][S];
		    sq = InBuf->params->rowidx[S][Q];
				  
		    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gsq][sq][rp];
		      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for sqrp sort.\n");
      dpd_error("buf4_sort_axpy", stderr);
    }
#ifdef DPD_TIMER
    timer_off("sqrp");
#endif
    break;

  case sqpr:
    fprintf(stderr,"\nDPD sort error: sqpr index ordering not yet coded.\n");
    dpd_error("buf_sort", stderr);
    break;

  case srqp:

#ifdef DPD_TIMER
    timer_on("srqp");
#endif

    /* s->p; r->q; q->r; p->s = srqp */
    if(incore) {
      for(h=0; h < nirreps; h++) {
	r_irrep = h^my_irrep;
	  
	for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	  p = OutBuf.params->roworb[h][pq][0];
	  q = OutBuf.params->roworb[h][pq][1];

	  col = InBuf->params->colidx[q][p];
	  
	  for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	    r = OutBuf.params->colorb[r_irrep][rs][0];
	    s = OutBuf.params->colorb[r_irrep][rs][1];

	    row = InBuf->params->rowidx[s][r];
		  
	    OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[r_irrep][row][col];

	  }
	}
      }
    }
    else {

      for(Gpq=0; Gpq < nirreps; Gpq++) {
	Grs = Gpq^my_irrep;

	/* determine how many rows of OutBuf we can store in half of the core */
	out_rows_per_bucket = dpd_memfree()/(2 * OutBuf.params->coltot[Grs]);
	if(out_rows_per_bucket > OutBuf.params->rowtot[Gpq])
	  out_rows_per_bucket = OutBuf.params->rowtot[Gpq];
	out_nbuckets = (int) ceil((double) OutBuf.params->rowtot[Gpq]/(double) out_rows_per_bucket);
	out_rows_left = OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

	/* determine how many rows of InBuf we can store in the other half of the core */
	in_rows_per_bucket = dpd_memfree()/(2 * InBuf->params->coltot[Gpq]);
	if(in_rows_per_bucket > InBuf->params->rowtot[Grs])
	  in_rows_per_bucket = InBuf->params->rowtot[Grs];
	in_nbuckets = (int) ceil((double) InBuf->params->rowtot[Grs]/(double) in_rows_per_bucket);
	in_rows_left = InBuf->params->rowtot[Grs] % in_rows_per_bucket;

	/* allocate space for the bucket of rows */
	dpd_buf4_mat_irrep_init_block(InBuf, Grs, in_rows_per_bucket);

	for(n=0; n < (out_rows_left ? out_nbuckets-1 : out_nbuckets); n++) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	  for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	    in_row_start = m*in_rows_per_bucket;
	    dpd_buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start, in_rows_per_bucket);

	    for(pq=0; pq < out_rows_per_bucket; pq++) {
	      p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
	      q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
	      qp = InBuf->params->colidx[q][p];

	      for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		r = OutBuf.params->colorb[Grs][rs][0];
		s = OutBuf.params->colorb[Grs][rs][1];

		sr = InBuf->params->rowidx[s][r] - in_row_start;

		if(sr >= 0 && sr < in_rows_per_bucket)
		  OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grs][sr][qp];
	      }
	    }

	  }
	  if(in_rows_left) {

	    in_row_start = m*in_rows_per_bucket;
	    dpd_buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start, in_rows_left);

	    for(pq=0; pq < out_rows_per_bucket; pq++) {
	      p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
	      q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
	      qp = InBuf->params->colidx[q][p];

	      for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		r = OutBuf.params->colorb[Grs][rs][0];
		s = OutBuf.params->colorb[Grs][rs][1];

		sr = InBuf->params->rowidx[s][r] - in_row_start;

		if(sr >= 0 && sr < in_rows_left)
		  OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grs][sr][qp];
	      }
	    }
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);

	}
	if(out_rows_left) {

	  out_row_start = n*out_rows_per_bucket;
	  dpd_buf4_mat_irrep_rd_block(&OutBuf, Gpq, out_row_start, out_rows_left);

	  for(m=0; m < (in_rows_left ? in_nbuckets-1 : in_nbuckets); m++) {

	    in_row_start = m*in_rows_per_bucket;
	    dpd_buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start, in_rows_per_bucket);

	    for(pq=0; pq < out_rows_left; pq++) {
	      p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
	      q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
	      qp = InBuf->params->colidx[q][p];

	      for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		r = OutBuf.params->colorb[Grs][rs][0];
		s = OutBuf.params->colorb[Grs][rs][1];

		sr = InBuf->params->rowidx[s][r] - in_row_start;

		if(sr >= 0 && sr < in_rows_per_bucket)
		  OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grs][sr][qp];
	      }
	    }

	  }
	  if(in_rows_left) {

	    in_row_start = m*in_rows_per_bucket;
	    dpd_buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start, in_rows_left);

	    for(pq=0; pq < out_rows_left; pq++) {
	      p = OutBuf.params->roworb[Gpq][pq+out_row_start][0];
	      q = OutBuf.params->roworb[Gpq][pq+out_row_start][1];
	      qp = InBuf->params->colidx[q][p];

	      for(rs=0; rs < OutBuf.params->coltot[Grs]; rs++) {
		r = OutBuf.params->colorb[Grs][rs][0];
		s = OutBuf.params->colorb[Grs][rs][1];

		sr = InBuf->params->rowidx[s][r] - in_row_start;

		if(sr >= 0 && sr < in_rows_left)
		  OutBuf.matrix[Gpq][pq][rs] += alpha * InBuf->matrix[Grs][sr][qp];
	      }
	    }
	  }

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_left);
	}

	dpd_buf4_mat_irrep_close_block(InBuf, Grs, in_rows_per_bucket);

	dpd_buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);

      } /* Gpq */

    }

#ifdef DPD_TIMER
  timer_off("srqp");
#endif
  break;

 case srpq:
#ifdef DPD_TIMER
   timer_on("srpq");
#endif

   /* s->p; r->q; p->r; q->s = rsqp */
   if(incore) {
     for(h=0; h < nirreps; h++) {
       r_irrep = h^my_irrep;

       for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	 p = OutBuf.params->roworb[h][pq][0];
	 q = OutBuf.params->roworb[h][pq][1];

	 col = InBuf->params->colidx[q][p];

	 for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	   r = OutBuf.params->colorb[r_irrep][rs][0];
	   s = OutBuf.params->colorb[r_irrep][rs][1];

	   row = InBuf->params->rowidx[r][s];

	   OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[h][row][col];

	 }
       }
     }
   }
   else {
     fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for srpq sort.\n");
     dpd_error("buf4_sort_axpy", stderr);
   }
#ifdef DPD_TIMER
   timer_off("srpq");
#endif
   break;

 case spqr:

#ifdef DPD_TIMER
   timer_on("spqr");
#endif

   /* s->p; p->q; q->r; r->s = qrsp */
   if(incore) {
     for(h=0; h < nirreps; h++) {
       r_irrep = h^my_irrep;
	  
       for(Gp=0; Gp < nirreps; Gp++) {
	 Gq = Gp^h;
	 for(Gr=0; Gr < nirreps; Gr++) {
	   Gs = Gr^r_irrep;
		  
	   Gqr = Gq^Gr;  Gsp = Gs^Gp;
		  
	   for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	     P = OutBuf.params->poff[Gp] + p;
	     for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	       Q = OutBuf.params->qoff[Gq] + q;
	       pq = OutBuf.params->rowidx[P][Q];
			  
	       for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		 R = OutBuf.params->roff[Gr] + r;
		 qr = InBuf->params->rowidx[Q][R];
			      			      
		 for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		   S = OutBuf.params->soff[Gs] + s;
		   rs = OutBuf.params->colidx[R][S];
		   sp = InBuf->params->colidx[S][P];
  
		   OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gqr][qr][sp];
		      
		 }
	       }
	     }
	   }
	 }
       }
     }
   }
   else {
     fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for spqr sort.\n");
     dpd_error("buf4_sort_axpy", stderr);
   }
#ifdef DPD_TIMER
   timer_off("spqr");
#endif
   break;

 case sprq:

#ifdef DPD_TIMER
   timer_on("sprq");
#endif

   /* s->p; p->q; r->r; q->s = qsrp */
   if(incore) {
     for(h=0; h < nirreps; h++) {
       r_irrep = h^my_irrep;
	  
       for(Gp=0; Gp < nirreps; Gp++) {
	 Gq = Gp^h;
	 for(Gr=0; Gr < nirreps; Gr++) {
	   Gs = Gr^r_irrep;
		  
	   Gqs = Gq^Gs;  Grp = Gr^Gp;
		  
	   for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	     P = OutBuf.params->poff[Gp] + p;
	     for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	       Q = OutBuf.params->qoff[Gq] + q;
	       pq = OutBuf.params->rowidx[P][Q];
			  
	       for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		 R = OutBuf.params->roff[Gr] + r;
		 rp = InBuf->params->colidx[R][P];
			      
		 for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		   S = OutBuf.params->soff[Gs] + s;
		   rs = OutBuf.params->colidx[R][S];
		   qs = InBuf->params->rowidx[Q][S];
				  
		   OutBuf.matrix[h][pq][rs] += alpha * InBuf->matrix[Gqs][qs][rp];
		      
		 }
	       }
	     }
	   }
	 }
       }
     }
   }
   else {
     fprintf(stderr, "LIBDPD: Out-of-core algorithm not yet coded for sprq sort.\n");
     dpd_error("buf4_sort_axpy", stderr);
   }
#ifdef DPD_TIMER
   timer_off("sprq");
#endif
   break;
}

#ifdef DPD_TIMER
  if(index == rspq) timer_on("axpy:alloc");
#endif
  if(incore) {
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_wrt(&OutBuf, h);
      dpd_buf4_mat_irrep_close(&OutBuf, h);
      dpd_buf4_mat_irrep_close(InBuf, h);
    }
  }
#ifdef DPD_TIMER
  if(index == rspq) timer_off("axpy:alloc");
#endif

  dpd_buf4_close(&OutBuf);

#ifdef DPD_TIMER
  timer_off("buf4_sort");
#endif
}

}
