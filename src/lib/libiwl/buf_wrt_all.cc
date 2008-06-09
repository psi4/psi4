/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include "iwl.h"
#include "iwl.hpp"

namespace psi {
  
void IWL::write_all(int nbfso, double *ints, 
    int *ioff, int printflg, FILE *outfile)
{
    int idx, p, q, r, s, smax, pq, rs, pqrs;
    Label *lblptr;
    Value *valptr;

    lblptr = labels_;
    valptr = values_;

    /* go through the lexical order and print to the output file */
    for (p=0; p<nbfso; p++) {
        for (q=0; q<=p; q++) {
            pq = ioff[p] + q;
            for (r=0; r<=p; r++) {
                smax = (p==r) ? (q+1) : (r+1);
                for (s=0; s < smax; s++) {
                    rs = ioff[r] + s;
                    pqrs = ioff[pq] + rs;
                    if (fabs(ints[pqrs]) > cutoff_) {
                        idx = 4 * idx_;
                        lblptr[idx++] = (Label) p;
                        lblptr[idx++] = (Label) q;
                        lblptr[idx++] = (Label) r;
                        lblptr[idx++] = (Label) s;
                        valptr[idx_] = (Value) ints[pqrs];
                        idx_++;
                        if (printflg) fprintf(outfile, "%d %d %d %d [%d] = %10.6f\n",
                            p, q, r, s, pqrs, ints[pqrs]) ;

                        if (idx_ == ints_per_buf_) {
                            lastbuf_ = 0;
                            inbuf_ = idx_;
                            put();
                            idx_ = 0;
                        } 
                    }
                }
            }
        }
    }
}

/*!
** iwl_buf_wrt_all()
**
** Write out two electron ints to IWL file.  Assume that the integrals
** are in ijkl canonical order (no spatial symmetry).
**
**    \param itap     = unit to write to
**    \param nbfso    = number of basis functions in symmetry orbitals
**    \param ints     = two electron integrals 
**    \param ioff     = the old ioff array for lexical ordering
**    \param printflg = print flag (1 or 0)
**    \param outfile  =  output file
**
** David Sherrill, 6/27/96
**
** NB: This routine will only write "standard" (pq|rs) indices to disk.
** The cints integral program marks certain integrals with negative index
** values to indicate the end of PK-matrix blocks.  This marking is used
** by the cscf code (only?).  Therefore, is this routine is used to write
** integrals to disk, the SCF code will most likely give incorrect data
** if it uses the resulting integrals.
** TDC 12/24/01
** \ingroup IWL
*/
void iwl_buf_wrt_all(struct iwlbuf *Buf, int nbfso, double *ints, int *ioff,
      int printflg, FILE *outfile)
{
  int idx, p, q, r, s, smax, pq, rs, pqrs;
  Label *lblptr;
  Value *valptr;

  lblptr = Buf->labels;
  valptr = Buf->values;
  
  /* go through the lexical order and print to the output file */
  for (p=0; p<nbfso; p++) {
    for (q=0; q<=p; q++) {
      pq = ioff[p] + q;
      for (r=0; r<=p; r++) {
	smax = (p==r) ? (q+1) : (r+1);
	for (s=0; s < smax; s++) {
	  rs = ioff[r] + s;
	  pqrs = ioff[pq] + rs;
	  if (fabs(ints[pqrs]) > Buf->cutoff) {
	    idx = 4 * Buf->idx;
	    lblptr[idx++] = (Label) p;
	    lblptr[idx++] = (Label) q;
	    lblptr[idx++] = (Label) r;
	    lblptr[idx++] = (Label) s;
	    valptr[Buf->idx] = (Value) ints[pqrs];
	    Buf->idx++;
	    if (printflg) fprintf(outfile, "%d %d %d %d [%d] = %10.6lf\n",
				  p, q, r, s, pqrs, ints[pqrs]) ;
	    
	    if (Buf->idx == Buf->ints_per_buf) {
	      Buf->lastbuf = 0;
	      Buf->inbuf = Buf->idx;
	      iwl_buf_put(Buf);
	      Buf->idx = 0;
	    } 
	  }
	}
      }
    }
  }
}

}

