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
  
void IWL::write(int p, int q, int pq, int pqsym,
    double *arr, int rmax, int *active, int *ioff, int *orbsym, int *firsti, 
    int *lasti, int sortby_rs, int printflag, FILE *outfile)
{
    int r, s, rs, rsym, ssym, smax, idx;
    double value;
    Label *lblptr;
    Value *valptr;

    lblptr = labels_;
    valptr = values_;

    for (r=0; r<rmax; r++) {
        rsym = orbsym[r];
        ssym = pqsym ^ rsym;
        smax = (rsym == ssym) ? r : lasti[ssym];

        for (s=firsti[ssym]; s<=smax; s++) {
            rs = ioff[r] + s;
            value = arr[rs];

            if (fabs(value) > cutoff_) {
                idx = 4 * idx_;
                lblptr[idx] = (Label) p;
                lblptr[idx+1] = (Label) q;
                lblptr[idx+2] = (Label) r;
                lblptr[idx+3] = (Label) s;
                valptr[idx_] = (Value) value;

                idx_++;

                if (idx_ == ints_per_buf_) {
                    inbuf_ = idx_;
                    lastbuf_ = 0;
                    put();
                    idx_ = 0;
                }

                if(printflag)
                    fprintf(outfile, "<%d %d %d %d [%d] [%d] = %20.10f\n",
                    p, q, r, s, pq, rs, value);

            } /* end if (fabs(value) > Buf->cutoff) ... */
        } /* end loop over s */
    } /* end loop over r */
}

/*!
** iwl_buf_wrt()
**
** Write to an Integrals With Labels formatted buffer.
** The buffer must have been initialized with iwl_buf_init().  Don't
** forget to call iwl_buf_flush() when finished with all writes to the
** buffer to ensure that all contents are written to disk.
** David Sherrill, March 1995
**
** Revised 6/27/96 by CDS for new format
** \ingroup IWL
*/
void iwl_buf_wrt(struct iwlbuf *Buf, int p, int q, int pq, int pqsym,
   double *arr, int rmax, int *active, int *ioff, int *orbsym, int *firsti, 
   int *lasti, int sortby_rs, int printflag, FILE *outfile)
{
  int r, s, rs, rsym, ssym, smax, idx;
  double value;
  Label *lblptr;
  Value *valptr;
  
  lblptr = Buf->labels;
  valptr = Buf->values;
  
  for (r=0; r<rmax; r++) {
    rsym = orbsym[r];
    ssym = pqsym ^ rsym;
    smax = (rsym == ssym) ? r : lasti[ssym];
    
    for (s=firsti[ssym]; s<=smax; s++) {
      rs = ioff[r] + s;
      value = arr[rs];
      
      if (fabs(value) > Buf->cutoff) {
	idx = 4 * Buf->idx;
	lblptr[idx] = (Label) p;
	lblptr[idx+1] = (Label) q;
	lblptr[idx+2] = (Label) r;
	lblptr[idx+3] = (Label) s;
	valptr[Buf->idx] = (Value) value;
	 
	Buf->idx++;
	
	if (Buf->idx == Buf->ints_per_buf) {
	  Buf->inbuf = Buf->idx;
	  Buf->lastbuf = 0;
	  iwl_buf_put(Buf);
	  Buf->idx = 0;
	}
	
	if(printflag)
	  fprintf(outfile, "<%d %d %d %d [%d] [%d] = %20.10lf\n",
		  p, q, r, s, pq, rs, value);
	 
      } /* end if (fabs(value) > Buf->cutoff) ... */
    } /* end loop over s */
  } /* end loop over r */
  
}

}

