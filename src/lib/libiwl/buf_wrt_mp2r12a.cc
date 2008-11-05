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

void IWL::write_mp2r12a(int p, int q, int pq, int pqsym, double **arr, 
    int rsym, int *firstr, int *lastr, int *firsts, int *lasts, 
    int *occ, int bra_ket_symm, int *ioff, int printflag, FILE *outfile)
{
    int idx, r, s, rs, ssym;
    int R,S,rnew,snew;
    double value;
    Label *lblptr;
    Value *valptr;

    lblptr = labels_;
    valptr = values_;

    ssym = pqsym ^ rsym;
    for (r=firstr[rsym],R=0; r <= lastr[rsym]; r++,R++) {
        rnew = occ[r];  /* r-index is in QTS-ordering, not Pitzer */
        for (s=firsts[ssym],S=0; s <=lasts[ssym]; s++,S++) {
            snew = s;     /* s-index is in Pitzer ordering */
            rs = ioff[rnew] + snew;
            /*---------------------------------------
            If bra_ket_symm != 0 -> we do not need
            integrals with rs > pq. rs can only
            increase here
            ---------------------------------------*/
            if (bra_ket_symm && rs > pq)
                return;
            value = arr[R][S];

            if (fabs(value) > cutoff_) {
                idx = 4 * idx_;
                lblptr[idx++] = (Label) p;
                lblptr[idx++] = (Label) q;
                lblptr[idx++] = (Label) rnew;
                lblptr[idx++] = (Label) snew;
                valptr[idx_] = (Value) value;

                idx_++;

                if (idx_ == ints_per_buf_) {
                    lastbuf_ = 0;
                    inbuf_ = idx_;
                    put();
                    idx_ = 0;
                }

                if(printflag)
                    fprintf(outfile, "<%d %d %d %d [%d] [%d] = %20.10f\n",
                    p, q, rnew, snew, pq, rs, value);

            } /* end if (fabs(value) > Buf->cutoff) ... */
        } /* end loop over s */
    } /* end loop over r */
}

/*!
** iwl_buf_wrt_mp2r12a()
**
** Write to an Integrals With Labels formatted PSI buffer.
** The buffer must have been initialized with iwl_buf_init().  Don't
** forget to call iwl_buf_flush() when finished with all writes to the
** buffer to ensure that all contents are written to disk.
** David Sherrill, March, 1995
**
** This routine is a modified form of iwl_buf_wrt() specific to mp2-type
** restricted transforms.  It's not general, but it should work.
** Daniel, 9/25/95
**
** This routine is a modified form of iwl_buf_wrt_mp2() specific to mp2r12a-type
** restricted transforms.
** Edward, 8/04/99
** \ingroup IWL
*/
void iwl_buf_wrt_mp2r12a(struct iwlbuf *Buf, int p, int q, int pq, int pqsym,
   double **arr, int rsym, int *firstr, int *lastr, int *firsts, int *lasts,
   int *occ, int bra_ket_symm, int *ioff, int printflag, FILE *outfile)
{
   int idx, r, s, rs, ssym;
   int R,S,rnew,snew;
   double value;
   Label *lblptr;
   Value *valptr;

   lblptr = Buf->labels;
   valptr = Buf->values;

   ssym = pqsym ^ rsym;
   for (r=firstr[rsym],R=0; r <= lastr[rsym]; r++,R++) {
     rnew = occ[r];  /* r-index is in QTS-ordering, not Pitzer */
     for (s=firsts[ssym],S=0; s <=lasts[ssym]; s++,S++) {
       snew = s;     /* s-index is in Pitzer ordering */
       rs = ioff[rnew] + snew;
       /*---------------------------------------
	 If bra_ket_symm != 0 -> we do not need
	 integrals with rs > pq. rs can only
	 increase here
	---------------------------------------*/
       if (bra_ket_symm && rs > pq)
	 return;
       value = arr[R][S];
       
       if (fabs(value) > Buf->cutoff) {
	 idx = 4 * Buf->idx;
	 lblptr[idx++] = (Label) p;
	 lblptr[idx++] = (Label) q;
	 lblptr[idx++] = (Label) rnew;
	 lblptr[idx++] = (Label) snew;
	 valptr[Buf->idx] = (Value) value;
	 
	 Buf->idx++;
	 
	 if (Buf->idx == Buf->ints_per_buf) {
	   Buf->lastbuf = 0;
	   Buf->inbuf = Buf->idx;
	   iwl_buf_put(Buf);
	   Buf->idx = 0;
	 }
	 
	 if(printflag)
	   fprintf(outfile, "<%d %d %d %d [%d] [%d] = %20.10f\n",
		   p, q, rnew, snew, pq, rs, value);
	 
       } /* end if (fabs(value) > Buf->cutoff) ... */
     } /* end loop over s */
   } /* end loop over r */
   
}

}

