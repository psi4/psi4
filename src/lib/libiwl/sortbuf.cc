/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include "iwl.h"
#include "iwl.hpp"

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

#define BIGNUM 40000

namespace psi {

void IWL::sort_buffer(IWL *Inbuf, IWL *Outbuf,
    double *ints, int fpq, int lpq, int *ioff, int *ioff2, 
    int nbfso, int elbert, int intermediate, int no_pq_perm, 
    int qdim, int add, int printflg, FILE *outfile)
{
    int i;
    Value *valptr;              /* array of integral values */
    Label *lblptr;              /* array of integral labels */
    int idx;                    /* index for curr integral (0..ints_per_buf) */
    int lastbuf;                /* last buffer flag */
    int p, q, qmax, qmin, r, rmin, rmax, s, smin, smax, pq, rs;
    long int pqrs, offset;
    int first_p, first_q, first_pq, last_p, last_q;
    int nbstri;

    if (printflg) {
        fprintf(outfile, "\nsortbuf for pq=%d to %d\n", fpq, lpq);
    }

    if (no_pq_perm && !intermediate) {
        fprintf(outfile,"(sortbuf): illegal parameter combination.\n");
        fprintf(stderr, "(sortbuf): illegal parameter combination.\n");
    }

    nbstri = nbfso * (nbfso + 1) / 2;

    /* figure out ranges on things */
    /* I believe this section works fine, even with different ioff arrays */
    i = 0;
    while (fpq >= ioff[i] && i < BIGNUM) i++;
    if (i == BIGNUM) {
        fprintf(outfile, "(sortbuf): parameter error\n") ;
        return;
    }
    first_p = i-1 ; first_q = fpq - ioff[i-1];
    first_pq = ioff[first_p] + first_q;
    if (first_pq != fpq) {
        fprintf(outfile, "(sortbuf): fpq != first_pq.\n");
        fprintf(stderr,  "(sortbuf): fpq != first_pq.\n");
    }

    if (!intermediate) {
        if (elbert) offset = ioff2[first_pq] + first_pq ;
        else offset = ioff[first_pq] ;
    }
    else offset = 0;

    i=0; 
    while (lpq >= ioff[i] && i < BIGNUM) i++ ;
    if (i == BIGNUM) {
        fprintf(outfile, "(sortbuf): parameter error\n") ;
        return ;
    }
    last_p = i-1 ; last_q = lpq - ioff[i-1] ;


    lblptr = Inbuf->labels_;
    valptr = Inbuf->values_;

    /* read a buffer at a time until we're done */

    do {
        Inbuf->fetch();
        lastbuf = Inbuf->lastbuf_;
        for (idx=4*Inbuf->idx_; Inbuf->idx_<Inbuf->inbuf_; Inbuf->idx_++) {
            p = (int) lblptr[idx++];
            q = (int) lblptr[idx++];
            r = (int) lblptr[idx++];
            s = (int) lblptr[idx++];

            /* if (no_pq_perm) ioff is the appropriate offset array for the left
               indices (ioff[p] = nvirt * p for MP2); ioff2 is then the usual
               ioff offset array, used for the right indices */
            if (no_pq_perm) {
                pq = ioff[p] + q;
                rs = ioff2[MAX0(r,s)] + MIN0(r,s);
            }
            else {
                pq = ioff[MAX0(p,q)] + MIN0(p,q);
                rs = ioff[MAX0(r,s)] + MIN0(r,s);
            }

            if (!intermediate) {
                if (elbert)
                    pqrs = ioff2[pq] + rs;
                else {
                    pqrs = ioff[MAX0(pq,rs)];
                    pqrs += MIN0(pq,rs);
                }
            }
            else {
                pqrs = (pq - first_pq);
                pqrs *= nbstri;
                pqrs += rs;
            }

            if (printflg && ints[pqrs-offset] != 0.0) 
                fprintf(outfile, "Adding %10.6f to el %d %d %d %d = %10.6lf\n", 
                valptr[Inbuf->idx_], p, q, r, s, ints[pqrs-offset]);

            if (add) ints[pqrs-offset] += valptr[Inbuf->idx_];
            else ints[pqrs-offset] += valptr[Inbuf->idx_];

            if (printflg) 
                fprintf(outfile, "<%d %d %d %d | %d %d [%ld] = %10.6f\n",
                p, q, r, s, pq, rs, pqrs, ints[pqrs-offset]) ;
        }
    } while (!lastbuf);

    /* now write them out again, in order */
    lblptr = Outbuf->labels_;
    valptr = Outbuf->values_;

    idx = 0;

    for (p=first_p; p<=last_p; p++) {
        qmax = (p==last_p) ? last_q : p ;
        qmin = (p==first_p) ? first_q : 0 ;
        if(no_pq_perm) {
            qmax = (p==last_p) ? last_q : (qdim-1);
            qmin = (p==first_p) ? first_q: 0;
        }
        for (q=qmin; q<=qmax; q++) {
            pq = ioff[p] + q ; /* This should be fine even with MP2 */

            if (!intermediate) {
                rmin = (elbert) ? p : 0 ; 
                rmax = (elbert) ? nbfso : p+1 ;
            }
            else {  /* This should be fine with MP2, also */
                rmin = 0;
                rmax = nbfso;
            }

            for (r=rmin; r<rmax; r++) {

                if (!intermediate) {
                    if (elbert) {
                        smax = r+1 ;
                        smin = (p==r) ? q : 0 ;
                    }
                    else {
                        smax = (p==r) ? (q+1) : (r+1) ;
                        smin = 0 ;
                    }
                }
                else { /* This should be fine with MP2, also */
                    smax = r + 1;
                    smin = 0;
                }

                for (s=smin; s < smax; s++) {
                    if(no_pq_perm) rs = ioff2[r] + s;
                    else rs = ioff[r] + s;

                    /* Again, this should be fine with MP2 */
                    if (elbert) pqrs = ioff2[pq] + rs ;
                    else if (intermediate) {
                        pqrs = (pq - first_pq);
                        pqrs *= nbstri;
                        pqrs += rs;
                    }
                    else pqrs = ioff[pq] + rs ;

                    if (fabs(ints[pqrs-offset]) > Outbuf->cutoff_) {
                        idx = 4*Outbuf->idx_;
                        lblptr[idx++] = p;
                        lblptr[idx++] = q;
                        lblptr[idx++] = r;
                        lblptr[idx++] = s;
                        valptr[Outbuf->idx_] = ints[pqrs-offset];
                        if (printflg) 
                            fprintf(outfile, ">%d %d %d %d | %d %d [%ld] = %10.6f\n",
                            p, q, r, s, pq, rs, pqrs, ints[pqrs-offset]) ;

                        Outbuf->idx_++;
                        if (Outbuf->idx_ == Outbuf->ints_per_buf_) {
                            Outbuf->lastbuf_ = 0;
                            Outbuf->inbuf_ = Outbuf->idx_;
                            Outbuf->put();
                            Outbuf->idx_ = 0;
                        } 
                    }
                }
            }
        }
    }
}

/*!
** sortbuf()
**
** Function reads a file of two-electron integrals into
** core and writes them back out again in canonical order.  Used in
** Yoshimine sorts where we have a file containing all rs for a few
** values of pq, but the ints are not in canonical order.  At the 
** very least, we need to sort to make sure that all (pq|rs) for a
** given pq are grouped together, since the transformation program
** wants to work with all rs for a single pq value at one time.
** We may or may not use the restriction pq >= rs (not used if 
** intermediate = 1, which is how this routine is always called
** right now).
**
** One interesting issue here is that the intermediate array ('ints')
** must be big enough to hold the integrals in the current buffer, but
** we don't generally want it to be much larger than necessary!  Thus
** we calculate an 'offset' which is the canonical index of the first
** integral in the buffer, and we use this so that the first integral
** in the buffer is stored in ints[0].  What's different for the
** Elbert ordering is that we have all rs for a given pq but rs >= pq!
** If we want our canonical indices to be consecutive WITHIN THE
** CURRENT BUFFER, we MUST use upper triangle ordering rather than 
** lower triangle!  That's what ioff2 is used for.  Obviously the
** offset must also be calculated with ioff2 for the Elbert order.
** Formula for ioff2: ioff2[0] = 0; ioff2[i] = ioff2[i-1] + n - i;
** Note that this is not the case when this routine is used for MP2
** sorts.  The definitions of the ioff arrays can be confusing, so
** care should be taken when using this routine.
**
**    \param Inbuf       = IWL buffer for input
**    \param Outbuf      = IWL buffer for output
**    \param ints        = array to hold integrals in
**    \param fpq         = first pq for this tape
**    \param lpq         = last pq for this tape 
**    \param ioff        = offset array for the left indices
**    \param ioff2       = offset array for Elbert sorts or for the right 
**                  indices when no_pq_perm=1 
**    \param nbfso       = number of basis functions in SO's
**    \param lastsort    = 1 if this is the last intape, 0 otherwise
**    \param elbert      = integrals obey rs >= pq.  Use ioff2 to get offset.
**    \param intermediate= 1 if sorting a intermediate in the transformation
**                  which is indexed as X[ij][kl] where ij runs from
**                  fpq to lpq and kl runs from 0 to nbstri
**    \param no_pq_perm  = don't use permutational symmetry to swap p and q
**                  (appropriate for MP2 where one is occ and one is virt)
**    \param qdim        = dimensions for the q index...nvirt for MP2
**    \param add         = add contributions to the same integral during sort
**    \param printflg    = 1 for printing, 0 otherwise
**    \param outfile     = output file pointer
**
** Returns: none
**
** Revised 6/27/96 by CDS for new IWL format
** N.B. Now need to iwl_flush the output buffer...not done in here!!
** \ingroup IWL
*/
void sortbuf(struct iwlbuf *Inbuf, struct iwlbuf *Outbuf,
      double *ints, int fpq, int lpq, int *ioff, int *ioff2, 
      int nbfso, int elbert, int intermediate, int no_pq_perm, 
      int qdim, int add, int printflg, FILE *outfile) 
{
   int i;
   Value *valptr;              /* array of integral values */
   Label *lblptr;              /* array of integral labels */
   int idx;                    /* index for curr integral (0..ints_per_buf) */
   int lastbuf;                /* last buffer flag */
   int p, q, qmax, qmin, r, rmin, rmax, s, smin, smax, pq, rs;
   long int pqrs, offset;
   int first_p, first_q, first_pq, last_p, last_q;
   int nbstri;

   if (printflg) {
     fprintf(outfile, "\nsortbuf for pq=%d to %d\n", fpq, lpq);
   }

   if (no_pq_perm && !intermediate) {
     fprintf(outfile,"(sortbuf): illegal parameter combination.\n");
     fprintf(stderr, "(sortbuf): illegal parameter combination.\n");
   }
   
   nbstri = nbfso * (nbfso + 1) / 2;
   
   /* figure out ranges on things */
   /* I believe this section works fine, even with different ioff arrays */
   i = 0;
   while (fpq >= ioff[i] && i < BIGNUM) i++;
   if (i == BIGNUM) {
     fprintf(outfile, "(sortbuf): parameter error\n") ;
     return;
   }
   first_p = i-1 ; first_q = fpq - ioff[i-1];
   first_pq = ioff[first_p] + first_q;
   if (first_pq != fpq) {
     fprintf(outfile, "(sortbuf): fpq != first_pq.\n");
     fprintf(stderr,  "(sortbuf): fpq != first_pq.\n");
   }
   
   if (!intermediate) {
     if (elbert) offset = ioff2[first_pq] + first_pq ;
     else offset = ioff[first_pq] ;
   }
   else offset = 0;

   i=0; 
   while (lpq >= ioff[i] && i < BIGNUM) i++ ;
   if (i == BIGNUM) {
     fprintf(outfile, "(sortbuf): parameter error\n") ;
     return ;
   }
   last_p = i-1 ; last_q = lpq - ioff[i-1] ;
   
   
   lblptr = Inbuf->labels;
   valptr = Inbuf->values;

   /* read a buffer at a time until we're done */

   do {
      iwl_buf_fetch(Inbuf);
      lastbuf = Inbuf->lastbuf;
      for (idx=4*Inbuf->idx; Inbuf->idx<Inbuf->inbuf; Inbuf->idx++) {
	p = (int) lblptr[idx++];
	q = (int) lblptr[idx++];
	r = (int) lblptr[idx++];
	s = (int) lblptr[idx++];

	/* if (no_pq_perm) ioff is the appropriate offset array for the left
           indices (ioff[p] = nvirt * p for MP2); ioff2 is then the usual
           ioff offset array, used for the right indices */
	if (no_pq_perm) {
	  pq = ioff[p] + q;
	  rs = ioff2[MAX0(r,s)] + MIN0(r,s);
	}
	else {
	  pq = ioff[MAX0(p,q)] + MIN0(p,q);
	  rs = ioff[MAX0(r,s)] + MIN0(r,s);
	}
	
	if (!intermediate) {
	  if (elbert)
	    pqrs = ioff2[pq] + rs;
	  else {
	    pqrs = ioff[MAX0(pq,rs)];
	    pqrs += MIN0(pq,rs);
	  }
	}
	else {
	  pqrs = (pq - first_pq);
          pqrs *= nbstri;
          pqrs += rs;
	}
	
        if (printflg && ints[pqrs-offset] != 0.0) 
	   fprintf(outfile, "Adding %10.6lf to el %d %d %d %d = %10.6lf\n", 
                   valptr[Inbuf->idx], p, q, r, s, ints[pqrs-offset]);

        if (add) ints[pqrs-offset] += valptr[Inbuf->idx];
        else ints[pqrs-offset] += valptr[Inbuf->idx];

	if (printflg) 
	  fprintf(outfile, "<%d %d %d %d | %d %d [%ld] = %10.6lf\n",
		  p, q, r, s, pq, rs, pqrs, ints[pqrs-offset]) ;
      }
   } while (!lastbuf);
   
   /* now write them out again, in order */
   lblptr = Outbuf->labels;
   valptr = Outbuf->values;

   idx = 0;

   for (p=first_p; p<=last_p; p++) {
     qmax = (p==last_p) ? last_q : p ;
     qmin = (p==first_p) ? first_q : 0 ;
     if(no_pq_perm) {
       qmax = (p==last_p) ? last_q : (qdim-1);
       qmin = (p==first_p) ? first_q: 0;
     }
     for (q=qmin; q<=qmax; q++) {
       pq = ioff[p] + q ; /* This should be fine even with MP2 */
       
       if (!intermediate) {
	 rmin = (elbert) ? p : 0 ; 
	 rmax = (elbert) ? nbfso : p+1 ;
       }
       else {  /* This should be fine with MP2, also */
	 rmin = 0;
	 rmax = nbfso;
       }
       
       for (r=rmin; r<rmax; r++) {
	 
	 if (!intermediate) {
	   if (elbert) {
	     smax = r+1 ;
	     smin = (p==r) ? q : 0 ;
	   }
	   else {
	     smax = (p==r) ? (q+1) : (r+1) ;
	     smin = 0 ;
	   }
	 }
	 else { /* This should be fine with MP2, also */
	   smax = r + 1;
	   smin = 0;
	 }
	 
	 for (s=smin; s < smax; s++) {
	   if(no_pq_perm) rs = ioff2[r] + s;
	   else rs = ioff[r] + s;
	   
	   /* Again, this should be fine with MP2 */
	   if (elbert) pqrs = ioff2[pq] + rs ;
	   else if (intermediate) {
	     pqrs = (pq - first_pq);
	     pqrs *= nbstri;
	     pqrs += rs;
	   }
	   else pqrs = ioff[pq] + rs ;
	   
	   if (fabs(ints[pqrs-offset]) > Outbuf->cutoff) {
	     idx = 4*Outbuf->idx;
	     lblptr[idx++] = p;
	     lblptr[idx++] = q;
	     lblptr[idx++] = r;
	     lblptr[idx++] = s;
	     valptr[Outbuf->idx] = ints[pqrs-offset];
	     if (printflg) 
	       fprintf(outfile, ">%d %d %d %d | %d %d [%ld] = %10.6lf\n",
		       p, q, r, s, pq, rs, pqrs, ints[pqrs-offset]) ;
	     
	     Outbuf->idx++;
	     if (Outbuf->idx == Outbuf->ints_per_buf) {
	       Outbuf->lastbuf = 0;
	       Outbuf->inbuf = Outbuf->idx;
	       iwl_buf_put(Outbuf);
	       Outbuf->idx = 0;
	     } 
	   }
	 }
       }
     }
   }
   
}

}

