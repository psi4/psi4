/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "iwl.h"
#include "iwl.hpp"
#include "psi4/libparallel/ParallelPrinter.h"
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

#define BIGNUM 40000

#if !defined( EXPLICIT_IOFF )
#   define EXPLICIT_IOFF(i) ( (i) * ((i) + 1) / 2 )
#endif

#if !defined( INDEX2 )
#   define INDEX2(i, j) ( (i) >= (j) ? EXPLICIT_IOFF(i) + (j) : EXPLICIT_IOFF(j) + (i) )
#endif
#include "psi4/libqt/qt.h"
namespace psi {

void IWL::sort_buffer(IWL *Inbuf, IWL *Outbuf,
    double *ints, int fpq, int lpq, int *ioff, int *ioff2,
    int nbfso, int elbert, int intermediate, int no_pq_perm,
    int qdim, int add, int printflg, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
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
        printer->Printf( "\nsortbuf for pq=%d to %d\n", fpq, lpq);
    }

    if (no_pq_perm && !intermediate) {
        printer->Printf("(sortbuf): illegal parameter combination.\n");
        outfile->Printf( "(sortbuf): illegal parameter combination.\n");
    }

    nbstri = nbfso * (nbfso + 1) / 2;

    /* figure out ranges on things */
    /* I believe this section works fine, even with different ioff arrays */
    i = 0;
    while (fpq >= ioff[i] && i < BIGNUM) i++;
    if (i == BIGNUM) {
        printer->Printf( "(sortbuf): parameter error\n") ;
        return;
    }
    first_p = i-1 ; first_q = fpq - ioff[i-1];
    first_pq = ioff[first_p] + first_q;
    if (first_pq != fpq) {
        printer->Printf( "(sortbuf): fpq != first_pq.\n");
        outfile->Printf(  "(sortbuf): fpq != first_pq.\n");
    }

    if (!intermediate) {
        if (elbert) offset = ioff2[first_pq] + first_pq ;
        else offset = ioff[first_pq] ;
    }
    else offset = 0;

    i=0;
    while (lpq >= ioff[i] && i < BIGNUM) i++ ;
    if (i == BIGNUM) {
        printer->Printf( "(sortbuf): parameter error\n") ;
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
                printer->Printf( "Adding %10.6f to el %d %d %d %d = %10.6f\n",
                valptr[Inbuf->idx_], p, q, r, s, ints[pqrs-offset]);

            if (add) ints[pqrs-offset] += valptr[Inbuf->idx_];
            else ints[pqrs-offset] += valptr[Inbuf->idx_];

            if (printflg)
                printer->Printf( "<%d %d %d %d | %d %d [%ld] = %10.6f\n",
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
                            printer->Printf( ">%d %d %d %d | %d %d [%ld] = %10.6f\n",
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
** sortbuf_pk()
**
** Function reads a file of two-electron integrals into
** core and writes them back out again in canonical order.  Used in
** Yoshimine PK sort where we have a file containing a few value of pq
** and all corresponding rs <= pq. We sort to make sure orbitals are
** directly readable to form Coulomb and exchange matrices, thus with
** pq >= rs, p >= q and r >= s. In addition, diagonal elements are
** multiplied by the appropriate factors.
** The integrals are written to the PK file directly, in the relevant
** entry, without labels.
**
** One interesting issue here is that the intermediate array ('ints')
** must be big enough to hold the integrals in the current buffer, but
** we don't generally want it to be much larger than necessary!  Thus
** we calculate an 'offset' which is the canonical index of the first
** integral in the buffer, and we use this so that the first integral
** in the buffer is stored in ints[0].
**
**    \param Inbuf       = IWL buffer for input
**    \param out_tape    = PK file number for output
**    \param is_exch     = 0 if we are sorting Coulomb, otherwise 1
**    \param ints        = array to hold integrals in
**    \param fpq         = first pq for this tape
**    \param lpq         = last pq for this tape
**    \param so2rel      = array mapping absolute basis function index to relative
**                         basis function index within an irrep, so2rel[abs] = rel
**    \param so2sym      = array mapping absolute basis function index to irrep
**                         number, so2sym[abs] = sym
**    \param pksymoff    = array containing the offset in each irrep to convert a
**                         pq index computed with relative indices to an absolute
**                         pq index, pqrel = ioff[prel] + qrel, pqabs = pqrel + pksymoff[psym]
**    \param ioff        = offset array for the left indices
**    \param num_so      = number of basis functions per irrep
**    \param qdim        = dimensions for the q index...nvirt for MP2
**    \param printflg    = 1 for printing, 0 otherwise
**    \param out     = output file pointer
**
** Returns: none
**
** N.B. No need to iwl_flush the output buffer...not done in here!!
** \ingroup IWL
*/
void IWL::sort_buffer_pk(IWL *Inbuf, int out_tape, int is_exch, double *ints,
                         unsigned int fpq, unsigned int lpq, int *so2ind, int *so2sym,
                         int *pksymoff, int printflg, std::string out)
{
   Value *valptr;              /* array of integral values */
   Label *lblptr;              /* array of integral labels */
   int idx;                    /* index for curr integral (0..ints_per_buf) */
   int lastbuf;                /* last buffer flag */
   int pabs, qabs, rabs, sabs;
   int prel, qrel, rrel, srel, psym, qsym, rsym, ssym;
   unsigned long int pq, rs, pqrs, offset, maxind;

   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));

   if (printflg) {
     printer->Printf( "\nsortbuf_pk for pq=%d to %d\n", fpq, lpq);
   }

   // Compute the index of the lowest pq for offset

   offset = (size_t)fpq * (fpq + 1L) / 2L;

   maxind = INDEX2(lpq, lpq);

//   outfile->Printf("Offset is %lu and maxind is %lu.\n", offset, maxind);
//
   lblptr = Inbuf->labels();
   valptr = Inbuf->values();

   /* Read all integrals from a bucket file,
    * one buffer at a time until we're done
      We sort them upon reading in the twoel array */

   do {
       timer_on("Reading the buckets");
       Inbuf->fetch();
       timer_off("Reading the buckets");
      lastbuf = Inbuf->last_buffer();
      for (idx=4*Inbuf->index(); Inbuf->index()<Inbuf->buffer_count(); Inbuf->index()++) {
          pabs = (int) lblptr[idx++];
          qabs = (int) lblptr[idx++];
          rabs = (int) lblptr[idx++];
          sabs = (int) lblptr[idx++];

          // Get indices within symmetry

          prel = so2ind[pabs];
          qrel = so2ind[qabs];
          rrel = so2ind[rabs];
          srel = so2ind[sabs];

          psym = so2sym[pabs];
          qsym = so2sym[qabs];
          rsym = so2sym[rabs];
          ssym = so2sym[sabs];

//DEBUG          outfile->Printf("INT <%d %d|%d %d>\n", pabs, qabs, rabs, sabs);

          if (!is_exch) {

              // We only stored integrals of the relevant symmetry
              pq = INDEX2(prel, qrel);
              pq += pksymoff[psym];
              rs = INDEX2(rrel, srel);
              rs += pksymoff[rsym];
              pqrs = INDEX2(pq, rs);
//DEBUG         if (pqrs > maxind || (pqrs < offset)) {
//DEBUG             outfile->Printf("pqrs is out of bounds for J\n");
//DEBUG         }
              ints[pqrs - offset] += valptr[Inbuf->index()];

          } else {
              // K (2nd sort, ILJK)
              if ((psym == qsym) && (rsym == ssym)) {
                  if ( (prel != qrel) && (rrel != srel)) {
                      if((psym == ssym) && (qsym == rsym)) {
                          pq = INDEX2(prel, srel);
                          pq += pksymoff[psym];
                          rs = INDEX2(qrel, rrel);
                          rs += pksymoff[qsym];
                          pqrs = INDEX2(pq, rs);
                          if ((pqrs <= maxind) && (pqrs >= offset)) {
//                              if(rs > pq) {
//                                outfile->Printf("rs > pq 2nd sort!!\n");
//                              }
                              if(prel == srel || qrel == rrel) {
                                  ints[pqrs - offset] += valptr[Inbuf->index()];
                              } else {
                                  ints[pqrs - offset] += 0.5 * valptr[Inbuf->index()];
                              }
                          }
                      }
                  }
              }
              // K (1st sort, IKJL)
              if ((psym == rsym) && (qsym == ssym)) {
                  pq = INDEX2(prel, rrel);
                  pq += pksymoff[psym];
                  rs = INDEX2(qrel, srel);
                  rs += pksymoff[qsym];
                  pqrs = INDEX2(pq, rs);
                  if ((pqrs <= maxind) && (pqrs >= offset) ) {
//                      if(rs > pq) {
//                        outfile->Printf("rs > pq 1st sort!!\n");
//                      }
                      if((prel == rrel) || (qrel == srel)) {
                          ints[pqrs - offset] += valptr[Inbuf->index()];
                      } else {
                          ints[pqrs - offset] += 0.5 * valptr[Inbuf->index()];
                      }
                  }
              }

          }


          if (printflg)
            printer->Printf( "<%d %d %d %d | %d %d [%ld] = %10.6f\n",
                pabs, qabs, rabs, sabs, pq, rs, pqrs, ints[pqrs-offset]) ;
      }
   } while (!lastbuf);

   for(pq = fpq; pq <= lpq; ++pq) {
       pqrs = INDEX2(pq, pq);
       ints[pqrs - offset] *= 0.5;
   }

//DEBUG   outfile->Printf("Final pqrs value %i\n", pqrs);
   /* That's all we do here. The integral array is now ready
    * to be written in the appropriate file. */

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
**    \param elbert      = integrals obey rs >= pq.  Use ioff2 to get offset.
**    \param no_pq_perm  = don't use permutational symmetry to swap p and q
**                  (appropriate for MP2 where one is occ and one is virt)
**    \param qdim        = dimensions for the q index...nvirt for MP2
**    \param add         = add contributions to the same integral during sort
**    \param printflg    = 1 for printing, 0 otherwise
**    \param out     = output file pointer
**
** Returns: none
**
** Revised 6/27/96 by CDS for new IWL format
** N.B. No need to iwl_flush the output buffer...not done in here!!
** \ingroup IWL
*/
void sortbuf(struct iwlbuf *Inbuf, struct iwlbuf *Outbuf,
      double *ints, int fpq, int lpq, int *ioff, int *ioff2,
      int nbfso, int elbert, int intermediate, int no_pq_perm,
      int qdim, int add, int printflg, std::string out)
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
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   if (printflg) {
     printer->Printf( "\nsortbuf for pq=%d to %d\n", fpq, lpq);
   }

   if (no_pq_perm && !intermediate) {
     printer->Printf("(sortbuf): illegal parameter combination.\n");
     outfile->Printf( "(sortbuf): illegal parameter combination.\n");
   }

   nbstri = nbfso * (nbfso + 1) / 2;

   /* figure out ranges on things */
   /* I believe this section works fine, even with different ioff arrays */
   i = 0;
   while (fpq >= ioff[i] && i < BIGNUM) i++;
   if (i == BIGNUM) {
     printer->Printf( "(sortbuf): parameter error\n") ;
     return;
   }
   first_p = i-1 ; first_q = fpq - ioff[i-1];
   first_pq = ioff[first_p] + first_q;
   if (first_pq != fpq) {
     printer->Printf( "(sortbuf): fpq != first_pq.\n");
     outfile->Printf(  "(sortbuf): fpq != first_pq.\n");
   }

   if (!intermediate) {
     if (elbert) offset = ioff2[first_pq] + first_pq ;
     else offset = ioff[first_pq] ;
   }
   else offset = 0;

   i=0;
   while (lpq >= ioff[i] && i < BIGNUM) i++ ;
   if (i == BIGNUM) {
     printer->Printf( "(sortbuf): parameter error\n") ;
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
	   printer->Printf( "Adding %10.6f to el %d %d %d %d = %10.6f\n",
                   valptr[Inbuf->idx], p, q, r, s, ints[pqrs-offset]);

        if (add) ints[pqrs-offset] += valptr[Inbuf->idx];
        else ints[pqrs-offset] += valptr[Inbuf->idx];

	if (printflg)
	  printer->Printf( "<%d %d %d %d | %d %d [%ld] = %10.6f\n",
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
	       printer->Printf( ">%d %d %d %d | %d %d [%ld] = %10.6f\n",
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

/*!
** sortbuf_pk()
**
** Function reads a file of two-electron integrals into
** core and writes them back out again in canonical order.  Used in
** Yoshimine PK sort where we have a file containing a few value of pq
** and all corresponding rs <= pq. We sort to make sure orbitals are
** directly readable to form Coulomb and exchange matrices, thus with
** pq >= rs, p >= q and r >= s. In addition, diagonal elements are
** multiplied by the appropriate factors.
** The integrals are written to the PK file directly, in the relevant
** entry, without labels.
**
** One interesting issue here is that the intermediate array ('ints')
** must be big enough to hold the integrals in the current buffer, but
** we don't generally want it to be much larger than necessary!  Thus
** we calculate an 'offset' which is the canonical index of the first
** integral in the buffer, and we use this so that the first integral
** in the buffer is stored in ints[0].
**
**    \param Inbuf       = IWL buffer for input
**    \param out_tape    = PK file number for output
**    \param is_exch     = 0 if we are sorting Coulomb, otherwise 1
**    \param ints        = array to hold integrals in
**    \param fpq         = first pq for this tape
**    \param lpq         = last pq for this tape
**    \param so2rel      = array mapping absolute basis function index to relative
**                         basis function index within an irrep, so2rel[abs] = rel
**    \param so2sym      = array mapping absolute basis function index to irrep
**                         number, so2sym[abs] = sym
**    \param pksymoff    = array containing the offset in each irrep to convert a
**                         pq index computed with relative indices to an absolute
**                         pq index, pqrel = ioff[prel] + qrel, pqabs = pqrel + pksymoff[psym]
**    \param ioff        = offset array for the left indices
**    \param num_so      = number of basis functions per irrep
**    \param qdim        = dimensions for the q index...nvirt for MP2
**    \param printflg    = 1 for printing, 0 otherwise
**    \param out     = output file pointer
**
** Returns: none
**
** N.B. No need to iwl_flush the output buffer...not done in here!!
** \ingroup IWL
*/
void sortbuf_pk(struct iwlbuf *Inbuf, int out_tape, int is_exch,
      double *ints, unsigned int fpq, unsigned int lpq, int* so2ind, int* so2sym,
      int* pksymoff, int printflg, std::string out)
{
   Value *valptr;              /* array of integral values */
   Label *lblptr;              /* array of integral labels */
   int idx;                    /* index for curr integral (0..ints_per_buf) */
   int lastbuf;                /* last buffer flag */
   int pabs, qabs, rabs, sabs;
   int prel, qrel, rrel, srel, psym, qsym, rsym, ssym;
   unsigned long int pq, rs, pqrs, offset, maxind;

   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));

   if (printflg) {
     printer->Printf( "\nsortbuf_pk for pq=%d to %d\n", fpq, lpq);
   }

   // Compute the index of the lowest pq for offset

   offset = (size_t)fpq * (fpq + 1L) / 2L;

   maxind = INDEX2(lpq, lpq);

//   outfile->Printf("Offset is %lu and maxind is %lu.\n", offset, maxind);
//
   lblptr = Inbuf->labels;
   valptr = Inbuf->values;

   /* Read all integrals from a bucket file,
    * one buffer at a time until we're done
      We sort them upon reading in the twoel array */

   do {
      timer_on("Reading the buckets");
      iwl_buf_fetch(Inbuf);
      timer_off("Reading the buckets");
      lastbuf = Inbuf->lastbuf;
      for (idx=4*Inbuf->idx; Inbuf->idx<Inbuf->inbuf; Inbuf->idx++) {
          pabs = (int) lblptr[idx++];
          qabs = (int) lblptr[idx++];
          rabs = (int) lblptr[idx++];
          sabs = (int) lblptr[idx++];

          // Get indices within symmetry

          prel = so2ind[pabs];
          qrel = so2ind[qabs];
          rrel = so2ind[rabs];
          srel = so2ind[sabs];

          psym = so2sym[pabs];
          qsym = so2sym[qabs];
          rsym = so2sym[rabs];
          ssym = so2sym[sabs];

//DEBUG          outfile->Printf("INT <%d %d|%d %d>\n", pabs, qabs, rabs, sabs);

          if (!is_exch) {

              // We only stored integrals of the relevant symmetry
              pq = INDEX2(prel, qrel);
              pq += pksymoff[psym];
              rs = INDEX2(rrel, srel);
              rs += pksymoff[rsym];
              pqrs = INDEX2(pq, rs);
//DEBUG         if (pqrs > maxind || (pqrs < offset)) {
//DEBUG             outfile->Printf("pqrs is out of bounds for J\n");
//DEBUG         }
              ints[pqrs - offset] += valptr[Inbuf->idx];

          } else {
              // K (2nd sort, ILJK)
              if ((psym == qsym) && (rsym == ssym)) {
                  if ( (prel != qrel) && (rrel != srel)) {
                      if((psym == ssym) && (qsym == rsym)) {
                          pq = INDEX2(prel, srel);
                          pq += pksymoff[psym];
                          rs = INDEX2(qrel, rrel);
                          rs += pksymoff[qsym];
                          pqrs = INDEX2(pq, rs);
                          if ((pqrs <= maxind) && (pqrs >= offset)) {
//                              if(rs > pq) {
//                                outfile->Printf("rs > pq 2nd sort!!\n");
//                              }
                              if(prel == srel || qrel == rrel) {
                                  ints[pqrs - offset] += valptr[Inbuf->idx];
                              } else {
                                  ints[pqrs - offset] += 0.5 * valptr[Inbuf->idx];
                              }
                          }
                      }
                  }
              }
              // K (1st sort, IKJL)
              if ((psym == rsym) && (qsym == ssym)) {
                  pq = INDEX2(prel, rrel);
                  pq += pksymoff[psym];
                  rs = INDEX2(qrel, srel);
                  rs += pksymoff[qsym];
                  pqrs = INDEX2(pq, rs);
                  if ((pqrs <= maxind) && (pqrs >= offset) ) {
//                      if(rs > pq) {
//                        outfile->Printf("rs > pq 1st sort!!\n");
//                      }
                      if((prel == rrel) || (qrel == srel)) {
                          ints[pqrs - offset] += valptr[Inbuf->idx];
                      } else {
                          ints[pqrs - offset] += 0.5 * valptr[Inbuf->idx];
                      }
                  }
              }

          }


          if (printflg)
            printer->Printf( "<%d %d %d %d | %d %d [%ld] = %10.6f\n",
                pabs, qabs, rabs, sabs, pq, rs, pqrs, ints[pqrs-offset]) ;
      }
   } while (!lastbuf);

   for(pq = fpq; pq <= lpq; ++pq) {
       pqrs = INDEX2(pq, pq);
       ints[pqrs - offset] *= 0.5;
   }

//DEBUG   outfile->Printf("Final pqrs value %i\n", pqrs);
   /* That's all we do here. The integral array is now ready
    * to be written in the appropriate file. */

}

}
