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
  
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

int IWL::read_all(double *ints, int *ioff_lt, int *ioff_rt, int no_pq_perm,
    int *ioff, int printflg, FILE *out)
{
    int lastbuf;
    Label *lblptr;
    Value *valptr;
    int idx, p, q, r, s, pq, rs, pqrs;

    lblptr = labels_;
    valptr = values_;

    lastbuf = lastbuf_;

    for (idx=4*idx_; idx_ < inbuf_; idx_++) {
        p = std::abs(lblptr[idx++]);
        q = (int) lblptr[idx++];
        r = (int) lblptr[idx++];
        s = (int) lblptr[idx++];

        if(no_pq_perm) { /*! I _think_ this will work */
            pq = ioff_lt[p] + q;
            rs = ioff_rt[r] + s;
        }
        else {
            pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
            rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
        }

        pqrs = INDEX(pq,rs);

        ints[pqrs] = (double) valptr[idx_];

        if (printflg) 
            fprintf(out, "<%2d %2d %2d %2d [%2d][%2d] [[%3d]] = %20.10f\n",
            p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;

    } /*! end loop through current buffer */

    /*! read new PSI buffers */
    while (!lastbuf) {
        fetch();
        lastbuf = lastbuf_;

        for (idx=4*idx_; idx_ < inbuf_; idx_++) {
            p = std::abs(lblptr[idx++]);
            q = (int) lblptr[idx++];
            r = (int) lblptr[idx++];
            s = (int) lblptr[idx++];

            if(no_pq_perm) { /*! I _think_ this will work */
                pq = ioff_lt[p] + q;
                rs = ioff_rt[r] + s;
            }
            else {
                pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
                rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
            }

            pqrs = INDEX(pq,rs);

            ints[pqrs] = (double) valptr[idx_];

            if (printflg) 
                fprintf(out, "<%d %d %d %d [%d][%d] [[%d]] = %20.10f\n",
                p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;

        } /*! end loop through current buffer */

    } /*! end loop over reading buffers */

    return(0); /*! we must have reached the last buffer at this point */
}

int IWL::read_all2(double **ints, int *ioff_lt, int *ioff_rt, int no_pq_perm, 
    int *, int printflg, FILE *out)
{
    int lastbuf;
    Label *lblptr;
    Value *valptr;
    int idx, p, q, r, s, pq, rs;

    lblptr = labels_;
    valptr = values_;

    lastbuf = lastbuf_;

    for (idx=4*idx_; idx_ < inbuf_; idx_++) {
        p = std::abs(lblptr[idx++]);
        q = (int) lblptr[idx++];
        r = (int) lblptr[idx++];
        s = (int) lblptr[idx++];

        if(no_pq_perm) { /*! I _think_ this will work */
            pq = ioff_lt[p] + q;
            rs = ioff_rt[r] + s;
        }
        else {
            pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
            rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
        }

        ints[pq][rs] = (double) valptr[idx_];

        if (printflg) 
            fprintf(out, "<%2d %2d %2d %2d [%2d][%2d] = %20.10f\n",
            p, q, r, s, pq, rs, ints[pq][rs]) ;

    } /*! end loop through current buffer */

     /*! read new PSI buffers */
    while (!lastbuf) {
        fetch();
        lastbuf = lastbuf_;

        for (idx=4*idx_; idx_ < inbuf_; idx_++) {
            p = std::abs(lblptr[idx++]);
            q = (int) lblptr[idx++];
            r = (int) lblptr[idx++];
            s = (int) lblptr[idx++];

            if(no_pq_perm) { /*! I _think_ this will work */
                pq = ioff_lt[p] + q;
                rs = ioff_rt[r] + s;
            }
            else {
                pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
                rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
            }

            ints[pq][rs] = (double) valptr[idx_];

            if (printflg) 
                fprintf(out, "<%d %d %d %d [%d][%d] = %20.10f\n",
                p, q, r, s, pq, rs, ints[pq][rs]) ;

        } /*! end loop through current buffer */

    } /*! end loop over reading buffers */

    return(0); /*! we must have reached the last buffer at this point */
}

/*!
** iwl_buf_rd_all()
**
** Read from an Integrals With Labels formatted buffer.
** The buffer must have been initialized with iwl_buf_init().
**
** Arguments:
**    \param Buf           =  IWL Buffer to read from (already initialized)
**    \param ints          =  memory buffer to put integrals into
**    \param ioff_lt       =  ioff array for the left pair of indices (p and q)
**    \param ioff_rt       =  ioff array for the right pair of indices (r and s)
**    \param no_pq_perm    =  if 1, do not use p/q or r/s permutational symmetry
**    \param ioff          =  the ioff array to figure the total index pqrs from
**                     the pair indices pq and rs
**    \param printflg      =  if 1, print integrals as they are read
**    \param out       =  pointer to output file for printing
**
** Returns: 0 if end of file, otherwise 1
** \ingroup IWL
*/
int iwl_buf_rd_all(struct iwlbuf *Buf, double *ints,
		   int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff,
                   int printflg, FILE *out)
{
  int lastbuf;
  Label *lblptr;
  Value *valptr;
  int idx, p, q, r, s, pq, rs, pqrs;
  
  lblptr = Buf->labels;
  valptr = Buf->values;
  
  lastbuf = Buf->lastbuf;
  
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
    p = std::abs(lblptr[idx++]);
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    if(no_pq_perm) { /*! I _think_ this will work */
      pq = ioff_lt[p] + q;
      rs = ioff_rt[r] + s;
    }
    else {
      pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
      rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
    }
    
    pqrs = INDEX(pq,rs);

    ints[pqrs] = (double) valptr[Buf->idx];
    
    if (printflg) 
      fprintf(out, "<%2d %2d %2d %2d [%2d][%2d] [[%3d]] = %20.10f\n",
	      p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;
    
  } /*! end loop through current buffer */
  
  /*! read new PSI buffers */
  while (!lastbuf) {
    iwl_buf_fetch(Buf);
    lastbuf = Buf->lastbuf;
    
    for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = std::abs(lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      if(no_pq_perm) { /*! I _think_ this will work */
	pq = ioff_lt[p] + q;
	rs = ioff_rt[r] + s;
      }
      else {
	pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
	rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
      }
      
      pqrs = INDEX(pq,rs);

      ints[pqrs] = (double) valptr[Buf->idx];
      
      if (printflg) 
	fprintf(out, "<%d %d %d %d [%d][%d] [[%d]] = %20.10f\n",
		p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;
      
    } /*! end loop through current buffer */
    
  } /*! end loop over reading buffers */
  
  return(0); /*! we must have reached the last buffer at this point */
}

/*!
** IWL_BUF_RD_ALL2(): This routine works exactly like
** iwl_buf_rd_all(), except that the integral list is not assumed to
** have bra-ket permutational symmetry.  The list is still required to
** have permutational symmetry WITHIN bra and ket, however, unless
** no_pq_perm is set.  This function requires that the input array be
** (double **) rather than (double *).  This routine is necessary, for
** example, for reading the alpha-beta two-electron integrals from the
** UHF transqt code.
**
** TDC, 6/01
** \ingroup IWL
*/

int iwl_buf_rd_all2(struct iwlbuf *Buf, double **ints,
		   int *ioff_lt, int *ioff_rt, int no_pq_perm, int *,
                   int printflg, FILE *out)
{
  int lastbuf;
  Label *lblptr;
  Value *valptr;
  int idx, p, q, r, s, pq, rs;
  
  lblptr = Buf->labels;
  valptr = Buf->values;
  
  lastbuf = Buf->lastbuf;
  
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
    p = std::abs(lblptr[idx++]);
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    if(no_pq_perm) { /*! I _think_ this will work */
      pq = ioff_lt[p] + q;
      rs = ioff_rt[r] + s;
    }
    else {
      pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
      rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
    }
    
    ints[pq][rs] = (double) valptr[Buf->idx];
    
    if (printflg) 
      fprintf(out, "<%2d %2d %2d %2d [%2d][%2d] = %20.10f\n",
	      p, q, r, s, pq, rs, ints[pq][rs]) ;
    
  } /*! end loop through current buffer */
  
   /*! read new PSI buffers */
  while (!lastbuf) {
    iwl_buf_fetch(Buf);
    lastbuf = Buf->lastbuf;
    
    for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = std::abs(lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      if(no_pq_perm) { /*! I _think_ this will work */
	pq = ioff_lt[p] + q;
	rs = ioff_rt[r] + s;
      }
      else {
	pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
	rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
      }
      
      ints[pq][rs] = (double) valptr[Buf->idx];
      
      if (printflg) 
	fprintf(out, "<%d %d %d %d [%d][%d] = %20.10f\n",
		p, q, r, s, pq, rs, ints[pq][rs]) ;
      
    } /*! end loop through current buffer */
    
  } /*! end loop over reading buffers */
  
  return(0); /*! we must have reached the last buffer at this point */
}

}

