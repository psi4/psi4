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
namespace psi {

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

int IWL::read_all_active(double *ints,
    int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff,
    int fstact, int lstact, int printflg,std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
    int lastbuf;
    Label *lblptr;
    Value *valptr;
    int idx, p, q, r, s, pq, rs, pqrs;

    lblptr = labels_;
    valptr = values_;

    lastbuf = lastbuf_;

    for (idx=4*idx_; idx_ < inbuf_; idx_++) {
        p = (int) lblptr[idx++];
        q = (int) lblptr[idx++];
        r = (int) lblptr[idx++];
        s = (int) lblptr[idx++];

        if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
        if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
        p -= fstact;
        q -= fstact;
        r -= fstact;
        s -= fstact;

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
            printer->Printf( "<%2d %2d %2d %2d [%2d][%2d] [[%3d]] = %20.10f\n",
            p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;

    } /*! end loop through current buffer */

     /*! read new PSI buffers */
    while (!lastbuf) {
        fetch();
        lastbuf = lastbuf_;

        for (idx=4*idx_; idx_ < inbuf_; idx_++) {
            p = (int) lblptr[idx++];
            q = (int) lblptr[idx++];
            r = (int) lblptr[idx++];
            s = (int) lblptr[idx++];

            if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
            if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
            p -= fstact;
            q -= fstact;
            r -= fstact;
            s -= fstact;

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
                printer->Printf( "<%d %d %d %d [%d][%d] [[%d]] = %20.10f\n",
                p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;

        } /*! end loop through current buffer */

    } /*! end loop over reading buffers */

    return(0); /*! we must have reached the last buffer at this point */
}

/*!
** iwl_buf_rd_all_act()
**
** Read from an Integrals With Labels formatted buffer.
** The buffer must have been initialized with iwl_buf_init().
** Same as function iwl_buf_rd_all() except that we only keep the
** integrals with all-active labels.
**
**    \param Buf           =  IWL Buffer to read from (already initialized)
**    \param ints          =  memory buffer to put integrals into
**    \param ioff_lt       =  ioff array for the left pair of indices (p and q)
**    \param ioff_rt       =  ioff array for the right pair of indices (r and s)
**    \param no_pq_perm    =  if 1, do not use p/q or r/s permutational symmetry
**    \param ioff          =  the ioff array to figure the total index pqrs from
**                     the pair indices pq and rs
**    \param fstact        =  index of first active orbital
**    \param lstact        =  index of last active orbital
**    \param printflg      =  if 1, print integrals as they are read
**    \param out       =  pointer to output file for printing
**
** Returns: 0 if end of file, otherwise 1
** \ingroup IWL
*/
int iwl_buf_rd_all_act(struct iwlbuf *Buf, double *ints,
                       int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff,
                       int fstact, int lstact, int printflg,std::string out)
{
  int lastbuf;
  Label *lblptr;
  Value *valptr;
  int idx, p, q, r, s, pq, rs, pqrs;

  lblptr = Buf->labels;
  valptr = Buf->values;

  lastbuf = Buf->lastbuf;
  std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
        std::shared_ptr<OutFile>(new OutFile(out)));
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
    p = (int) lblptr[idx++];
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
    if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
    p -= fstact;
    q -= fstact;
    r -= fstact;
    s -= fstact;

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
      printer->Printf( "<%2d %2d %2d %2d [%2d][%2d] [[%3d]] = %20.10f\n",
	      p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;

  } /*! end loop through current buffer */

   /*! read new PSI buffers */
  while (!lastbuf) {
    iwl_buf_fetch(Buf);
    lastbuf = Buf->lastbuf;

    for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
      if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
      p -= fstact;
      q -= fstact;
      r -= fstact;
      s -= fstact;

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
	printer->Printf( "<%d %d %d %d [%d][%d] [[%d]] = %20.10f\n",
		p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;

    } /*! end loop through current buffer */

  } /*! end loop over reading buffers */

  return(0); /*! we must have reached the last buffer at this point */
}

}
