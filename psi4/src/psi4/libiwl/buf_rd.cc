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

namespace psi {

int IWL::read(int target_pq, double *ints, int *ioff_lt, int *ioff_rt,
    int mp2, int printflg,std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
    int lastbuf;
    Value *valptr;
    Label *lblptr;
    int idx, p, q, r, s, pq, rs;

    lblptr = labels_;
    valptr = values_;

    lastbuf = lastbuf_;

    for (idx=4*idx_; idx_<inbuf_; idx_++) {
        p = (int) lblptr[idx++];
        q = (int) lblptr[idx++];
        r = (int) lblptr[idx++];
        s = (int) lblptr[idx++];

        if(mp2) { /*! I _think_ this will work */
            pq = ioff_lt[p] + q;
            rs = ioff_rt[r] + s;
        }
        else {
            pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
            rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
        }

      /*!      if (pq < target_pq) continue; */
        if (pq != target_pq) return(1);

        ints[rs] = (double) valptr[idx_];

        if (printflg)
            printer->Printf( "<%d %d %d %d [%d][%d] = %20.10f\n",
            p, q, r, s, pq, rs, ints[rs]) ;

    } /*! end loop through current buffer */

    /*! read new buffers */
    while (!lastbuf) {
        fetch();
        lastbuf = lastbuf_;

        for (idx=4*idx_; idx_ < inbuf_; idx_++) {
            p = (int) lblptr[idx++];
            q = (int) lblptr[idx++];
            r = (int) lblptr[idx++];
            s = (int) lblptr[idx++];

            if(mp2) { /*! I _think_ this will work */
                pq = ioff_lt[p] + q;
                rs = ioff_rt[r] + s;
            }
            else {
                pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
                rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
            }

            if (pq < target_pq) continue;
            if (pq > target_pq) return(1);

            ints[rs] = (double) valptr[idx_];

            if (printflg)
                printer->Printf( "<%d %d %d %d [%d][%d] = %20.10f\n",
                p, q, r, s, pq, rs, ints[rs]) ;

        } /*! end loop through current buffer */

    } /*! end loop over reading buffers */

    return(0); /*! we must have reached the last buffer at this point */
}

/*!
** iwl_buf_rd(struct iwlbuf *Buf, int target_pq, double *ints,
**	       int *ioff_lt, int *ioff_rt, int mp2, int printflg,
**	       FILE *out)
**
**
** Read from an Integrals With Labels formatted PSI buffer.
** The buffer must have been initialized with iwl_buf_init().
** David Sherrill, March 1995
**
** Returns: 0 if end of file, otherwise 1
**
** Altered such that if the current pq value does not equal the target_pq
** then routine returns.  This may be dangerous in that if you don't know
** the order of pq's in the iwl_buf, you may skip integrals!
** -Daniel, November 9, 1995
**
** Revised 6/26/96 by CDS for new format
** \ingroup IWL
*/
int iwl_buf_rd(struct iwlbuf *Buf, int target_pq, double *ints,
	       int *ioff_lt, int *ioff_rt, int mp2, int printflg,
	       std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
  int lastbuf;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s, pq, rs;

  lblptr = Buf->labels;
  valptr = Buf->values;

  lastbuf = Buf->lastbuf;

  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
    p = (int) lblptr[idx++];
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    if(mp2) { /*! I _think_ this will work */
      pq = ioff_lt[p] + q;
      rs = ioff_rt[r] + s;
    }
    else {
      pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
      rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
    }

    /*!      if (pq < target_pq) continue; */
    if (pq != target_pq) return(1);

    ints[rs] = (double) valptr[Buf->idx];

    if (printflg)
      printer->Printf( "<%d %d %d %d [%d][%d] = %20.10f\n",
	      p, q, r, s, pq, rs, ints[rs]) ;

  } /*! end loop through current buffer */

  /*! read new buffers */
  while (!lastbuf) {
    iwl_buf_fetch(Buf);
    lastbuf = Buf->lastbuf;

    for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      if(mp2) { /*! I _think_ this will work */
	pq = ioff_lt[p] + q;
	rs = ioff_rt[r] + s;
      }
      else {
	pq = ioff_lt[MAX0(p,q)] + MIN0(p,q);
	rs = ioff_rt[MAX0(r,s)] + MIN0(r,s);
      }

      if (pq < target_pq) continue;
      if (pq > target_pq) return(1);

      ints[rs] = (double) valptr[Buf->idx];

      if (printflg)
	printer->Printf( "<%d %d %d %d [%d][%d] = %20.10f\n",
		p, q, r, s, pq, rs, ints[rs]) ;

    } /*! end loop through current buffer */

  } /*! end loop over reading buffers */

  return(0); /*! we must have reached the last buffer at this point */
}

}
