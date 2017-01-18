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

int IWL::read_array2(double *ints, int *plist, int *qlist, int *rlist, int *slist,
    int *size, int *ioff, int printflg, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
    int lastbuf;
    Label *lblptr;
    Value *valptr;
    int idx, p, q, pq, r, s;
    double value;

    lblptr = labels_;
    valptr = values_;

    lastbuf = lastbuf_;

    *size = 0;

    for (idx=4*idx_; idx_ < inbuf_; idx_++) {
        p = (int) lblptr[idx++];
        q = (int) lblptr[idx++];
        r = (int) lblptr[idx++];
        s = (int) lblptr[idx++];

        pq = INDEX(p,q);

        value = (double) valptr[idx_];
        *ints++ = value;
        *plist++ = p;
        *qlist++ = q;
        *rlist++ = r;
        *slist++ = s;
        *size= *size + 1;

        if (printflg)
            printer->Printf( "<%d %d %d %d [%d] = %20.10f\n", p, q, r, s,
            pq, value);
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

            pq = INDEX(p,q);

            value = (double) valptr[idx_];
            *ints++ = value;
            *plist++ = p;
            *qlist++ = q;
            *rlist++ = r;
            *slist++ = s;
            *size = *size + 1;

            if (printflg)
                printer->Printf( "<%d %d %d %d [%d] = %20.10f\n", p, q, r, s,
                pq, value);
        } /*! end loop through current buffer */
    } /*! end loop over reading buffers */

    return(0); /*! we must have reached the last buffer at this point */
}

/*!
** iwl_buf_rd_arr2()
**
** Read from an Integrals With Labels formatted PSI buffer.
** The buffer must have been initialized with iwl_buf_init().  The
** integrals and their labels are returned in the arrays ints, plist,
** qlist, rlist, and slist, and the size of these arrays is returned in
** 'size.'
**
** Returns: 0 if end of file, otherwise 1
**
** Revised 6/27/96 by CDS for new format
**
*/
int iwl_buf_rd_arr2(struct iwlbuf *Buf, double *ints, int *plist,
      int *qlist, int *rlist, int *slist, int *size, int *ioff,
      int printflg, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
   int lastbuf;
  Label *lblptr;
  Value *valptr;
  int idx, p, q, pq, r, s;
  double value;

  lblptr = Buf->labels;
  valptr = Buf->values;

  lastbuf = Buf->lastbuf;

  *size = 0;

  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
    p = (int) lblptr[idx++];
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    pq = INDEX(p,q);

    value = (double) valptr[Buf->idx];
    *ints++ = value;
    *plist++ = p;
    *qlist++ = q;
    *rlist++ = r;
    *slist++ = s;
    *size= *size + 1;

    if (printflg)
      printer->Printf( "<%d %d %d %d [%d] = %20.10f\n", p, q, r, s,
	      pq, value);

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

      pq = INDEX(p,q);

      value = (double) valptr[Buf->idx];
      *ints++ = value;
      *plist++ = p;
      *qlist++ = q;
      *rlist++ = r;
      *slist++ = s;
      *size = *size + 1;

      if (printflg)
	printer->Printf( "<%d %d %d %d [%d] = %20.10f\n", p, q, r, s,
		pq, value);

    } /*! end loop through current buffer */

  } /*! end loop over reading buffers */

  return(0); /*! we must have reached the last buffer at this point */
}

}
