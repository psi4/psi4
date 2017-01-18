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

void IWL::write_mp2(int p, int q, int pq,
    int pqsym, double **arr, int rsym, int *firstr, int *lastr,
    int *firsts, int *lasts, int *occ, int *vir, int *ioff,
    int printflag, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
    int idx, r, s, rs, ssym;
    int R,S,rnew,snew;
    double value;
    Label *lblptr;
    Value *valptr;

    lblptr = labels_;
    valptr = values_;

    ssym = pqsym ^ rsym;
    for (r=firstr[rsym],R=0; r <= lastr[rsym]; r++,R++) {
        rnew = occ[r];
        for (s=firsts[ssym],S=0; s <=lasts[ssym]; s++,S++) {
            snew = vir[s];
            rs = ioff[rnew] + snew;
            /*------------------------------------------
            We do not need integrals with rs > pq
            rs can only increase, hence if rs > pq -
            it is time to leave
            ------------------------------------------*/
            if (rs > pq)
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
                    printer->Printf( "<%d %d %d %d [%d] [%d] = %20.10f\n",
                    p, q, rnew, snew, pq, rs, value);

            } /* end if (fabs(value) > Buf->cutoff) ... */
        } /* end loop over s */
    } /* end loop over r */
}

/*!
** iwl_buf_wrt_mp2()
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
** \ingroup IWL
*/
void iwl_buf_wrt_mp2(struct iwlbuf *Buf, int p, int q, int pq, int pqsym,
   double **arr, int rsym, int *firstr, int *lastr, int *firsts, int *lasts,
   int *occ, int *vir, int *ioff, int printflag, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
   int idx, r, s, rs, ssym;
   int R,S,rnew,snew;
   double value;
   Label *lblptr;
   Value *valptr;

   lblptr = Buf->labels;
   valptr = Buf->values;

   ssym = pqsym ^ rsym;
   for (r=firstr[rsym],R=0; r <= lastr[rsym]; r++,R++) {
     rnew = occ[r];
     for (s=firsts[ssym],S=0; s <=lasts[ssym]; s++,S++) {
       snew = vir[s];
       rs = ioff[rnew] + snew;
       /*------------------------------------------
	 We do not need integrals with rs > pq
	 rs can only increase, hence if rs > pq -
	 it is time to leave
	------------------------------------------*/
       if (rs > pq)
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
	   printer->Printf( "<%d %d %d %d [%d] [%d] = %20.10f\n",
		   p, q, rnew, snew, pq, rs, value);

       } /* end if (fabs(value) > Buf->cutoff) ... */
     } /* end loop over s */
   } /* end loop over r */

}

}
