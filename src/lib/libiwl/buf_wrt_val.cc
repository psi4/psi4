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
  
void IWL::write_value(int p, int q, int r, int s, double value, int printflag,
    FILE *outfile, int dirac)
{
    int idx;
    Label *lblptr;
    Value *valptr;

    lblptr = labels_;
    valptr = values_;

    if (fabs(value) > cutoff_) {
        idx = 4 * idx_;
        if(dirac) {
            lblptr[idx++] = (Label) p;
            lblptr[idx++] = (Label) r;
            lblptr[idx++] = (Label) q;
            lblptr[idx++] = (Label) s;
        }
        else {
            lblptr[idx++] = (Label) p;
            lblptr[idx++] = (Label) q;
            lblptr[idx++] = (Label) r;
            lblptr[idx++] = (Label) s;
        }
        valptr[idx_] = (Value) value;

        idx_++;

        if (idx_ == ints_per_buf_) {
            lastbuf_ = 0;
            inbuf_ = idx_;
            put();
            idx_ = 0;
        }

        if (printflag) {
            if(dirac) {
                fprintf(outfile, ">%d %d %d %d = %20.10f\n",
                    p, r, q, s, value);
            }
            else {
                fprintf(outfile, ">%d %d %d %d = %20.10f\n",
                    p, q, r, s, value);
            }
        }
    }
}

/*!
** iwl_buf_wrt_val()
**
** Write to an Integrals With Labels formatted buffer.
** The buffer must have been initialized with iwl_buf_init().  Don't
** forget to call iwl_buf_flush() when finished with all writes to the
** buffer to ensure that all contents are written to disk.
**
** This function writes only a particular value and its indices to the
** given iwl buffer.  This is useful when index rearragements are
** necessary (e.g. conversion from Mulliken to Dirac notation).  This
** is not as nice as being able to write entire arrays of values to the
** buffer, but may be necessary at times.
** Daniel Crawford, Novemeber 1995
** \ingroup IWL
*/
void iwl_buf_wrt_val(struct iwlbuf *Buf, int p, int q, int r, int s,
                     double value, int printflag, FILE *outfile, int dirac)
{
  int idx;
  Label *lblptr;
  Value *valptr;

  lblptr = Buf->labels;
  valptr = Buf->values;

  if (fabs(value) > Buf->cutoff) {
    idx = 4 * Buf->idx;
    if(dirac) {
      lblptr[idx++] = (Label) p;
      lblptr[idx++] = (Label) r;
      lblptr[idx++] = (Label) q;
      lblptr[idx++] = (Label) s;
    }
    else {
      lblptr[idx++] = (Label) p;
      lblptr[idx++] = (Label) q;
      lblptr[idx++] = (Label) r;
      lblptr[idx++] = (Label) s;
    }
    valptr[Buf->idx] = (Value) value;

    Buf->idx++;
     
    if (Buf->idx == Buf->ints_per_buf) {
      Buf->lastbuf = 0;
      Buf->inbuf = Buf->idx;
      iwl_buf_put(Buf);
      Buf->idx = 0;
    }
     
    if (printflag) {
      if(dirac) {
	fprintf(outfile, ">%d %d %d %d = %20.10lf\n",
		p, r, q, s, value);
      }
      else {
	fprintf(outfile, ">%d %d %d %d = %20.10lf\n",
		p, q, r, s, value);
      }
    }
  } 
}

}

