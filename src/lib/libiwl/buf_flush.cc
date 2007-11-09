/*!
  \file buf_flush.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <libciomr/libciomr.h>
#include "iwl.h"

extern "C" {
	
/*!
** iwl_buf_flush()
**
**	\param Buf     To be flushed buffer
**	\param lastbuf Flag for the last buffer
**
** Flush an Integrals With Labels Buffer
** All flushing should be done through this routine!
** David Sherrill, March 1995
** \ingroup (IWL)
*/
void iwl_buf_flush(struct iwlbuf *Buf, int lastbuf)
{
  int idx;
  Label *lblptr;
  Value *valptr;
  
  Buf->inbuf = Buf->idx;
  lblptr = Buf->labels;
  valptr = Buf->values;
  
  idx = 4 * Buf->idx;

  while (Buf->idx < Buf->ints_per_buf) {
    lblptr[idx++] = 0;
    lblptr[idx++] = 0;
    lblptr[idx++] = 0;
    lblptr[idx++] = 0;
    valptr[Buf->idx] = 0.0;
    Buf->idx++;
  }
  
  if (lastbuf) Buf->lastbuf = 1;
  else Buf->lastbuf = 0;

  iwl_buf_put(Buf);
  Buf->idx = 0;
}

} /* extern "C" */