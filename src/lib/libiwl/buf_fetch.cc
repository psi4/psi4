/*!
  \file buf_fetch.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <libpsio/psio.h>
#include "iwl.h"

extern "C" {
	
/*!
** iwl_buf_fetch()
**
** Fetch an IWL buffer from disk
** David Sherrill, 26 June 1996
** \ingroup (IWL)
*/
void iwl_buf_fetch(struct iwlbuf *Buf)
{
  psio_read(Buf->itap, IWL_KEY_BUF, (char *) &(Buf->lastbuf), sizeof(int),
	    Buf->bufpos, &Buf->bufpos);
  psio_read(Buf->itap, IWL_KEY_BUF, (char *) &(Buf->inbuf), sizeof(int),
	    Buf->bufpos, &Buf->bufpos);
  psio_read(Buf->itap, IWL_KEY_BUF, (char *) Buf->labels, Buf->ints_per_buf * 
	    4 * sizeof(Label), Buf->bufpos, &Buf->bufpos);
  psio_read(Buf->itap, IWL_KEY_BUF, (char *) Buf->values, Buf->ints_per_buf *
	    sizeof(Value), Buf->bufpos, &Buf->bufpos);
  Buf->idx = 0;
}

} /* extern "C" */