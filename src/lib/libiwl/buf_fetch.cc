/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <libpsio/psio.h>
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

void IWL::fetch()
{
    psio_->read(itap_, IWL_KEY_BUF, (char *) &(lastbuf_), sizeof(int),
  	    bufpos_, &bufpos_);
    psio_->read(itap_, IWL_KEY_BUF, (char *) &(inbuf_), sizeof(int),
  	    bufpos_, &bufpos_);
    psio_->read(itap_, IWL_KEY_BUF, (char *) labels_, ints_per_buf_ * 
  	    4 * sizeof(Label), bufpos_, &bufpos_);
    psio_->read(itap_, IWL_KEY_BUF, (char *) values_, ints_per_buf_ *
  	    sizeof(Value), bufpos_, &bufpos_);
    idx_ = 0;
}

/*!
** iwl_buf_fetch()
**
** Fetch an IWL buffer from disk
** David Sherrill, 26 June 1996
** \ingroup IWL
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

}

