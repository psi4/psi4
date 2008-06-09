/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <libpsio/psio.h>
#include "iwl.h"
#include "iwl.hpp"

namespace psi {
 
void IWL::put()
{
    psio_->write(itap_, IWL_KEY_BUF, (char *) &(lastbuf_), sizeof(int),
        bufpos_, &(bufpos_));
    psio_->write(itap_, IWL_KEY_BUF, (char *) &(inbuf_), sizeof(int),
        bufpos_, &(bufpos_));
    psio_->write(itap_, IWL_KEY_BUF, (char *) labels_, ints_per_buf_ * 
        4 * sizeof(Label), bufpos_, &(bufpos_));
    psio_->write(itap_, IWL_KEY_BUF, (char *) values_, ints_per_buf_ *
        sizeof(Value), bufpos_, &(bufpos_));
}  

/*!
** iwl_buf_put(struct iwlbuf *Buf)
**
** Put an IWL buffer to disk
** David Sherrill, 26 June 1996
** \ingroup IWL
*/
void iwl_buf_put(struct iwlbuf *Buf)
{
  psio_write(Buf->itap, IWL_KEY_BUF, (char *) &(Buf->lastbuf), sizeof(int),
	     Buf->bufpos, &(Buf->bufpos));
  psio_write(Buf->itap, IWL_KEY_BUF, (char *) &(Buf->inbuf), sizeof(int),
	     Buf->bufpos, &(Buf->bufpos));
  psio_write(Buf->itap, IWL_KEY_BUF, (char *) Buf->labels, Buf->ints_per_buf * 
	     4 * sizeof(Label), Buf->bufpos, &(Buf->bufpos));
  psio_write(Buf->itap, IWL_KEY_BUF, (char *) Buf->values, Buf->ints_per_buf *
	     sizeof(Value), Buf->bufpos, &(Buf->bufpos));
}

}

