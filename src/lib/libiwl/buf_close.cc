/*! \defgroup IWL libiwl: I/O Library for Integrals with Labels */

/*!
  \file buf_close.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <stdlib.h>
#include <libpsio/psio.h>
#include "iwl.h"

extern "C" {

/*!
** IWL_BUF_CLOSE()
** 
**	\param Buf      Buffer to be closed
**	\param keep    Do not delete if keep==1
**
** Close a Integrals With Labels Buffer
** \ingroup (IWL)
*/
void iwl_buf_close(struct iwlbuf *Buf, int keep)
{

   psio_close(Buf->itap, keep ? 1 : 0);
   free(Buf->labels);
   free(Buf->values);
}

} /* extern "C" */