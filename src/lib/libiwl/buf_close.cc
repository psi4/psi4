/*! \defgroup IWL libiwl: I/O Library for Integrals with Labels */

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cstdlib>
#include <libpsio/psio.h>
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

IWL::~IWL()
{
    close();
}

void IWL::close()
{
    if (psio_->open_check(itap_))
        psio_->close(itap_, keep_);
    Communicator::world->barrier();
    if (labels_)
        delete[](labels_);
    if (values_)
        delete[](values_);
    labels_ = NULL;
    values_ = NULL;
}

/*!
** IWL_BUF_CLOSE()
** 
**	\param Buf      Buffer to be closed
**	\param keep    Do not delete if keep==1
**
** Close a Integrals With Labels Buffer
** \ingroup IWL
*/
void iwl_buf_close(struct iwlbuf *Buf, int k)
{

   psio_close(Buf->itap, k ? 1 : 0);
   free(Buf->labels);
   free(Buf->values);
}

}

