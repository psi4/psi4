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

/*! \defgroup IWL libiwl: I/O Library for Integrals with Labels */

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libpsio/psio.h"
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
void iwl_buf_close(struct iwlbuf *Buf, int keep)
{

   psio_close(Buf->itap, keep ? 1 : 0);
   free(Buf->labels);
   free(Buf->values);
}

}
