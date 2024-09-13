/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
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

IWL::~IWL() { close(); }

void IWL::close() {
    if (psio_->open_check(itap_)) psio_->close(itap_, keep_);
    if (labels_) delete[](labels_);
    if (values_) delete[](values_);
    labels_ = nullptr;
    values_ = nullptr;
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
void PSI_API iwl_buf_close(struct iwlbuf *Buf, int keep) {
    psio_close(Buf->itap, keep ? 1 : 0);
    free(Buf->labels);
    free(Buf->values);
}
}
