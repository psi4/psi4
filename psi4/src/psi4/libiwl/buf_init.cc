/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "iwl.h"
#include "iwl.hpp"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {

namespace {

// Constant on-disk bucket size in bytes. Matches the layout written by
// buf_put.cc / read by buf_fetch.cc.
inline int bucket_bytes(int ints_per_buf) {
    return 2 * static_cast<int>(sizeof(int)) + ints_per_buf * 4 * static_cast<int>(sizeof(Label)) +
           ints_per_buf * static_cast<int>(sizeof(Value));
}

// Open the underlying psio file for an IWL buffer and verify it has the
// expected TOC entry when reopening an existing file. Throws on failure
// instead of leaving the buffer in a half-initialized state.
void open_iwl_unit(PSIO *psio, int itap, bool oldfile) {
    psio->open(itap, oldfile ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
    if (oldfile && psio->tocscan(itap, IWL_KEY_BUF) == nullptr) {
        psio->close(itap, 0);
        throw PSIEXCEPTION("iwl_buf_init: file " + std::to_string(itap) +
                           " does not contain an IWL buffer (TOC key '" + IWL_KEY_BUF + "' missing)");
    }
}

}  // namespace

IWL::IWL()
    : itap_(-1),
      bufpos_(PSIO_ZERO),
      ints_per_buf_(IWL_INTS_PER_BUF),
      bufszc_(bucket_bytes(IWL_INTS_PER_BUF)),
      cutoff_(1.e-14),
      lastbuf_(0),
      inbuf_(0),
      idx_(0),
      labels_(nullptr),
      values_(nullptr),
      psio_(nullptr),
      keep_(true) {}

IWL::IWL(PSIO *psio, int it, double coff, int oldfile, int readflag) : IWL() { init(psio, it, coff, oldfile, readflag); }

void IWL::init(PSIO *psio, int it, double coff, int oldfile, int readflag) {
    psio_ = psio;
    itap_ = it;
    bufpos_ = PSIO_ZERO;
    ints_per_buf_ = IWL_INTS_PER_BUF;
    cutoff_ = coff;
    bufszc_ = bucket_bytes(ints_per_buf_);
    lastbuf_ = 0;
    inbuf_ = 0;
    idx_ = 0;

    labels_ = new Label[4 * ints_per_buf_];
    values_ = new Value[ints_per_buf_];

    open_iwl_unit(psio_, itap_, oldfile != 0);

    if (readflag) fetch();
}

/*!
** iwl_buf_init()
**
**	\param Buf               Buffer to be initialised
**	\param itape		Filenumber
**	\param cutoff           Cutoff for keeping integral
**	\param oldfile		If ==0 create file
**	\param readflag		If ==1 fetch buffer
**
** Prepare a PSI Buffer according to the Integrals
** With Labels format for reading or writing.  Important to set
** readflag=1 if opening for reading, since other IWL buffer read
** routines anticipate that there is already data in the buffer.
**
** David Sherrill, March 1995
** Revised 6/26/96 by CDS for new format
** \ingroup IWL
*/
void PSI_API iwl_buf_init(struct iwlbuf *Buf, int itape, double cutoff, int oldfile, int readflag) {
    Buf->itap = itape;
    Buf->bufpos = PSIO_ZERO;
    Buf->ints_per_buf = IWL_INTS_PER_BUF;
    Buf->cutoff = cutoff;
    Buf->bufszc = bucket_bytes(Buf->ints_per_buf);
    Buf->lastbuf = 0;
    Buf->inbuf = 0;
    Buf->idx = 0;

    Buf->labels = new Label[4 * Buf->ints_per_buf];
    Buf->values = new Value[Buf->ints_per_buf];

    open_iwl_unit(_default_psio_lib_.get(), Buf->itap, oldfile != 0);

    if (readflag) iwl_buf_fetch(Buf);
}
}
