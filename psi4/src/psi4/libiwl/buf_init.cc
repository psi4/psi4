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

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libpsio/psio.h"
#include "iwl.h"
#include "iwl.hpp"
#include "psi4/psi4-dec.h" //need outfile

namespace psi {

IWL::IWL()
{
    /*! set up buffer info */
    itap_ = -1;
    bufpos_ = PSIO_ZERO;
    ints_per_buf_ = IWL_INTS_PER_BUF;
    cutoff_ = 1.e-14;
    bufszc_ = 2 * sizeof(int) + ints_per_buf_ * 4 * sizeof(Label) +
        ints_per_buf_ * sizeof(Value);
    lastbuf_ = 0;
    inbuf_ = 0;
    idx_ = 0;
}

IWL::IWL(PSIO *psio, int it, double coff, int oldfile, int readflag):
    keep_(true)
{
    init(psio, it, coff, oldfile, readflag);
}

void IWL::init(PSIO *psio, int it, double coff, int oldfile, int readflag)
{
    psio_ = psio;

    /*! set up buffer info */
    itap_ = it;
    bufpos_ = PSIO_ZERO;
    ints_per_buf_ = IWL_INTS_PER_BUF;
    cutoff_ = coff;
    bufszc_ = 2 * sizeof(int) + ints_per_buf_ * 4 * sizeof(Label) +
        ints_per_buf_ * sizeof(Value);
    lastbuf_ = 0;
    inbuf_ = 0;
    idx_ = 0;

    /*! make room in the buffer */
    // labels_ = (Label *) malloc (4 * ints_per_buf_ * sizeof(Label));
    // values_ = (Value *) malloc (ints_per_buf_ * sizeof(Value));
    labels_ = new Label[4 * ints_per_buf_];
    values_ = new Value[ints_per_buf_];

    /*! open the output file */
    /*! Note that we assume that if oldfile isn't set, we O_CREAT the file */
    psio_->open(itap_, oldfile ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
    if (oldfile && (psio_->tocscan(itap_, IWL_KEY_BUF) == NULL)) {
        outfile->Printf("iwl_buf_init: Can't open file %d\n", itap_);
        psio_->close(itap_,0);
        return;
    }

    /*! go ahead and read a buffer */
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
void iwl_buf_init(struct iwlbuf *Buf, int itape, double cutoff,
      int oldfile, int readflag)
{

  /*! set up buffer info */
  Buf->itap = itape;
  Buf->bufpos = PSIO_ZERO;
  Buf->ints_per_buf = IWL_INTS_PER_BUF;
  Buf->cutoff = cutoff;
  Buf->bufszc = 2 * sizeof(int) + Buf->ints_per_buf * 4 * sizeof(Label) +
    Buf->ints_per_buf * sizeof(Value);
  Buf->lastbuf = 0;
  Buf->inbuf = 0;
  Buf->idx = 0;

  /*! make room in the buffer */
  Buf->labels = (Label *) malloc (4 * Buf->ints_per_buf * sizeof(Label));
  Buf->values = (Value *) malloc (Buf->ints_per_buf * sizeof(Value));

  /*! open the output file */
  /*! Note that we assume that if oldfile isn't set, we O_CREAT the file */
  psio_open(Buf->itap, oldfile ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
  if (oldfile && (psio_tocscan(Buf->itap, IWL_KEY_BUF) == NULL)) {
    outfile->Printf("iwl_buf_init: Can't open file %d\n", Buf->itap);
    psio_close(Buf->itap,0);
    return;
  }

  /*! go ahead and read a buffer */
  if (readflag) iwl_buf_fetch(Buf);

}

}
