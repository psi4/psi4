/*!
  \file buf_init.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <stdlib.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include "iwl.hpp"
#include "iwl.h"

extern "C" {
	
extern FILE *outfile;

};

using namespace psi;

IWL::IWL(PSIO* psio_obj, int itape, double cutoff, int oldfile, int readflag)
  : keep(true)
{
  psio = psio_obj;
  
  /*! set up buffer info */
  Buf.itap = itape;
  Buf.bufpos = PSIO_ZERO;
  Buf.ints_per_buf = IWL_INTS_PER_BUF;
  Buf.cutoff = cutoff;
  Buf.bufszc = 2 * sizeof(int) + Buf.ints_per_buf * 4 * sizeof(Label) +
    Buf.ints_per_buf * sizeof(Value);
  Buf.lastbuf = 0;
  Buf.inbuf = 0;
  Buf.idx = 0;

  /*! make room in the buffer */
  Buf.labels = (Label *) malloc (4 * Buf.ints_per_buf * sizeof(Label));
  Buf.values = (Value *) malloc (Buf.ints_per_buf * sizeof(Value));

  /*! open the output file */
  /*! Note that we assume that if oldfile isn't set, we O_CREAT the file */
  psio->open(Buf.itap, oldfile ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
  if (oldfile && (psio_tocscan(Buf.itap, IWL_KEY_BUF) == NULL)) {
    fprintf(outfile,"iwl_buf_init: Can't open file %d\n", Buf.itap);
    psio->close(Buf.itap,0);
    return;
  } 

  /*! go ahead and read a buffer */
  if (readflag) fetch();
}

extern "C" {
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
** \ingroup (IWL)
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
    fprintf(outfile,"iwl_buf_init: Can't open file %d\n", Buf->itap);
    psio_close(Buf->itap,0);
    return;
  } 

  /*! go ahead and read a buffer */
  if (readflag) iwl_buf_fetch(Buf);
}

} /* extern "C" */