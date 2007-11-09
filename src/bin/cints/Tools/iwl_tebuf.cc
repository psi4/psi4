/*! \file iwl_tebuf.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<libiwl/iwl.h>
#include<libint/libint.h>
#include"defines.h"
#include"data_structs.h"


namespace psi { namespace CINTS {

/*!
** IWL_BUF_WRT_STRUCT_nocut()
**
** This function writes out an array of two-electron
** integrals using the Integrals With Labels file format
** with integrals stored in the array of struct tebuf. It DOES NOT
** use Buf->Cutoff when writing.
** Ed Valeev, May 1999
**
*/
void iwl_buf_wrt_struct_nocut(struct iwlbuf *Buf, struct tebuf *Tebuf, int size)
{

  int i,j,idx;
  Label *lblptr;
  Value *valptr;

  lblptr = Buf->labels;
  valptr = Buf->values;

  for (i=0; i<size; i++) {
    idx = 4 * Buf->idx;
    lblptr[idx++] = (Label) Tebuf[i].i;
    lblptr[idx++] = (Label) Tebuf[i].j;
    lblptr[idx++] = (Label) Tebuf[i].k;
    lblptr[idx++] = (Label) Tebuf[i].l;
    valptr[Buf->idx] = (Value) Tebuf[i].val;
      
    Buf->idx++;

    if (Buf->idx == Buf->ints_per_buf) {
      Buf->lastbuf = 0;
      Buf->inbuf = Buf->idx;
      iwl_buf_put(Buf);
      Buf->idx = 0;
    }
  } /* end loop over i */

  return;
}


/*!
** IWL_BUF_WRT_STRUCT()
**
** This function writes out an array of two-electron
** integrals using the Integrals With Labels file format
** with integrals stored in the array of struct tebuf. It
** uses specified cutoff when writing.
** Ed Valeev, May 1999
**
*/
void iwl_buf_wrt_struct(struct iwlbuf *Buf, struct tebuf *Tebuf, int size, double cutoff)
{

  int i,j,idx;
  Label *lblptr;
  Value *valptr;

  lblptr = Buf->labels;
  valptr = Buf->values;

  for (i=0; i<size; i++) {
    if (fabs((double) Tebuf[i].val) > cutoff) {
      idx = 4 * Buf->idx;
      lblptr[idx++] = (Label) Tebuf[i].i;
      lblptr[idx++] = (Label) Tebuf[i].j;
      lblptr[idx++] = (Label) Tebuf[i].k;
      lblptr[idx++] = (Label) Tebuf[i].l;
      valptr[Buf->idx] = (Value) Tebuf[i].val;
      
      Buf->idx++;

      if (Buf->idx == Buf->ints_per_buf) {
	Buf->lastbuf = 0;
	Buf->inbuf = Buf->idx;
	iwl_buf_put(Buf);
	Buf->idx = 0;
      }
    }
  } /* end loop over i */
  
}

};};
