/*!
  \file buf_wrt_arr_SI_nocut.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include "iwl.h"

extern "C" {
	
/*!
** IWL_BUF_WRT_ARR_SI_nocut()
**
** This function writes out an array of two-electron
** integrals using the Integrals With Labels file format
** with indices stored in arrays of short int's. It DOES NOT
** use Buf->Cutoff when writing.
** Ed Valeev, February 1999
** \ingroup (IWL)
*/
void iwl_buf_wrt_arr_SI_nocut(struct iwlbuf *Buf, double *arr, short int *p, 
		     short int *q, short int *r, short int *s, int size)
{

  int i,idx;
  double value;
  Label *lblptr;
  Value *valptr;

  if (size < 0) {
    printf("(iwl_buf_wrt_arr_SI_nocut): Called with size = %d\n",
	   size);
    return;
  }
  
  if (Buf == NULL || arr == NULL || p == NULL || q == NULL || r == NULL 
      || s == NULL) {
    printf("(iwl_buf_wrt_arr_SI_nocut): Called with null pointer argument\n");
    return;
  }
  
  lblptr = Buf->labels;
  valptr = Buf->values;

  for (i=0; i<size; i++) {
    value = *arr++;
    idx = 4 * Buf->idx;
    lblptr[idx++] = (Label) p[i];
    lblptr[idx++] = (Label) q[i];
    lblptr[idx++] = (Label) r[i];
    lblptr[idx++] = (Label) s[i];
    valptr[Buf->idx] = (Value) value;
      
    Buf->idx++;

    if (Buf->idx == Buf->ints_per_buf) {
      Buf->lastbuf = 0;
      Buf->inbuf = Buf->idx;
      iwl_buf_put(Buf);
      Buf->idx = 0;
    }
  } /* end loop over i */
  
}

} /* extern "C" */