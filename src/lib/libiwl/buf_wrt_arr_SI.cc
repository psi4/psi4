/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

void IWL::write_array_SI(double *arr, short int *p, short int *q, short int *r, 
    short int *s, int size)
{
    int i,idx;
    double value;
    Label *lblptr;
    Value *valptr;

    if (size < 0) {
        printf("(iwl_buf_wrt_arr_SI): Called with size = %d\n",
            size);
        return;
    }

    if (arr == NULL || p == NULL || q == NULL || r == NULL 
    || s == NULL) {
        printf("(iwl_buf_wrt_arr_SI): Called with null pointer argument\n");
        return;
    }

    lblptr = labels_;
    valptr = values_;

    for (i=0; i<size; i++) {
        value = *arr++;
        if (fabs(value) > cutoff_) {
            idx = 4 * idx_;
            lblptr[idx++] = (Label) p[i];
            lblptr[idx++] = (Label) q[i];
            lblptr[idx++] = (Label) r[i];
            lblptr[idx++] = (Label) s[i];
            valptr[idx_] = (Value) value;

            idx_++;

            if (idx_ == ints_per_buf_) {
                lastbuf_ = 0;
                inbuf_ = idx_;
                put();
                idx_ = 0;
            }
        } /* end if cutoff */
    } /* end loop over i */
}

/*!
** IWL_BUF_WRT_ARR_SI()
**
** This function writes out an array of two-electron
** integrals using the Integrals With Labels file format
** with indices stored in arrays of short int's
** Ed Valeev, February 1999
** \ingroup IWL
*/
void iwl_buf_wrt_arr_SI(struct iwlbuf *Buf, double *arr, short int *p, 
		     short int *q, short int *r, short int *s, int size)
{

  int i,idx;
  double value;
  Label *lblptr;
  Value *valptr;

  if (size < 0) {
    printf("(iwl_buf_wrt_arr_SI): Called with size = %d\n",
	   size);
    return;
  }
  
  if (Buf == NULL || arr == NULL || p == NULL || q == NULL || r == NULL 
      || s == NULL) {
    printf("(iwl_buf_wrt_arr_SI): Called with null pointer argument\n");
    return;
  }
  
  lblptr = Buf->labels;
  valptr = Buf->values;

  for (i=0; i<size; i++) {
    value = *arr++;
    if (fabs(value) > Buf->cutoff) {
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
      
    } /* end if cutoff */
  } /* end loop over i */
  
}

}

