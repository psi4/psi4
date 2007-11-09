#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"
#include <libciomr/libciomr.h>

extern "C" {

/* sets file pointer for unit to sector address */

void rsetsa(int unit,int address)
   {
      PSI_FPTR ipos;

      ipos = (PSI_FPTR) sec2i(--address);
      ptr.wptr[unit]=ipos;
    }

} /* extern "C" */
