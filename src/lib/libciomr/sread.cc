#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"
#include <libciomr/libciomr.h>

extern "C" {

/* reads nlen bytes from itape into array starting at current pointer */
/* location */

void sread(int itape,char* array,int nlen)
   {
      PSI_FPTR start, end;

      start=ptr.wptr[itape];

      wreadw(itape,(char *) array,nlen,start,&end);
      ptr.wptr[itape]=((ptr.wptr[itape]-1+4096)/4096)*4096;
   }

} /* extern "C" */
