#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"
#include <libciomr/libciomr.h>

extern "C" {

/* writes nlen bytes from array to file itape starting at current */
/* pointer location */

void swrit(int itape,char* array,int nlen)
   {
      PSI_FPTR start, end;

      start=ptr.wptr[itape];

      wwritw(itape,(char *) array,nlen,start,&end);
      ptr.wptr[itape]=((ptr.wptr[itape]-1+4096)/4096)*4096;
   }

} /* extern "C" */
