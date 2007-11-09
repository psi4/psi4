#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"
#include <libciomr/libciomr.h>

extern "C" {

/* reads nlen bytes into array starting at irec */

void rread(int itape,char* array,int nlen,int irec)
   {
      PSI_FPTR ipos, junk;

      ipos = sec2i(--irec);
      wreadw(itape,(char *) array,nlen,ipos,&junk);
   }

} /* extern "C" */
