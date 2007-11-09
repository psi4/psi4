#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"
#include <libciomr/libciomr.h>

extern "C" {

/* writes nlen bytes from array to file itape starting at sector irec */

void rwrit(int itape,char* array,int nlen,int irec)
   {
    PSI_FPTR ipos, junk;

    ipos = sec2i(--irec);
    wwritw(itape,(char *) array,nlen,ipos,&junk);
   }

} /* extern "C" */
