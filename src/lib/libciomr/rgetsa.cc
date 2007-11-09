#include <psifiles.h>
#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern "C" {

void rgetsa(int unit,int* iadr)
   {
      PSI_FPTR ipos,irec,test;
    
      ipos=ptr.wptr[unit];
      irec=(ipos/sizeof(PSI_FPTR)*sector)+1;
      test=sizeof(PSI_FPTR)*sector*(irec-1);

      if (ipos != test) {
         fprintf(stderr,"error encountered in rgetsa for file %d",unit);
         fprintf(stderr,"ipos,test = %lu %lu",ipos,test);
         exit(PSI_RETURN_FAILURE);
         }

      *iadr=irec;
  }

} /* extern "C" */
