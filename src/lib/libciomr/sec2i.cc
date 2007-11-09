#include "includes.h"
#include "pointers.h"

extern "C" {

/* converts sector pointer to integer pointer */
/* returns address of 1st element of sector n */

PSI_FPTR sec2i(int n)
   {
      PSI_FPTR num;
  
      num=(PSI_FPTR) (sizeof(int)*sector*n);
      return(num);
    }

} /* extern "C" */
