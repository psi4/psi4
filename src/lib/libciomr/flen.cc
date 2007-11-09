/*!
  \file flen.cc
  \ingroup (CIOMR)
*/
 
#include "includes.h"
#include "pointers.h"
#include "types.h"

extern "C" {

extern PSI_FPTR iosize_(int *);

/*!
** flen(): Get the number of bytes for file number 'itape'.
** \ingroup (CIOMR)
*/
PSI_FPTR flen(int itape)
   {
      PSI_FPTR length;
      length = iosize_(&itape);
      return(length);
   }

} /* extern "C" */
