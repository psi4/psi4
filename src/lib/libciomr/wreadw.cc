/*!
** \file wreadw.cc
** \ingroup (CIOMR)
*/

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"
#include <libciomr/libciomr.h>

extern "C" {

/*!
** wreadw: reads size bytes from tape into buffer starting at fword.
** nxtwrd is modified to give the new current pointer location after the
** read operation is completed.
**
**   \param tape   = file number
**   \param buffer = buffer to store read information
**   \param size   = number of bytes to be read
**   \param fword  = first byte of buffer to be read
**   \param nxtwrd = pointer to hold file pointer position after read
**
** \ingroup (CIOMR)
*/  
void wreadw(int tape, char *buffer, int size, PSI_FPTR fword, PSI_FPTR *nxtwrd)
   {
      iordr_(&tape,buffer,&fword,&size);
      *nxtwrd = fword + size;
      ptr.wptr[tape] = *nxtwrd;
   }

} /* extern "C" */
