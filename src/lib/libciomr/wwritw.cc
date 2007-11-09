/*!
** \file wwritw.cc
** \ingroup (CIOMR)
*/

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"
#include <libciomr/libciomr.h>

extern "C" {

/*!
** wwritw: writes nwords bytes to unit from buffer starting at fword
** and returns nxtwrd which is the current pointer location.
**
** \param unit   = file number
** \param buffer = buffer holding info to write
** \param nwords = number of bytes to write
** \param fword  = file pointer to where to put first byte written
** \param nxtwrd = file pointer to next byte on disk (returned value)
**
** \ingroup (CIMOR)
*/
void wwritw(int unit, char *buffer, int nwords, PSI_FPTR fword, 
            PSI_FPTR *nxtwrd)
{
  iowrr_(&unit,buffer,&fword,&nwords);
  *nxtwrd = fword + nwords;
  ptr.wptr[unit]= *nxtwrd;
}

} /* extern "C" */
