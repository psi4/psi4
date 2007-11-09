/*!
** \file rclose.cc
** \ingroup (CIOMR)
*/

#include <libciomr/libciomr.h>

extern "C" {

/*!
** rclose: close a binary file
**
** \param unit = file number
** \param status = 3 to keep file, 4 to erase it
**
** \ingroup (CIOMR)
*/
void rclose(int unit, int status)
{
      ioclos_(&unit,&status);
   }

} /* extern "C" */
